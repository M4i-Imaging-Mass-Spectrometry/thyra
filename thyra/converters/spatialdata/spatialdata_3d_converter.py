# thyra/converters/spatialdata/spatialdata_3d_converter.py

import logging
from typing import Any, Dict, Tuple

import numpy as np
from numpy.typing import NDArray

from .base_spatialdata_converter import SPATIALDATA_AVAILABLE, BaseSpatialDataConverter

if SPATIALDATA_AVAILABLE:
    import xarray as xr
    from anndata import AnnData
    from spatialdata.models import Image2DModel, TableModel
    from spatialdata.transformations import Identity


class SpatialData3DConverter(BaseSpatialDataConverter):
    """Converter for MSI data to SpatialData format as true 3D volume or single 2D slice."""

    def __init__(self, *args, **kwargs):
        """Initialize 3D converter with handle_3d=True."""
        kwargs["handle_3d"] = True  # Force 3D mode
        super().__init__(*args, **kwargs)

    def _create_data_structures(self) -> Dict[str, Any]:
        """Create data structures for 3D volume format.

        Returns:
            Dict containing tables, shapes, images, and data arrays for
            3D volume
        """
        # Return dictionaries to store tables, shapes, and images
        tables: Dict[str, Any] = {}
        shapes: Dict[str, Any] = {}
        images: Dict[str, Any] = {}

        if self._dimensions is None:
            raise ValueError("Dimensions are not initialized")
        if self._common_mass_axis is None:
            raise ValueError("Common mass axis is not initialized")

        n_x, n_y, n_z = self._dimensions

        return {
            "mode": "3d_volume",
            "sparse_matrix": self._create_sparse_matrix(),
            "coords_df": self._create_coordinates_dataframe(),
            "var_df": self._create_mass_dataframe(),
            "tables": tables,
            "shapes": shapes,
            "images": images,
            "tic_values": np.zeros((n_y, n_x, n_z), dtype=np.float64),
            "total_intensity": np.zeros(len(self._common_mass_axis), dtype=np.float64),
            "pixel_count": 0,
        }

    def _process_single_spectrum(
        self,
        data_structures: Dict[str, Any],
        coords: Tuple[int, int, int],
        mzs: NDArray[np.float64],
        intensities: NDArray[np.float64],
    ) -> None:
        """Process a single spectrum for 3D volume format.

        Delegates to parent's resampling-aware processing.

        Args:
            data_structures: Data structures for storing processed data
            coords: (x, y, z) pixel coordinates
            mzs: Array of m/z values
            intensities: Array of intensity values
        """
        # Delegate to parent's resampling-aware processing
        super()._process_single_spectrum(data_structures, coords, mzs, intensities)

    def _process_resampled_spectrum(
        self,
        data_structures: Dict[str, Any],
        coords: Tuple[int, int, int],
        mz_indices: NDArray[np.int_],
        intensities: NDArray[np.float64],
    ) -> None:
        """Process a spectrum with resampled intensities for 3D volume format.

        Args:
            data_structures: Data structures for storing processed data
            coords: (x, y, z) pixel coordinates
            mz_indices: Indices in the common mass axis (all indices for
            resampled)
            intensities: Resampled intensity values
        """
        if self._dimensions is None:
            raise ValueError("Dimensions are not initialized")

        x, y, z = coords

        # Calculate TIC for this pixel
        tic_value = float(np.sum(intensities))

        # Update total intensity for average spectrum calculation
        # Handle both resampled and non-resampled cases
        if len(intensities) == len(data_structures["total_intensity"]):
            # Resampled case - intensities match common mass axis length
            data_structures["total_intensity"] += intensities
        else:
            # Non-resampled case - need to map to indices
            for i, intensity in enumerate(intensities):
                if i < len(mz_indices) and mz_indices[i] < len(
                    data_structures["total_intensity"]
                ):
                    data_structures["total_intensity"][mz_indices[i]] += intensity
        data_structures["pixel_count"] += 1

        # Get pixel index for 3D volume
        pixel_idx = self._get_pixel_index(x, y, z)

        # Store TIC value for this pixel
        data_structures["tic_values"][y, x, z] = tic_value

        # Add to sparse matrix
        self._add_to_sparse_matrix(
            data_structures["sparse_matrix"],
            pixel_idx,
            mz_indices,
            intensities,
        )

        self._non_empty_pixel_count += 1

    def _finalize_data(self, data_structures: Dict[str, Any]) -> None:
        """Finalize 3D volume data by creating tables, shapes, and images.

        Args:
            data_structures: Data structures containing processed data
        """
        if not SPATIALDATA_AVAILABLE:
            raise ImportError("SpatialData dependencies not available")

        try:
            # Store pixel count for metadata
            self._non_empty_pixel_count = data_structures["pixel_count"]

            # Convert COO arrays to sparse matrix (CSC or CSR based on config)
            format_name = "CSC" if self._sparse_format == "csc" else "CSR"
            logging.info(f"Converting COO arrays to {format_name} format...")
            coo_arrays = data_structures["sparse_matrix"]
            current_idx = coo_arrays["current_idx"]

            # Trim arrays to actual size
            from scipy import sparse

            coo = sparse.coo_matrix(
                (
                    coo_arrays["data"][:current_idx],
                    (
                        coo_arrays["rows"][:current_idx],
                        coo_arrays["cols"][:current_idx],
                    ),
                ),
                shape=(coo_arrays["n_rows"], coo_arrays["n_cols"]),
                dtype=np.float64,
            )
            # Convert to configured sparse format
            if self._sparse_format == "csc":
                sparse_matrix = coo.tocsc()
            else:
                sparse_matrix = coo.tocsr()

            logging.info(
                f"Converted sparse matrix: {sparse_matrix.nnz:,} non-zero entries ({format_name})"
            )

            # Create AnnData
            adata = AnnData(
                X=sparse_matrix,
                obs=data_structures["coords_df"],
                var=data_structures["var_df"],
            )

            # Add average spectrum to .uns (use total_intensity to match
            # original behavior)
            adata.uns["average_spectrum"] = data_structures["total_intensity"]

            # Make sure region column exists and is correct
            region_key = f"{self.dataset_id}_pixels"
            if "region" not in adata.obs.columns:
                adata.obs["region"] = region_key

            # Ensure instance_key is a string column
            adata.obs["instance_key"] = adata.obs.index.astype(str)

            # Create table model
            table = TableModel.parse(
                adata,
                region=region_key,
                region_key="region",
                instance_key="instance_key",
            )

            # Add to tables and create shapes
            data_structures["tables"][self.dataset_id] = table
            data_structures["shapes"][region_key] = self._create_pixel_shapes(
                adata, is_3d=True
            )

            # Create TIC image
            self._create_tic_image(data_structures)

        except Exception as e:
            logging.error(f"Error processing 3D volume: {e}")
            import traceback

            logging.debug(f"Detailed traceback:\n{traceback.format_exc()}")
            raise

    def _create_tic_image(self, data_structures: Dict[str, Any]) -> None:
        """Create TIC image for 3D volume or 2D slice.

        Args:
            data_structures: Data structures containing TIC values
        """
        if self._dimensions is None:
            raise ValueError("Dimensions are not initialized")

        n_x, n_y, n_z = self._dimensions

        if n_z > 1:
            # True 3D TIC image
            tic_values = data_structures["tic_values"]
            z_size, y_size, x_size = tic_values.shape

            # Add channel dimension for 3D image
            tic_values_with_channel = tic_values.reshape(1, z_size, y_size, x_size)

            tic_image = xr.DataArray(
                tic_values_with_channel,
                dims=("c", "z", "y", "x"),
                coords={
                    "c": [0],  # Single channel
                    "z": np.arange(z_size) * self.pixel_size_um,
                    "y": np.arange(y_size) * self.pixel_size_um,
                    "x": np.arange(x_size) * self.pixel_size_um,
                },
            )

            # Create Image model for 3D image
            transform = Identity()
            try:
                from spatialdata.models import Image3DModel

                data_structures["images"][f"{self.dataset_id}_tic"] = (
                    Image3DModel.parse(
                        tic_image,
                        transformations={
                            self.dataset_id: transform,
                            "global": transform,
                        },
                    )
                )
            except (ImportError, AttributeError):
                # Fallback if Image3DModel is not available
                logging.warning("Image3DModel not available, using generic image model")
                from spatialdata.models import ImageModel

                data_structures["images"][f"{self.dataset_id}_tic"] = ImageModel.parse(
                    tic_image,
                    transformations={
                        self.dataset_id: transform,
                        "global": transform,
                    },
                )
        else:
            # Single 2D slice
            tic_values = data_structures["tic_values"]

            # Handle both 3D array with single z-slice and 2D array
            if len(tic_values.shape) == 3:
                tic_values = tic_values[:, :, 0]

            y_size, x_size = tic_values.shape

            # Add channel dimension to make it (c, y, x)
            tic_values_with_channel = tic_values.reshape(1, y_size, x_size)

            tic_image = xr.DataArray(
                tic_values_with_channel,
                dims=("c", "y", "x"),
                coords={
                    "c": [0],  # Single channel
                    "y": np.arange(y_size) * self.pixel_size_um,
                    "x": np.arange(x_size) * self.pixel_size_um,
                },
            )

            # Create Image2DModel for the TIC image
            transform = Identity()
            data_structures["images"][f"{self.dataset_id}_tic"] = Image2DModel.parse(
                tic_image,
                transformations={
                    self.dataset_id: transform,
                    "global": transform,
                },
            )

        # Add optical images if available
        self._add_optical_images(data_structures)
