"""Converter plugins for MSI data output formats."""

import logging

logger = logging.getLogger(__name__)

try:
    from . import spatialdata  # noqa: F401

    logger.debug("Successfully imported spatialdata package")
except (ImportError, NotImplementedError) as e:
    # Skip if spatialdata dependencies not available or incompatible
    logger.error(f"SpatialData converter not available due to dependency issues: {e}")
    import traceback

    logger.error(f"Full traceback: {traceback.format_exc()}")
except Exception as e:
    logger.error(f"Unexpected error importing spatialdata package: {e}")
    import traceback

    logger.error(f"Full traceback: {traceback.format_exc()}")
