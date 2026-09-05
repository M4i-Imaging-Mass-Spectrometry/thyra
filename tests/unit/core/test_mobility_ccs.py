"""Collision cross section from 1/K0: the pure-Python Mason-Schamp route.

The SDK-backed route is pinned to this one in
``tests/integration/test_bruker_tdf_synthetic.py``, where the vendor
library is loaded. Here: the prefactor is the published one, the
formula behaves as physics says it should, and a charge is never
invented.
"""

import numpy as np
import pytest

from thyra.core.mobility import (
    CCS_DRIFT_GAS_MASS_DA,
    CCS_TEMPERATURE_K,
    MASON_SCHAMP_PREFACTOR,
    MobilityAxis,
    ccs_from_one_over_k0,
    mason_schamp_ccs,
)

# The constant the vendor's documentation and the open reimplementations
# carry, in A^2 sqrt(Da K) cm^2 V^-1 s^-1.
PUBLISHED_PREFACTOR = 18509.8632163405


class TestMasonSchamp:
    def test_prefactor_is_the_published_constant(self):
        # Built from CODATA 2018 constants; the published value predates
        # them, hence the 2e-7 rather than exact equality.
        assert MASON_SCHAMP_PREFACTOR == pytest.approx(PUBLISHED_PREFACTOR, rel=1e-6)
        assert CCS_DRIFT_GAS_MASS_DA == 28.0134
        assert CCS_TEMPERATURE_K == 305.0

    def test_matches_the_closed_form(self):
        ook0 = np.array([0.8, 1.0, 1.2, 1.5, 1.9])
        mz = np.array([300.0, 500.0, 760.5, 1000.0, 1500.0])
        reduced = mz * 28.0134 / (mz + 28.0134)
        expected = PUBLISHED_PREFACTOR * ook0 / np.sqrt(reduced * 305.0)
        np.testing.assert_allclose(mason_schamp_ccs(ook0, mz, 1), expected, rtol=1e-6)

    def test_linear_in_one_over_k0(self):
        a = mason_schamp_ccs(1.0, 500.0, 1)
        b = mason_schamp_ccs(2.0, 500.0, 1)
        assert b == pytest.approx(2 * a)

    def test_charge_enters_through_the_ion_mass_too(self):
        # Doubling the charge doubles the prefactor but also doubles the
        # ion mass, so the result is a little under twice as large.
        one = mason_schamp_ccs(1.0, 500.0, 1)
        two = mason_schamp_ccs(1.0, 500.0, 2)
        assert 1.9 * one < two < 2.0 * one

    def test_broadcasts(self):
        out = mason_schamp_ccs(np.array([1.0, 1.1, 1.2]), 500.0, 1)
        assert out.shape == (3,)
        assert np.all(np.diff(out) > 0)

    def test_scalar_inputs_give_a_scalar_array(self):
        out = mason_schamp_ccs(1.0, 500.0, 1)
        assert np.ndim(out) == 0
        assert 150.0 < float(out) < 250.0

    @pytest.mark.parametrize("charge", [0, -1])
    def test_charge_must_be_a_positive_state(self, charge):
        with pytest.raises(ValueError, match="charge"):
            mason_schamp_ccs(1.0, 500.0, charge)

    @pytest.mark.parametrize("charge", [1.0, "1", True, None])
    def test_charge_must_be_an_integer(self, charge):
        with pytest.raises(TypeError, match="charge"):
            mason_schamp_ccs(1.0, 500.0, charge)


class TestDispatch:
    def test_without_an_sdk_it_is_the_formula(self):
        ook0 = np.array([1.0, 1.5])
        mz = np.array([400.0, 900.0])
        np.testing.assert_array_equal(
            ccs_from_one_over_k0(ook0, mz, 1), mason_schamp_ccs(ook0, mz, 1)
        )

    def test_an_sdk_without_the_export_falls_back(self):
        class NoExport:
            file_type = "tdf"
            _bound_conversions = {"tims_oneoverk0_to_ccs_for_mz": False}

        np.testing.assert_array_equal(
            ccs_from_one_over_k0(1.0, 500.0, 1, sdk=NoExport()),
            mason_schamp_ccs(1.0, 500.0, 1),
        )

    def test_a_tsf_sdk_is_ignored(self):
        class Tsf:
            file_type = "tsf"

            def oneoverk0_to_ccs(self, *args):  # pragma: no cover
                raise AssertionError("must not be called for TSF")

        np.testing.assert_array_equal(
            ccs_from_one_over_k0(1.0, 500.0, 1, sdk=Tsf()),
            mason_schamp_ccs(1.0, 500.0, 1),
        )

    def test_a_tdf_sdk_with_the_export_is_used(self):
        calls = []

        class WithExport:
            file_type = "tdf"
            _bound_conversions = {"tims_oneoverk0_to_ccs_for_mz": True}

            def oneoverk0_to_ccs(self, ook0, charge, mz):
                calls.append((ook0, charge, mz))
                return np.full(np.shape(ook0), 123.0)

        out = ccs_from_one_over_k0(
            np.array([1.0, 1.1]), np.array([5.0, 6.0]), 2, sdk=WithExport()
        )
        np.testing.assert_array_equal(out, [123.0, 123.0])
        assert calls[0][1] == 2

    def test_charge_is_checked_before_the_sdk_is_touched(self):
        class Explodes:
            file_type = "tdf"
            _bound_conversions = {"tims_oneoverk0_to_ccs_for_mz": True}

            def oneoverk0_to_ccs(self, *args):  # pragma: no cover
                raise AssertionError("must not be reached")

        with pytest.raises(ValueError, match="charge"):
            ccs_from_one_over_k0(1.0, 500.0, 0, sdk=Explodes())


class TestAxisUnsContract:
    def test_uns_carries_the_contract_keys(self):
        axis = MobilityAxis(
            kind_accession="MS:1002815",
            kind_name="inverse reduced ion mobility",
            unit_accession="MS:1002814",
            unit_name="volt-second per square centimeter",
            values=np.array([1.3, 1.2, 1.1]),
            acq_range=(1.0, 1.29),
            calibration={"model_type": 2, "coefficients": [1.0, 2.0]},
            source="bruker_tdf",
        )
        block = axis.to_uns()
        assert block["present"] is True
        assert block["type_accession"] == "MS:1002815"
        assert block["type_name"] == "inverse reduced ion mobility"
        assert block["unit_accession"] == "MS:1002814"
        assert block["unit_name"] == "volt-second per square centimeter"
        assert block["n_scans"] == 3
        assert block["values"].dtype == np.float64
        np.testing.assert_array_equal(block["acq_range"], [1.0, 1.29])
        assert block["acq_range"].dtype == np.float64
        assert block["calibration"] == {"model_type": 2, "coefficients": [1.0, 2.0]}
        assert block["source"] == "bruker_tdf"
        # The older spellings stay for stores written before the contract.
        assert block["n_values"] == 3
        assert (block["range_lower"], block["range_upper"]) == (1.0, 1.29)
        assert not any(":" in key for key in block)

    def test_range_falls_back_to_the_values(self):
        axis = MobilityAxis(
            kind_accession="MS:1002815",
            kind_name="inverse reduced ion mobility",
            unit_accession=None,
            unit_name=None,
            values=np.array([0.95, 1.35, 1.10]),
        )
        block = axis.to_uns()
        np.testing.assert_array_equal(block["acq_range"], [0.95, 1.35])
        assert "calibration" not in block and "source" not in block
