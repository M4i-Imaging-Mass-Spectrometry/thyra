"""Tests for thyra.convert streaming auto-detection.

The standard SpatialData converter pre-allocates the whole dataset's
sparse matrix in memory, so ``convert_msi(streaming="auto")`` must switch
to the streaming converter once the estimated footprint would claim too
much RAM. These tests pin that decision and the RAM-aware budget.
"""

from types import SimpleNamespace

import pytest

from thyra import convert
from thyra.convert import _should_use_streaming, _streaming_memory_budget_gb


def _reader_with_pixels(n_pixels: int) -> SimpleNamespace:
    """A stub reader exposing get_essential_metadata().dimensions."""
    meta = SimpleNamespace(dimensions=(n_pixels, 1, 1))
    return SimpleNamespace(get_essential_metadata=lambda: meta)


class TestShouldUseStreaming:
    def test_explicit_true_always_streams(self):
        assert _should_use_streaming(True, _reader_with_pixels(1)) is True

    def test_explicit_false_never_streams(self):
        assert _should_use_streaming(False, _reader_with_pixels(10**9)) is False

    def test_auto_small_dataset_stays_in_memory(self, monkeypatch):
        monkeypatch.setattr(convert, "_streaming_memory_budget_gb", lambda: 4.0)
        # 10k pixels -> ~0.07 GB estimated, below the 4 GB budget.
        assert _should_use_streaming("auto", _reader_with_pixels(10_000)) is False

    def test_auto_large_dataset_streams(self, monkeypatch):
        monkeypatch.setattr(convert, "_streaming_memory_budget_gb", lambda: 4.0)
        # 100k pixels -> ~7.5 GB estimated, above the 4 GB budget.
        assert _should_use_streaming("auto", _reader_with_pixels(100_000)) is True

    def test_auto_estimation_error_defaults_to_in_memory(self, monkeypatch):
        monkeypatch.setattr(convert, "_streaming_memory_budget_gb", lambda: 4.0)

        def boom():
            raise RuntimeError("no metadata")

        reader = SimpleNamespace(get_essential_metadata=boom)
        # A metadata failure must not raise out of the decision helper.
        assert _should_use_streaming("auto", reader) is False


class TestStreamingMemoryBudget:
    def test_env_override_caps_budget(self, monkeypatch):
        import psutil

        monkeypatch.setenv("THYRA_STREAMING_MAX_GB", "2")
        monkeypatch.setattr(
            psutil,
            "virtual_memory",
            lambda: SimpleNamespace(available=999 * 1024**3),
        )
        # Huge available RAM -> the 2 GB env cap wins the min().
        assert _streaming_memory_budget_gb() == pytest.approx(2.0)

    def test_ram_fraction_can_win(self, monkeypatch):
        import psutil

        monkeypatch.delenv("THYRA_STREAMING_MAX_GB", raising=False)
        monkeypatch.setattr(
            psutil,
            "virtual_memory",
            lambda: SimpleNamespace(available=4 * 1024**3),
        )
        # 0.5 * 4 GB = 2 GB < 4 GB cap -> the RAM fraction wins.
        assert _streaming_memory_budget_gb() == pytest.approx(2.0)

    def test_invalid_env_falls_back_to_cap(self, monkeypatch):
        import psutil

        monkeypatch.setenv("THYRA_STREAMING_MAX_GB", "not-a-number")
        monkeypatch.setattr(
            psutil,
            "virtual_memory",
            lambda: SimpleNamespace(available=999 * 1024**3),
        )
        # Bad override is ignored; the 4 GB default cap applies.
        assert _streaming_memory_budget_gb() == pytest.approx(4.0)
