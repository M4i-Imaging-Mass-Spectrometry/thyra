"""The mass-axis guard's limits: the count ceiling, tables, regions and the process cap.

Issue #528 found three ways the guard admitted an axis that then failed:

* past 2^28 bins the count array of the first pass raised a bare
  ``MemoryError`` with a traceback, however much memory was free (F59);
* each region and each extra plane table added memory the projection did
  not count (F148, F116);
* it read the machine's free memory, so under a job or cgroup limit it
  admitted an axis the limit then failed (F194).

Nothing here allocates anything large: the limits are read from fake
cgroup trees and, on Windows, from a real job object in a child process.
"""

from __future__ import annotations

import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from thyra.resampling import memory_guard
from thyra.resampling.memory_guard import (
    AXIS_BYTES_PER_BIN,
    AXIS_REFUSE_FRACTION,
    MAX_COUNT_KEYS,
    REGION_BYTES_PER_BIN,
    TABLE_BYTES_PER_BIN,
    axis_bytes_per_bin,
    mass_axis_refusal,
)

GIB = 1024**3


class TestTheCountCeiling:
    def test_it_is_two_to_the_28(self):
        assert MAX_COUNT_KEYS == 2**28

    def test_it_matches_what_the_count_array_refuses(self):
        """One number, so the guard cannot admit what the engine then refuses."""
        from thyra.converters.spatialdata.csc_assembly import count_refusal

        assert count_refusal(MAX_COUNT_KEYS) is None
        assert count_refusal(MAX_COUNT_KEYS + 1) is not None

    def test_the_last_countable_axis_passes_on_a_big_machine(self):
        assert mass_axis_refusal(MAX_COUNT_KEYS, available_gb=1e6) is None

    def test_one_more_bin_is_refused_however_much_is_free(self):
        refusal = mass_axis_refusal(MAX_COUNT_KEYS + 1, available_gb=1e6)

        assert refusal is not None
        assert "268,435,457 bins" in refusal
        assert "268,435,456" in refusal
        assert "--resample-bins" in refusal

    def test_a_caller_with_other_levers_names_its_own(self):
        """``plan_mass_axis`` has no CLI options to name (issue #473)."""
        refusal = mass_axis_refusal(
            MAX_COUNT_KEYS + 1, available_gb=1e6, levers="plan a coarser axis"
        )

        assert refusal is not None
        assert refusal.endswith("plan a coarser axis")
        assert "--resample-bins" not in refusal


class TestTablesAndRegions:
    def test_one_table_and_one_region_is_the_measured_figure(self):
        assert axis_bytes_per_bin() == AXIS_BYTES_PER_BIN
        # One region builds no per-region arrays (D57).
        assert axis_bytes_per_bin(1, 1) == AXIS_BYTES_PER_BIN

    def test_each_table_and_each_region_adds_its_own(self):
        assert axis_bytes_per_bin(3, 6) == (
            AXIS_BYTES_PER_BIN + 2 * TABLE_BYTES_PER_BIN + 6 * REGION_BYTES_PER_BIN
        )

    def test_the_boundary_moves_with_the_tables(self):
        free = 8.0
        per_bin = axis_bytes_per_bin(10, 0)
        bins = int(free * AXIS_REFUSE_FRACTION * GIB / per_bin)

        assert mass_axis_refusal(bins, available_gb=free, n_tables=10) is None
        assert mass_axis_refusal(bins + 1, available_gb=free, n_tables=10)
        # The same axis on one table is far inside the budget.
        assert mass_axis_refusal(bins + 1, available_gb=free) is None

    def test_the_refusal_says_how_the_figure_is_made_up(self):
        refusal = mass_axis_refusal(
            18_000_000, available_gb=2.0, n_tables=3, n_regions=6
        )

        assert refusal is not None
        assert f"{axis_bytes_per_bin(3, 6)} bytes per bin" in refusal
        assert f"{TABLE_BYTES_PER_BIN} for each of 2 more tables" in refusal
        assert f"{REGION_BYTES_PER_BIN} for each of 6 regions" in refusal

    def test_the_warning_says_it_too(self, thyra_logs):
        free = 8.0
        bins = int(free * 0.3 * GIB / axis_bytes_per_bin(1, 4))

        with thyra_logs("thyra.resampling.memory_guard") as records:
            assert mass_axis_refusal(bins, available_gb=free, n_regions=4) is None

        assert any("for each of 4 regions" in r.getMessage() for r in records)


def _cgroup_v2(tmp_path: Path, levels: dict) -> Path:
    """A fake cgroup v2 tree; ``levels`` maps a path to (max, current, inactive_file)."""
    root = tmp_path / "sys_fs_cgroup"
    for path, (limit, current, cache) in levels.items():
        folder = root.joinpath(*[p for p in path.strip("/").split("/") if p])
        folder.mkdir(parents=True, exist_ok=True)
        (folder / "memory.max").write_text(f"{limit}\n")
        (folder / "memory.current").write_text(f"{current}\n")
        (folder / "memory.stat").write_text(f"anon 1\ninactive_file {cache}\n")
    return root


def _proc(tmp_path: Path, text: str) -> Path:
    path = tmp_path / "proc_self_cgroup"
    path.write_text(text)
    return path


class TestCgroupLimits:
    """Linux containers and cluster jobs, read from fake trees on any OS."""

    def test_a_v2_limit_is_the_limit_less_what_is_used(self, tmp_path):
        root = _cgroup_v2(tmp_path, {"/job": (4 * GIB, 1 * GIB, 0)})
        proc = _proc(tmp_path, "0::/job\n")

        assert memory_guard._cgroup_headroom(proc, root) == 3 * GIB

    def test_reclaimable_page_cache_is_not_counted_as_used(self, tmp_path):
        """A source read once fills the cache; the kernel gives it back."""
        root = _cgroup_v2(tmp_path, {"/job": (4 * GIB, 4 * GIB, 3 * GIB)})
        proc = _proc(tmp_path, "0::/job\n")

        assert memory_guard._cgroup_headroom(proc, root) == 3 * GIB

    def test_no_limit_is_none(self, tmp_path):
        root = _cgroup_v2(tmp_path, {"/job": ("max", 1 * GIB, 0)})
        proc = _proc(tmp_path, "0::/job\n")

        assert memory_guard._cgroup_headroom(proc, root) is None

    def test_a_tighter_parent_wins(self, tmp_path):
        root = _cgroup_v2(
            tmp_path,
            {
                "/": ("max", 9 * GIB, 0),
                "/slurm": (2 * GIB, 1 * GIB, 0),
                "/slurm/job": (8 * GIB, 1 * GIB, 0),
            },
        )
        proc = _proc(tmp_path, "0::/slurm/job\n")

        assert memory_guard._cgroup_headroom(proc, root) == 1 * GIB

    def test_a_container_reports_a_path_its_mount_does_not_have(self, tmp_path):
        """Inside a container the mount root is the container's own cgroup."""
        root = _cgroup_v2(tmp_path, {"/": (2 * GIB, GIB // 2, 0)})
        proc = _proc(tmp_path, "0::/kubepods/pod1/abc\n")

        assert memory_guard._cgroup_headroom(proc, root) == GIB + GIB // 2

    def test_a_v1_memory_controller(self, tmp_path):
        root = tmp_path / "sys_fs_cgroup"
        folder = root / "memory" / "docker" / "abc"
        folder.mkdir(parents=True)
        (folder / "memory.limit_in_bytes").write_text(f"{3 * GIB}\n")
        (folder / "memory.usage_in_bytes").write_text(f"{2 * GIB}\n")
        (folder / "memory.stat").write_text(f"total_inactive_file {GIB}\n")
        proc = _proc(tmp_path, "4:cpu,cpuacct:/docker/abc\n7:memory:/docker/abc\n")

        assert memory_guard._cgroup_headroom(proc, root) == 2 * GIB

    def test_a_v1_unlimited_value_is_no_limit(self, tmp_path):
        root = tmp_path / "sys_fs_cgroup"
        folder = root / "memory"
        folder.mkdir(parents=True)
        (folder / "memory.limit_in_bytes").write_text("9223372036854771712\n")
        (folder / "memory.usage_in_bytes").write_text(f"{GIB}\n")
        proc = _proc(tmp_path, "7:memory:/\n")

        assert memory_guard._cgroup_headroom(proc, root) is None

    def test_no_cgroup_file_is_none(self, tmp_path):
        assert memory_guard._cgroup_headroom(tmp_path / "absent", tmp_path) is None


class TestTheProcessLimitIsWhatIsFree:
    def test_a_tighter_process_limit_wins_over_the_machine(self, monkeypatch):
        monkeypatch.setattr(memory_guard, "process_memory_headroom", lambda: 2 * GIB)

        assert memory_guard.available_memory_gb() == pytest.approx(2.0)

    def test_no_limit_leaves_the_machine_figure(self, monkeypatch):
        import psutil

        monkeypatch.setattr(memory_guard, "process_memory_headroom", lambda: None)

        free = psutil.virtual_memory().available / GIB
        assert memory_guard.available_memory_gb() == pytest.approx(free, rel=0.2)

    def test_a_probe_that_fails_is_no_limit(self, monkeypatch):
        """A probe must never be the thing that stops a conversion."""

        def broken():
            raise OSError("no such file")

        monkeypatch.setattr(memory_guard, "_job_object_headroom", broken)
        monkeypatch.setattr(memory_guard, "_cgroup_headroom", lambda: 5 * GIB)
        monkeypatch.setattr(memory_guard, "_address_space_headroom", lambda: None)

        assert memory_guard.process_memory_headroom() == 5 * GIB

    def test_the_smallest_limit_wins(self, monkeypatch):
        monkeypatch.setattr(memory_guard, "_job_object_headroom", lambda: 7 * GIB)
        monkeypatch.setattr(memory_guard, "_cgroup_headroom", lambda: 3 * GIB)
        monkeypatch.setattr(memory_guard, "_address_space_headroom", lambda: None)

        assert memory_guard.process_memory_headroom() == 3 * GIB


#: Run in a child, because a process cannot leave a job once it is in one.
_JOB_CHILD = textwrap.dedent("""
    import ctypes, ctypes.wintypes as wt, sys
    import psutil
    from thyra.resampling.memory_guard import available_memory_gb, mass_axis_refusal

    k = ctypes.WinDLL("kernel32")
    k.CreateJobObjectW.restype = k.GetCurrentProcess.restype = wt.HANDLE
    k.SetInformationJobObject.argtypes = [wt.HANDLE, ctypes.c_int, ctypes.c_void_p, wt.DWORD]
    k.AssignProcessToJobObject.argtypes = [wt.HANDLE, wt.HANDLE]

    class Ext(ctypes.Structure):
        _fields_ = [
            ("t", ctypes.c_longlong * 2), ("flags", wt.DWORD), ("ws", ctypes.c_size_t * 2),
            ("apl", wt.DWORD), ("aff", ctypes.c_size_t), ("pc", wt.DWORD * 2),
            ("io", ctypes.c_ulonglong * 6), ("proc_lim", ctypes.c_size_t),
            ("job_lim", ctypes.c_size_t), ("pk", ctypes.c_size_t * 2),
        ]

    flag = int(sys.argv[1])
    info = Ext(flags=flag)
    limit = psutil.Process().memory_info().private + 2 * 1024**3
    if flag == 0x200:
        info.job_lim = limit
    else:
        info.proc_lim = limit
    job = k.CreateJobObjectW(None, None)
    assert k.SetInformationJobObject(job, 9, ctypes.byref(info), ctypes.sizeof(info))
    assert k.AssignProcessToJobObject(job, k.GetCurrentProcess())
    print(round(available_memory_gb(), 2))
    print(mass_axis_refusal(18_000_000))
    """)


@pytest.mark.skipif(sys.platform != "win32", reason="Windows job objects")
class TestAWindowsJobLimit:
    """F194's own case: a job that lets the process commit 2 GB more."""

    @pytest.mark.parametrize("flag", [0x100, 0x200], ids=["process", "job"])
    def test_the_guard_reads_the_job_and_refuses(self, flag):
        repo = Path(__file__).resolve().parents[3]
        result = subprocess.run(
            [sys.executable, "-c", _JOB_CHILD, str(flag)],
            capture_output=True,
            text=True,
            cwd=repo,
            timeout=120,
        )

        assert result.returncode == 0, result.stderr
        free, refusal = result.stdout.strip().splitlines()[-2:]
        assert float(free) == pytest.approx(2.0, abs=0.1)
        assert "18,000,000 bins" in refusal and "GB free" in refusal

    def test_outside_a_limiting_job_there_is_no_job_headroom(self):
        """This process may sit in a job (a terminal's), but not one with a cap."""
        headroom = memory_guard._job_object_headroom()
        assert headroom is None or headroom > 0
