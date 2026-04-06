# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
"""NUMA auto-bind helper.

Detects multi-node NUMA topology via ``/sys/devices/system/node/`` and binds
the calling process (and its forked children) to the **least-loaded** NUMA
node.  On single-node or non-NUMA systems this is a no-op that produces no
output.

Binding mechanism tried in priority order
-----------------------------------------
1. ``ctypes`` → ``libnuma.so.1`` — calls ``numa_run_on_node()`` and
   ``numa_set_preferred()`` directly.  No compilation required; works on
   Python 3.12 + where the ``numa`` PyPI package no longer builds.  Requires
   the system package ``libnuma1`` (Debian/Ubuntu) or ``numactl-libs`` (RHEL).
2. ``numa`` Python library (``pip install numa``) — older Python (≤ 3.11) only.
3. ``numactl --cpunodebind=N --localalloc`` — re-exec via ``os.execvpe``
   (replaces the process image; covers child subprocesses too; sets memory
   policy correctly).  Only attempted when ``numactl`` is found in PATH.
4. ``os.sched_setaffinity()`` — Python stdlib, zero external deps.  Sets CPU
   affinity only; memory policy is left to the OS first-touch heuristic.
   CPU affinity *is* inherited by ``fork()``-ed child processes.

HPC / job-scheduler integration
---------------------------------
Autobind is **silently skipped** when any of the following environment
variables are set (they indicate the scheduler has already configured
placement):

* ``GOMP_CPU_AFFINITY``  (GCC OpenMP)
* ``KMP_AFFINITY``        (Intel MKL / OpenMP)
* ``SLURM_CPU_BIND``      (Slurm)
* ``SGE_BINDING``         (SGE / UGE)
* ``_NUMA_AUTOBIND_DONE`` (our own re-exec loop-guard)

``OMP_PROC_BIND`` is honoured: ``close``/``master``/``true`` → *local*;
``spread`` → *spread*.
"""
# This work © 2025-2026 by Jiří Zahradník and Martin Mokrejš
# (First Medical Faculty - Charles University in Prague) is licensed under
# Creative Commons Attribution 4.0 International.

from __future__ import annotations

import ctypes
import ctypes.util
import glob
import os
import shutil
import sys
from typing import NamedTuple

# ---------------------------------------------------------------------------
# Internal data type
# ---------------------------------------------------------------------------


class _NUMANode(NamedTuple):
    """Per-node NUMA information read from sysfs."""

    index:      int
    cpuset:     frozenset[int]
    free_bytes: int


# ---------------------------------------------------------------------------
# sysfs helpers
# ---------------------------------------------------------------------------

def _parse_cpulist(cpulist_str: str) -> frozenset[int]:
    """Parse a kernel cpulist string (e.g. ``'0-3,8-11'``) into a frozenset."""
    cpus: set[int] = set()
    for part in cpulist_str.strip().split(','):
        if '-' in part:
            lo, hi = part.split('-', 1)
            cpus.update(range(int(lo), int(hi) + 1))
        elif part.strip():
            cpus.add(int(part))
    return frozenset(cpus)


def _read_node_free_bytes(node_idx: int) -> int:
    """Return free memory in bytes for *node_idx* from sysfs ``meminfo``."""
    path = f'/sys/devices/system/node/node{node_idx}/meminfo'
    try:
        with open(path, encoding='ascii') as fh:
            for line in fh:
                # "Node 0 MemFree:       7938564 kB"
                if 'MemFree' in line:
                    return int(line.split(':')[1].split()[0]) * 1024
    except (OSError, ValueError):
        pass
    return 0


def detect_nodes() -> list[_NUMANode]:
    """Return a list of :class:`_NUMANode` for each online NUMA node.

    Returns an **empty list** when there is only one node (or the sysfs path
    does not exist) — meaning NUMA optimisation is not applicable.
    """
    node_dirs = sorted(glob.glob('/sys/devices/system/node/node[0-9]*'))
    if len(node_dirs) <= 1:
        return []

    nodes: list[_NUMANode] = []
    for d in node_dirs:
        idx = int(os.path.basename(d)[4:])          # "node0" → 0
        cpulist_path = os.path.join(d, 'cpulist')
        try:
            with open(cpulist_path, encoding='ascii') as fh:
                cpulist_str = fh.read()
        except OSError:
            cpulist_str = ''
        cpuset = _parse_cpulist(cpulist_str)
        if not cpuset:
            continue                                  # memory-only node: skip
        free = _read_node_free_bytes(idx)
        nodes.append(_NUMANode(index=idx, cpuset=cpuset, free_bytes=free))
    return nodes


# ---------------------------------------------------------------------------
# Environment / scheduler checks
# ---------------------------------------------------------------------------

#: Variables set by HPC schedulers or the user to configure CPU placement.
#: If any of these are non-empty we must not override the layout.
_SCHEDULER_BIND_VARS = (
    'GOMP_CPU_AFFINITY',   # GCC OpenMP
    'KMP_AFFINITY',         # Intel MKL / OpenMP
    'SLURM_CPU_BIND',       # Slurm
    'SGE_BINDING',          # SGE / UGE
    '_NUMA_AUTOBIND_DONE',  # our own re-exec loop-guard
)


def _scheduler_already_bound() -> bool:
    """Return ``True`` if a job scheduler or user already set CPU affinity."""
    return any(os.environ.get(v) for v in _SCHEDULER_BIND_VARS)


def _omp_proc_bind_override() -> str | None:
    """Return ``'local'`` or ``'spread'`` from ``OMP_PROC_BIND``, or ``None``."""
    val = os.environ.get('OMP_PROC_BIND', '').lower()
    if val in ('close', 'master', 'true'):
        return 'local'
    if val == 'spread':
        return 'spread'
    return None


# ---------------------------------------------------------------------------
# Binding back-ends
# ---------------------------------------------------------------------------

def _bind_via_libnuma_ctypes(node: _NUMANode) -> bool:
    """In-process bind via ``ctypes`` → ``libnuma.so.1`` (no compilation needed).

    Works on Python 3.12 + where the ``numa`` PyPI extension module no longer
    builds.  Requires the system package ``libnuma1`` (Debian/Ubuntu) or
    ``numactl-libs`` (RHEL/Fedora).
    """
    libname = ctypes.util.find_library('numa') or 'libnuma.so.1'
    try:
        lib = ctypes.CDLL(libname, use_errno=True)
    except OSError:
        return False
    try:
        if lib.numa_available() < 0:
            return False          # kernel has no NUMA support (CONFIG_NUMA not set)
        lib.numa_run_on_node(ctypes.c_int(node.index))
        lib.numa_set_preferred(ctypes.c_int(node.index))
        return True
    except (AttributeError, OSError):
        return False


def _bind_via_numa_pylib(node: _NUMANode) -> bool:
    """Try in-process bind via the ``numa`` PyPI package (Python ≤ 3.11 only).

    Sets both CPU affinity and memory allocation policy.
    """
    try:
        import numa  # type: ignore[import-untyped]  # optional dependency
    except ImportError:
        return False
    try:
        if not numa.available():
            return False
        numa.run_on_node(node.index)
        numa.set_preferred(node.index)
        return True
    except Exception:  # pylint: disable=broad-except
        return False


def _bind_via_numactl_reexec(node: _NUMANode, verbose: bool) -> None:
    """Re-exec the current process under ``numactl`` (sets memory policy too).

    Calls ``os.execvpe`` which **replaces** the process image — this function
    never returns on success.  The guard variable ``_NUMA_AUTOBIND_DONE``
    prevents an infinite re-exec loop.
    """
    numactl = shutil.which('numactl')
    if numactl is None:
        return                                        # numactl not in PATH; caller falls back

    # Reconstruct the full command:
    # • installed entry point (no .py suffix): sys.argv[0] is the wrapper
    #   executable, which already has the correct shebang.
    # • standalone .py script: we must prepend sys.executable.
    if sys.argv[0].endswith('.py'):
        cmd = [sys.executable] + sys.argv
    else:
        cmd = list(sys.argv)

    env = {**os.environ, '_NUMA_AUTOBIND_DONE': '1'}
    exec_argv = [numactl, f'--cpunodebind={node.index}', '--localalloc'] + cmd

    if verbose:
        free_mib = node.free_bytes // (1024 * 1024)
        print(
            f"[numa_bind] re-exec under numactl --cpunodebind={node.index} "
            f"--localalloc  (node free: {free_mib:,} MiB)",
            file=sys.stderr,
        )
    os.execvpe(numactl, exec_argv, env)
    # os.execvpe replaces the process image; the lines below are unreachable.


def _bind_via_sched_affinity(node: _NUMANode) -> bool:
    """Set CPU affinity via ``os.sched_setaffinity()`` (Python stdlib).

    Sets CPU affinity only (no memory policy).  Affinity is inherited by
    child processes created with ``fork()``.
    """
    try:
        os.sched_setaffinity(0, node.cpuset)
        return True
    except (AttributeError, OSError):
        return False


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def autobind(mode: str = 'local', verbose: bool = True) -> int | None:
    """Detect NUMA topology and bind the process to the best node.

    Parameters
    ----------
    mode:
        ``'local'``  — compact: all threads/processes on the single NUMA node
        with the most free RAM (default).

        ``'spread'`` — leave OS to use all available nodes (no binding applied;
        useful for aggregate-bandwidth experiments).

        ``'none'``   — disable autobind entirely.

    verbose:
        When ``True`` (default) print a one-line summary to *stderr* when
        binding is actually applied.  Nothing is printed on single-node systems
        or when skipped due to scheduler variables.

    Returns
    -------
    int | None
        The bound node index, or ``None`` when no binding was applied.
    """
    if mode == 'none':
        return None

    # Honour OMP_PROC_BIND if explicitly set by the scheduler/user.
    omp_mode = _omp_proc_bind_override()
    if omp_mode is not None:
        mode = omp_mode

    if _scheduler_already_bound():
        return None

    nodes = detect_nodes()
    if not nodes:
        return None                   # single-node / non-NUMA: silent no-op

    if mode == 'spread':
        if verbose:
            total_cpus = sum(len(nd.cpuset) for nd in nodes)
            print(
                f"[numa_bind] {len(nodes)} NUMA nodes, {total_cpus} CPUs — "
                f"spread policy: no binding applied.",
                file=sys.stderr,
            )
        return None

    # ── local / compact ── pick node with most free memory ──────────────
    best = max(nodes, key=lambda nd: nd.free_bytes)

    # Back-end 1: ctypes → libnuma.so.1 (in-process, CPU + memory policy).
    bound_method: str | None = None
    if _bind_via_libnuma_ctypes(best):
        bound_method = 'ctypes/libnuma'
    elif _bind_via_numa_pylib(best):
        bound_method = 'numa-pylib'
    else:
        # Back-end 2: numactl re-exec (CPU + memory policy, replaces image).
        # Does NOT return on success (os.execvpe).
        _bind_via_numactl_reexec(best, verbose=verbose)
        # Back-end 3: os.sched_setaffinity (CPU affinity only, inherited by fork).
        if _bind_via_sched_affinity(best):
            bound_method = 'sched_setaffinity (CPU only)'

    if bound_method is not None:
        if verbose:
            _print_bound_msg(best, nodes, method=bound_method)
        return best.index
    # Nothing worked — only print advisory on genuine multi-NUMA systems.
    if verbose:
        print(
            f"[numa_bind] {len(nodes)} NUMA nodes detected but no binding method "
            f"is available (install `libnuma1` + `numactl`). "
            f"Using OS default placement.",
            file=sys.stderr,
        )
    return None


def _print_bound_msg(best: _NUMANode, nodes: list[_NUMANode], method: str) -> None:
    free_mib = best.free_bytes // (1024 * 1024)
    total_mib = sum(nd.free_bytes for nd in nodes) // (1024 * 1024)
    cpus_str = ','.join(str(c) for c in sorted(best.cpuset))
    print(
        f"[numa_bind] Bound to NUMA node {best.index} via {method}  "
        f"(free: {free_mib:,}/{total_mib:,} MiB, CPUs: {cpus_str})",
        file=sys.stderr,
    )
