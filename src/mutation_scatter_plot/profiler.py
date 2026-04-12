# -*- coding: utf-8 -*-
# vim: set fileencoding=utf-8 ts=4 sw=4 expandtab :
import datetime
import os
import subprocess
import threading
import time


class ResourceProfiler(threading.Thread):
    """Background thread to sample CPU and RAM (RSS) usage of this script and its children."""

    def __init__(self, interval: float = 5.0):
        super().__init__(daemon=True)
        self.interval = interval
        self.pid = os.getpid()
        self.lock = threading.Lock()
        self._has_started = False

        # State
        self.active_phase = None
        self.phase_start_time = None

        # Metrics for current phase
        self.phase_cpu_samples = []
        self.phase_ram_samples = []

        # Current instant load
        self.current_cpu = 0.0
        self.current_ram_gb = 0.0

        self.psutil = None
        try:
            import psutil
            self.psutil = psutil
        except ImportError:
            pass

    def _get_metrics_psutil(self) -> tuple[float, float]:
        try:
            main_proc = self.psutil.Process(self.pid)
            procs = [main_proc]
            procs.extend(main_proc.children(recursive=True))

            # Baseline trigger
            for p in procs:
                try:
                    p.cpu_percent(interval=None)
                except (self.psutil.NoSuchProcess, self.psutil.AccessDenied):
                    pass

            time.sleep(0.5)  # half-second measurement window

            total_cpu = 0.0
            total_rss = 0.0
            for p in procs:
                try:
                    total_cpu += p.cpu_percent(interval=None)
                    total_rss += p.memory_info().rss
                except (self.psutil.NoSuchProcess, self.psutil.AccessDenied):
                    pass
            return total_cpu, total_rss / (1024**3)
        except Exception:
            # Fallback defensively if psutil tree topology evaporates mid-sweep
            return 0.0, 0.0

    def _get_metrics_ps(self) -> tuple[float, float]:
        try:
            # -p PID (parent) and --ppid PID (children)
            cmd = ["ps", "--no-headers", "-o", "pcpu,rss", "-p", str(self.pid), "--ppid", str(self.pid)]
            out = subprocess.check_output(cmd, encoding='ascii', stderr=subprocess.DEVNULL)

            total_cpu = 0.0
            total_rss_kb = 0.0
            for line in out.strip().split('\n'):
                parts = line.strip().split()
                if len(parts) == 2:
                    total_cpu += float(parts[0])
                    total_rss_kb += float(parts[1])
            return total_cpu, total_rss_kb / (1024**2)  # ps RSS is in KB
        except (subprocess.CalledProcessError, FileNotFoundError, ValueError):
            return 0.0, 0.0

    def start(self):
        """Start the profiler thread (idempotent; safe to call multiple times)."""
        with self.lock:
            if self._has_started:
                return
            self._has_started = True
            super().start()

    def run(self):
        """Background loop: sample CPU and RAM at ``self.interval`` seconds."""
        while True:
            t0 = time.time()
            if self.psutil:
                cpu, ram_gb = self._get_metrics_psutil()
            else:
                cpu, ram_gb = self._get_metrics_ps()

            with self.lock:
                # Discard 0.0 readings if the OS jittered
                if cpu > 0 or ram_gb > 0:
                    self.current_cpu = cpu
                    self.current_ram_gb = ram_gb
                if self.active_phase is not None:
                    self.phase_cpu_samples.append(self.current_cpu)
                    self.phase_ram_samples.append(self.current_ram_gb)

            # Sleep remaining interval (account for the 0.5s psutil measurement window)
            elapsed = time.time() - t0
            if elapsed < self.interval:
                time.sleep(self.interval - elapsed)

    def mark_phase_start(self, phase_name: str):
        """Begin a new named phase, resetting sample accumulators."""
        with self.lock:
            self.active_phase = phase_name
            self.phase_start_time = datetime.datetime.now()
            self.phase_cpu_samples = []
            self.phase_ram_samples = []

    def get_current_load_str(self) -> str:
        """Return a human-readable string of the latest CPU and RAM readings."""
        with self.lock:
            return f"  [Load: {self.current_ram_gb:.1f} GB RAM | {self.current_cpu:.0f}% CPU]"

    def pop_phase_summary(self) -> str | None:
        """End the current phase and return a formatted summary string.

        Returns ``None`` if no phase was active.
        """
        with self.lock:
            if self.active_phase is None or not self.phase_start_time:
                return None

            name = self.active_phase
            duration = datetime.datetime.now() - self.phase_start_time

            secs = int(duration.total_seconds())
            hours, remainder = divmod(secs, 3600)
            minutes, seconds = divmod(remainder, 60)
            if hours > 0:
                dur_str = f"{hours}h {minutes:02d}m"
            elif minutes > 0:
                dur_str = f"{minutes}m {seconds:02d}s"
            else:
                dur_str = f"{seconds}s"

            if self.phase_cpu_samples:
                peak_cpu = max(self.phase_cpu_samples)
                avg_cpu = sum(self.phase_cpu_samples) / len(self.phase_cpu_samples)
                peak_ram = max(self.phase_ram_samples)
            else:
                peak_cpu = avg_cpu = self.current_cpu
                peak_ram = self.current_ram_gb

            summary = (f"         ↳ Profiler Summary ({name}): Peak RAM = {peak_ram:.1f} GB | "
                       f"Peak CPU = {peak_cpu:.0f}% | Avg CPU = {avg_cpu:.0f}% | Duration = {dur_str}")

            self.active_phase = None
            return summary


# Global instance imported by downstream tools
PROFILER = ResourceProfiler(interval=5.0)
