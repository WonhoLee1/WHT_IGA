"""V06 runner for G2 evidence: out_mode drop/ghost, per-run RSS sampling.

Usage: python run_v06_drop.py <level> <p> [out_mode] [n_geo]
level 0/1/2 -> n_cells 10/20/40; n_geo default 24/96/192 (as the ghost-mode evidence: v06_l0_momfix, v06_l1_r2, v06_l2).
Writes <evid>/L<level>_p<p>_<mode>.json with the error row and RSS peak.
"""
import json, os, sys, threading, time

import psutil

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, "..", "..", "..", "..")))  # WHT_IGA

from tests.kratos_spikes.test_p2mesh_v06 import _csv_tols, run_v06  # noqa: E402

def main() -> None:
    level, p = int(sys.argv[1]), int(sys.argv[2])
    mode = sys.argv[3] if len(sys.argv) > 3 else "drop"
    n = {0: 10, 1: 20, 2: 40}[level]
    n_geo = int(sys.argv[4]) if len(sys.argv) > 4 else {0: 24, 1: 96, 2: 192}[level]

    proc = psutil.Process()
    peak = {"rss": 0}
    stop = threading.Event()


    def sample():
        while not stop.is_set():
            rss = proc.memory_info().rss
            peak["rss"] = max(peak["rss"], rss)
            print(f"[mem] rss={rss / 2**30:.2f}GB", flush=True)
            stop.wait(15)


    threading.Thread(target=sample, daemon=True).start()
    row = run_v06((n, n, n), p, n_geo, ext_theta=0.2, out_mode=mode)
    stop.set()
    mi = proc.memory_info()
    row["rss_peak_gb_sampled"] = peak["rss"] / 2**30
    row["peak_wset_gb"] = getattr(mi, "peak_wset", 0) / 2**30
    tols = _csv_tols()
    row["tol"] = {q: tols[(p, level, q)] for q in ("u", "ur", "srr", "stt")}
    row["pass"] = {q: bool(row[q] <= row["tol"][q]) for q in ("u", "ur", "srr", "stt")}
    out = os.path.join(HERE, f"L{level}_p{p}_{mode}.json")
    with open(out, "w", encoding="utf-8") as f:
        json.dump(row, f, indent=1)
    print(f"RES L{level} p={p} {mode}: " + " ".join(f"{q}={row[q]:.3e}({'ok' if row['pass'][q] else 'FAIL'})"
          for q in ("u", "ur", "srr", "stt")) + f" peak_wset={row['peak_wset_gb']:.2f}GB t={row['runtime_s']:.0f}s",
          flush=True)


if __name__ == "__main__":
    main()
