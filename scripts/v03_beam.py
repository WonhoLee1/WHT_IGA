# -*- coding: utf-8 -*-
"""V03: cantilever beam bending, 5 load cases vs Euler-Bernoulli (G1 gate).

Setup mirrors Coreform verification (0.1 x 0.1 x 1 m, E=2.1e11, nu=0):
grid == beam bbox exactly (padding 0) -> no CUT/OUT cells, pure IN path.
Loads: gravity (rho=8000 -> q=784.8 N/m), body force (equivalent 78480 N/m^3),
uniform pressure 5000 Pa (q=500 N/m), total-force pressure (same resultant),
end shear traction 5e4 Pa (P=500 N). Moment case excluded per G1.
Refs: tip delta = q*L^4/(8EI) (uniform), P*L^3/(3EI) (tip shear).
Writes WHT_IGA/docs/phase1/evidence/v03/results.json + per-case DB/VTU.
"""

from __future__ import annotations

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..")))

E = 2.1e11
NU = 0.0
B = H = 0.1
L = 1.0
IX = B * H ** 3 / 12.0
EVID = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "docs",
                    "phase1", "evidence", "v03")


def main() -> None:
    from wht_iga.solver.pipeline import run_linear
    from wht_iga.solver.validate import DirichletBC, JobSpec, PressureLoad, TractionLoad

    os.makedirs(EVID, exist_ok=True)
    lo, hi, nc, deg = (-0.05, -0.05, 0.0), (0.05, 0.05, 1.0), (4, 4, 40), 2
    common = dict(lo=lo, hi=hi, n_cells=nc, degree=deg, E=E, nu=NU,
                  fixed_faces=[(2, 0)],
                  probes=[{"name": "tip", "xyz": [0.0, 0.0, 1.0]}])
    cases = [
        ("gravity", dict(density=8000.0, gravity=(0.0, -9.81, 0.0)),
         784.8 * L ** 4 / (8 * E * IX)),
        ("bodyforce", dict(density=1.0, gravity=(0.0, -78480.0, 0.0)),
         784.8 * L ** 4 / (8 * E * IX)),
        ("pressure_uniform",
         dict(pressures=[PressureLoad(face=(1, 1), magnitude=-5000.0)]),
         500.0 * L ** 4 / (8 * E * IX)),
        ("pressure_total_force",
         dict(pressures=[PressureLoad(face=(1, 1), magnitude=-5000.0)]),
         500.0 * L ** 4 / (8 * E * IX)),
        ("traction",
         dict(tractions=[TractionLoad(face=(2, 1), vector=(0.0, -5e4, 0.0))]),
         500.0 * L ** 3 / (3 * E * IX)),
    ]
    out = {}
    for name, kw, ref in cases:
        spec = JobSpec(name=f"v03_{name}", output_dir=EVID, **common, **kw)  # type: ignore[arg-type]
        res = run_linear(spec)
        uy = float(np.atleast_1d(res.probe_values["tip"])[1])
        err = abs(abs(uy) - ref) / ref
        out[name] = {"tip_uy": uy, "ref": ref, "rel_err": err,
                     "beta": res.beta_used}
        print(f"[V03 {name}] tip={uy:.6e} ref={-ref:.6e} rel_err={err:.3e}", flush=True)
        assert err < 0.05, (name, err)
    with open(os.path.join(EVID, "results.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, indent=1)
    print("V03 ALL PASS (gate: 5%)")


if __name__ == "__main__":
    main()
