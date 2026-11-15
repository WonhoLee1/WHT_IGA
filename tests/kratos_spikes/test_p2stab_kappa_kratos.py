"""P2-STAB Kratos spike K-A2: kappa tool anchored to Kratos GetSystemMatrix.

Builds a small unconstrained V01-reduced Kratos model, extracts the system
matrix (GetSystemMatrix + scipy_conversion_tools, v2 §5.6), drops the
strong-Dirichlet rows/cols, and checks (a) agreement with the reference
assembler (≤ 1e-10, S04 pattern) and (b) cond_eigsh against the dense
reference (≤ 5%). This anchors wht_iga.stabilize.cond to real Kratos output.
"""

from __future__ import annotations

import numpy as np
import pytest

from wht_iga.kratos_backend import load_kratos
from wht_iga.kratos_backend.elements import (
    attach_cell_elements,
    equation_permutation,
    solve_linear,
    system_matrix_csr,
)
from wht_iga.kratos_backend.geometry import make_nurbs_volume
from wht_iga.reference.assemble_elastic import assemble_stiffness
from wht_iga.reference.box import BoxSpline
from wht_iga.stabilize import cond as Q

pytestmark = pytest.mark.kratos


def test_kratos_kappa_matches_reference():
    mods = load_kratos()
    KM = mods["KratosMultiphysics"]
    SM = mods["StructuralMechanicsApplication"]
    Emod, nu = 2.5, 0.3
    box = BoxSpline(lo=(0.0, 0.0, 0.0), hi=(1.0, 1.0, 1.0),
                    n_cells=(2, 2, 2), degree=1)
    model = KM.Model()
    mp = model.CreateModelPart("P2STAB_K")
    mp.AddNodalSolutionStepVariable(KM.DISPLACEMENT)
    mp.AddNodalSolutionStepVariable(KM.REACTION)
    mp.ProcessInfo.SetValue(KM.DOMAIN_SIZE, 3)
    vol, _ = make_nurbs_volume(KM, mp, box)
    cl = SM.LinearElastic3DLaw()
    props = mp.GetProperties()[1]
    props.SetValue(KM.YOUNG_MODULUS, Emod)
    props.SetValue(KM.POISSON_RATIO, nu)
    props.SetValue(KM.CONSTITUTIVE_LAW, cl)
    attach_cell_elements(KM, mp, vol, box, props)
    # Strong Dirichlet on x=0 face (dropped after extraction).
    ncp = box.n_cp[0]
    for node in mp.Nodes:
        I = node.Id - 1
        a = I % ncp
        if a == 0:
            for var in (KM.DISPLACEMENT_X, KM.DISPLACEMENT_Y,
                        KM.DISPLACEMENT_Z):
                node.Fix(var)
    solve_linear(KM, mp)
    Kk = system_matrix_csr(solve_linear(KM, mp)).tocsr()
    perm = equation_permutation(KM, mp, box)
    inv = np.empty_like(perm)
    inv[perm] = np.arange(len(perm))
    Kk_ref = Kk[inv][:, inv]
    Kref = assemble_stiffness(box, Emod, nu)
    free_cps = [I for I in range(box.n_ctrl) if I % ncp != 0]
    free = np.sort(np.concatenate([np.array(free_cps) * 3 + d
                                   for d in range(3)]))
    np.testing.assert_allclose(Kk_ref[free][:, free].toarray(),
                               Kref[free][:, free].toarray(),
                               rtol=0, atol=1e-10)
    kappa_k = Q.cond_eigsh(Kk_ref[free][:, free])[0]
    kappa_r = Q.cond_dense(Kref[free][:, free].toarray())
    print(f"K-A2: kratos kappa={kappa_k:.4e} reference kappa={kappa_r:.4e}")
    assert abs(kappa_k - kappa_r) / kappa_r < 0.05
