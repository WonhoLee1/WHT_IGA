import os, psutil, numpy as np, scipy.sparse as sp
from wht_iga.kratos_backend import load_kratos
mods = load_kratos()
from wht_iga.linalg import pardiso_solve
A = sp.diags([4.0] * 50) + sp.eye(50, k=1) + sp.eye(50, k=-1)
x = pardiso_solve(A.tocsr(), np.ones(50), kind="auto")
print("pardiso ok", float(np.abs(A @ x - 1).max()))
seen = sorted({m.path for m in psutil.Process().memory_maps()
               if any(k in os.path.basename(m.path).lower() for k in ("mkl_rt", "mkl_core", "iomp", "libomp", "vcomp"))})
print("\n".join(seen))
