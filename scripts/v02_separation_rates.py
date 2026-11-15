"""V02 분리 수렴률(e_u, e_s, peak_A) 사후분석 스크립트 (#8, 솔버 재실행 없음).

기존 ``v02_h*_U.npz`` 와 ``results.json`` 만을 읽어 크라운 단면 선분 ``(0, 10, 0) -> (0, 20, 0)``
위 동일 복합 가우스 점에서의 변위 L2 상대오차(``e_u``)와 최대주응력 L2 상대오차(``e_s``)를
분리 산출하고, CSV / JSON / log-log 수렴 플롯을 생성한다.

주의: ``qoi.py`` (사다리꼴) vs ``probes.py`` (복합 가우스) 적분 정의가 별존하므로,
게이트 정의(``probes.py`` 복합 가우스)에 맞춰 일관 평가한다.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

from wht_iga.reference.box import BoxSpline
from wht_iga.solver import probes as PR

EVID_DEFAULT = os.path.join(
    os.path.dirname(os.path.abspath(__file__)),
    "..", "docs", "phase2", "evidence", "v02_meshfile_drop_ext",
)
PROFILE_P0 = (0.0, 10.0, 0.0)
PROFILE_P1 = (0.0, 20.0, 0.0)
REF_ROARK = 4.3175e-6
REF_ABAQUS = 4.34796947956784e-6
FIGSIZE = (7.2, 4.4)
CAVEAT_TEXT = (
    "주의: rate>2는 관측된 자기수렴 기울기이며 초수렴·이론차수 보장을 주장하지 않음. "
    "기준해는 최세 IGA 해(h=0.375), Roark/Abaqus는 peak 참조용."
)


@dataclass
class LevelField:
    """단일 격자 레벨의 저장 해와 요약 지표."""

    h: float
    box: BoxSpline
    U: np.ndarray
    peak_A: float
    err_A_Roark: float
    err_A_Abaqus: float
    n_cut: int
    n_fallback: int
    crown_fallback: int
    u_line: np.ndarray
    s_line: np.ndarray


def load_v02_levels(evid_dir: str, p0: tuple[float, float, float],
                    p1: tuple[float, float, float]) -> tuple[list[LevelField], np.ndarray]:
    """``results.json`` 및 ``v02_h*_U.npz`` 를 읽어 동일 선분 가우스점에서 변위·응력을 평가한다."""
    res_path = os.path.join(evid_dir, "results.json")
    with open(res_path, "r", encoding="utf-8") as fh:
        meta = json.load(fh)
    Pp, wp = PR.line_gauss_points(p0, p1)
    levels: list[LevelField] = []
    for row in meta["rows"]:
        h_val = float(row["h"])
        npz_path = os.path.join(evid_dir, f"v02_h{h_val:g}_U.npz")
        data = np.load(npz_path)
        box = BoxSpline(
            lo=tuple(float(v) for v in data["lo"]),
            hi=tuple(float(v) for v in data["hi"]),
            n_cells=tuple(int(v) for v in data["n_cells"]),
            degree=int(data["degree"]),
        )
        U = np.asarray(data["U"], dtype=float)
        u_line = box.evaluate(U.reshape(-1, 3), Pp)
        s_line = PR.max_principal_profile(box, U, 1.0, 0.0, Pp)
        if "profile" in row:
            np.testing.assert_allclose(s_line, np.asarray(row["profile"], dtype=float), rtol=1e-11, atol=1e-18)
        peak_A = float(row["peak_A"])
        levels.append(
            LevelField(
                h=h_val,
                box=box,
                U=U,
                peak_A=peak_A,
                err_A_Roark=abs(peak_A - REF_ROARK) / REF_ROARK,
                err_A_Abaqus=abs(peak_A - REF_ABAQUS) / REF_ABAQUS,
                n_cut=int(row["n_cut"]),
                n_fallback=int(row["n_fallback"]),
                crown_fallback=int(row["crown_fallback"]),
                u_line=u_line,
                s_line=s_line,
            )
        )
    levels.sort(key=lambda lv: lv.h, reverse=True)
    return levels, wp


def _vector_rel_l2(a: np.ndarray, b: np.ndarray, w: np.ndarray) -> float:
    """벡터장 ``(N, 3)`` 의 가중 L2 상대오차 ``||a - b||_L2 / ||b||_L2``."""
    num = float(np.sum(w * np.sum((a - b) ** 2, axis=1)))
    den = float(np.sum(w * np.sum(b**2, axis=1)))
    return float(np.sqrt(num / den))


def compute_separation_table(levels: list[LevelField], wp: np.ndarray) -> tuple[list[dict], dict]:
    """변위 분리율 ``e_u`` 와 응력 분리율 ``e_s`` 및 log-log 기울기(all / last3)를 계산한다."""
    ref = levels[-1]
    hs = np.array([lv.h for lv in levels[:-1]], dtype=float)
    e_u = [_vector_rel_l2(lv.u_line, ref.u_line, wp) for lv in levels[:-1]] + [0.0]
    e_s = [PR.rel_l2(lv.s_line, ref.s_line, wp) for lv in levels[:-1]] + [0.0]

    eu_arr = np.asarray(e_u[:-1], dtype=float)
    es_arr = np.asarray(e_s[:-1], dtype=float)
    rows: list[dict] = []
    for lv, u_err, s_err in zip(levels, e_u, e_s, strict=True):
        rows.append(
            {
                "h": lv.h,
                "e_u": u_err,
                "e_s": s_err,
                "peak_A": lv.peak_A,
                "err_A_Roark": lv.err_A_Roark,
                "err_A_Abaqus": lv.err_A_Abaqus,
                "n_cut": lv.n_cut,
                "n_fallback": lv.n_fallback,
                "crown_fallback": lv.crown_fallback,
            }
        )
    summary = {
        "n_levels": len(levels),
        "finest_reference_h": ref.h,
        "gate_criteria": {"peak_A_err_max": 0.01, "e_s_rate_min": 1.8},
        "e_u": {
            "rate_all": float(np.polyfit(np.log(hs), np.log(eu_arr), 1)[0]),
            "rate_last3": float(np.polyfit(np.log(hs[-3:]), np.log(eu_arr[-3:]), 1)[0]),
        },
        "e_s": {
            "rate_all": float(np.polyfit(np.log(hs), np.log(es_arr), 1)[0]),
            "rate_last3": float(np.polyfit(np.log(hs[-3:]), np.log(es_arr[-3:]), 1)[0]),
        },
        "finest_peak_A_err_Roark": ref.err_A_Roark,
        "finest_peak_A_err_Abaqus": ref.err_A_Abaqus,
        "integrator_note": "probes.py composite Gauss (qoi.py trapezoidal left untouched per #8 scope)",
        "caveat": CAVEAT_TEXT,
    }
    return rows, summary


def write_csv(rows: list[dict], out_csv: str) -> None:
    """표준 스키마 CSV를 UTF-8(BOM 없음)로 저장한다."""
    cols = ["h", "e_u", "e_s", "peak_A", "err_A_Roark", "err_A_Abaqus",
            "n_cut", "n_fallback", "crown_fallback"]
    with open(out_csv, "w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=cols)
        writer.writeheader()
        writer.writerows(rows)


def plot_separation_rates(rows: list[dict], summary: dict, out_png: str) -> None:
    """``e_u`` 와 ``e_s`` 의 log-log 수렴 플롯을 생성한다."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    try:
        import koreanize_matplotlib  # noqa: F401
    except ImportError:
        plt.rcParams["font.family"] = "Malgun Gothic"
        plt.rcParams["axes.unicode_minus"] = False

    plt.rcParams["font.size"] = 9
    hs = np.array([r["h"] for r in rows[:-1]], dtype=float)
    eu = np.array([r["e_u"] for r in rows[:-1]], dtype=float)
    es = np.array([r["e_s"] for r in rows[:-1]], dtype=float)

    fig, ax = plt.subplots(figsize=FIGSIZE)
    ax.loglog(
        hs, eu, "o-",
        label=f"변위 L2 (e_u): rate_all={summary['e_u']['rate_all']:.2f}, last3={summary['e_u']['rate_last3']:.2f}",
    )
    ax.loglog(
        hs, es, "s-",
        label=f"최대주응력 L2 (e_s): rate_all={summary['e_s']['rate_all']:.2f}, last3={summary['e_s']['rate_last3']:.2f}",
    )
    ax.invert_xaxis()
    ax.grid(True, which="both", linestyle="--", alpha=0.4)
    ax.set_xlabel("셀 크기 h (mm)")
    ax.set_ylabel("최세 IGA 해 (h=0.375) 대비 상대 L2 오차")
    ax.set_title("V02 Roark 원공 평판 분리 수렴률 (솔버 재실행 없음)")
    ax.legend(loc="upper left", bbox_to_anchor=(1.02, 1.0), ncol=1)
    fig.text(0.01, 0.01, CAVEAT_TEXT, fontsize=8)
    fig.tight_layout(rect=[0.0, 0.04, 1.0, 1.0])
    fig.savefig(out_png, dpi=180, bbox_inches="tight")
    plt.close(fig)


def main(evid_dir: str = EVID_DEFAULT) -> dict:
    """V02 분리 수렴률 CSV, JSON, PNG를 생성한다."""
    levels, wp = load_v02_levels(evid_dir, PROFILE_P0, PROFILE_P1)
    rows, summary = compute_separation_table(levels, wp)
    out_csv = os.path.join(evid_dir, "v02_separation_rates.csv")
    out_json = os.path.join(evid_dir, "v02_separation_summary.json")
    out_png = os.path.join(evid_dir, "v02_separation_rates.png")
    write_csv(rows, out_csv)
    with open(out_json, "w", encoding="utf-8") as fh:
        json.dump({"rows": rows, "summary": summary}, fh, indent=1, ensure_ascii=False)
    plot_separation_rates(rows, summary, out_png)
    print(f"[V02 #8] e_u rate_all={summary['e_u']['rate_all']:.3f} rate_last3={summary['e_u']['rate_last3']:.3f}")
    print(f"[V02 #8] e_s rate_all={summary['e_s']['rate_all']:.3f} rate_last3={summary['e_s']['rate_last3']:.3f}")
    print(f"[V02 #8] peak_A finest Roark err={summary['finest_peak_A_err_Roark']:.3%}, Abaqus err={summary['finest_peak_A_err_Abaqus']:.3%}")
    return summary


if __name__ == "__main__":
    main()
