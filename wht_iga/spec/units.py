"""Consistent unit systems and physical-range plausibility checks for materials.

A study declares one consistent unit system; material numbers are interpreted
in it. The check converts Young's modulus and density to SI and rejects values
outside wide physical bounds. That catches the common factor-1e6/1e12 mix-ups
(SI steel typed into a mm-N-tonne-s study, or the reverse for density) without
rejecting real soft or dense materials.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

UnitSystemName = Literal["mm-N-tonne-s", "m-N-kg-s", "mm-kN-kg-ms"]


@dataclass(frozen=True)
class UnitSystem:
    """Scale factors from one consistent unit system to SI."""

    name: str
    length_to_m: float
    stress_to_pa: float
    density_to_kg_m3: float


UNIT_SYSTEMS: dict[str, UnitSystem] = {
    # mm, N, tonne, s: stress MPa, density tonne/mm^3 (= 1e12 kg/m^3)
    "mm-N-tonne-s": UnitSystem("mm-N-tonne-s", 1e-3, 1e6, 1e12),
    # SI: stress Pa, density kg/m^3
    "m-N-kg-s": UnitSystem("m-N-kg-s", 1.0, 1.0, 1.0),
    # mm, kN, kg, ms: stress kN/mm^2 = GPa, density kg/mm^3 (= 1e9 kg/m^3)
    "mm-kN-kg-ms": UnitSystem("mm-kN-kg-ms", 1e-3, 1e9, 1e9),
}

#: Young's modulus bounds in Pa: soft gels (~1 kPa) up to above diamond (~1.2 TPa).
E_BOUNDS_PA = (1.0e3, 2.0e12)
#: density bounds in kg/m^3 (when > 0): aerogels (~1) up to above osmium (~22.6e3).
RHO_BOUNDS_KG_M3 = (1.0e-1, 3.0e4)


def check_material_units(units: str, E: float, density: float, label: str = "material") -> list[str]:
    """Return human-readable problems (empty list = plausible) for E and density."""
    us = UNIT_SYSTEMS[units]
    problems: list[str] = []
    e_pa = E * us.stress_to_pa
    if not E_BOUNDS_PA[0] <= e_pa <= E_BOUNDS_PA[1]:
        problems.append(
            f"{label}: E={E:g} in {units} is {e_pa:.3e} Pa, outside the physical range "
            f"[{E_BOUNDS_PA[0]:.0e}, {E_BOUNDS_PA[1]:.0e}] Pa -- check the unit system")
    if density > 0.0:
        rho = density * us.density_to_kg_m3
        if not RHO_BOUNDS_KG_M3[0] <= rho <= RHO_BOUNDS_KG_M3[1]:
            problems.append(
                f"{label}: density={density:g} in {units} is {rho:.3e} kg/m^3, outside the "
                f"physical range [{RHO_BOUNDS_KG_M3[0]:g}, {RHO_BOUNDS_KG_M3[1]:g}] kg/m^3 "
                f"-- check the unit system")
    return problems
