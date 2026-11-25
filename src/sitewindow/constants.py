"""Physical constants (CODATA 2018) and unit factors used throughout ``sitewindow``.

All energies internal to the package are in electronvolts per adsorbed entity unless a
function states otherwise; temperatures in kelvin; pressures in bar (absolute).
"""

from __future__ import annotations

# CODATA 2018 exact or recommended values
K_B_J: float = 1.380649e-23  # J K^-1
K_B_EV: float = 8.617333262e-5  # eV K^-1
R_GAS: float = 8.314462618  # J mol^-1 K^-1
N_A: float = 6.02214076e23  # mol^-1
H_PLANCK: float = 6.62607015e-34  # J s
HBAR: float = 1.054571817e-34  # J s
E_CHARGE: float = 1.602176634e-19  # C
C_LIGHT: float = 299792458.0  # m s^-1
M_E_KG: float = 9.1093837015e-31  # kg
M_E_C2_EV: float = 510998.95  # eV
R_E_M: float = 2.8179403262e-15  # classical electron radius, m
ALPHA_FS: float = 7.2973525693e-3  # fine-structure constant
AMU_KG: float = 1.66053906660e-27  # kg

# Masses (amu)
M_H_AMU: float = 1.00794
M_D_AMU: float = 2.014102
M_H2_AMU: float = 2.01588
M_C_AMU: float = 12.0107
M_O_AMU: float = 15.9994
M_C12_AMU: float = 12.0
M_C13_AMU: float = 13.003355

# Unit factors
EV_TO_KJ_PER_MOL: float = 96.48533212  # kJ mol^-1 per eV
KJ_PER_MOL_TO_EV: float = 1.0 / EV_TO_KJ_PER_MOL
EV_TO_J: float = E_CHARGE
BAR_TO_PA: float = 1.0e5
ATM_TO_BAR: float = 1.01325
MPA_TO_BAR: float = 10.0
KPA_TO_BAR: float = 0.01
TORR_TO_BAR: float = 1.01325 / 760.0
PSI_TO_BAR: float = 0.0689475729

# Graphene structural constants
GRAPHENE_CC_BOND_A: float = 1.42  # angstrom
GRAPHENE_LATTICE_A: float = 2.46  # angstrom
GRAPHENE_AREAL_DENSITY_NM2: float = 38.2  # carbon atoms per nm^2 (one layer)
GRAPHENE_SSA_M2_G: float = 2630.0  # both faces of a single layer

# Molar volume at STP (0 degC, 1 atm) for cm3(STP) conversions
V_MOLAR_STP_CM3: float = 22413.97
