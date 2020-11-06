import numpy as np
import matplotlib.pyplot  as plt

time_step = 1E-9 # 1ns

C_ini = 1E-3  # mol/L
Ce = 1E-6  # Mol/L
C = C_ini  # C(t)

Na = 6.022E23
V_mL = 1 * 1E-3  # 1 mL
n_monomer_ini = C_ini * V_mL * Na
n_monomer = n_monomer_ini  # n_mono(t)

# packing density
p = 0.7
# v : volume of the monomer (not the aggregate)
R = 0.59  # radius in nm (!= R_hydro)
v = 4 / 3 * np.pi * (R * 1E-9) ** 3
print("v :", v)

"""
17 cal / mol / A²
17/4.18 = 4.06 J /mol / A²
4.06E20 J /mol / m²
4/6.02 E23 = 0.66E-3 J / m²
soit 0.66E-3 J / m²
"""
gamma_eau_air = 0.07
gamma_sl = 0.66E-2  # N/m (eau/air) -> (monomer/solvant) 19.5 (±0.6), 18.5 (±0.6) and 17.0 (±0.6) cal mol-1  Å-2 in the solvent condition with a DMSO mole fraction of 0.21, 0.26 and 0.32



T = 293  # K
eta = 1E-3  # visosity Pa.s

R_hydro = 0.63  # nm

k_B = 1.38E-23
D = k_B * T / (6 * np.pi * eta * R_hydro * 1E-9)  # Diffusion
print("D : ", D)
n = np.linspace(0.1, 1E6, 1000)

S = 1000
pas_tps = 1E-6

Ce_per_m3 = Ce * 6.02E23 / 1E3
print("Ce_per_m3  : ", Ce_per_m3)
delta_n = 4 * np.pi * D * Ce_per_m3 * np.power(3*n*v/(4*np.pi*p), 1.0/3) * (S-1) * pas_tps




plt.plot(n, delta_n)
plt.xlabel("n")
plt.ylabel("delta_n")
plt.title("Evolution delta n with Ce=" + str(Ce) + " and S=" + str(S) + " in " + str(pas_tps*1E6) + "µs")
plt.show()
