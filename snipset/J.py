import numpy as np
import matplotlib.pyplot  as plt

time_step_s = 1E-7 # 1ns
Na = 6.02E23

C_ini = 1E-3  # mol/L
Ce = 1E-6  # Mol/L
Ce_m3 = Ce * Na / 1E3
C = C_ini  # C(t)

Na = 6.022E23
V_solvent_L = 1 * 1E-5
n_monomer_ini = C_ini * V_solvent_L * Na
n_monomer = n_monomer_ini  # n_mono(t)

# packing density
p = 0.7
# v : volume of the monomer (not the aggregate)
R = 0.59  # radius in nm (!= R_hydro)
v_monomer = 4 / 3 * np.pi * (R * 1E-9) ** 3
v_monomer_nm3 = 4 / 3 * np.pi * (R) ** 3
print("v :", v_monomer)
print("v_nm3 :", v_monomer_nm3)

"""
17 cal / mol / A²
17/4.18 = 4.06 J /mol / A²
4.06E20 J /mol / m²
4/6.02 E23 = 0.66E-3 J / m²
soit 0.66E-3 J / m²
"""
gamma_eau_air = 0.07
gamma_sl = 0.66E-2  # N/m (eau/air) -> (monomer/solvant) 19.5 (±0.6), 18.5 (±0.6) and 17.0 (±0.6) cal mol-1  Å-2 in the solvent condition with a DMSO mole fraction of 0.21, 0.26 and 0.32
gamma_sl = 0.006


T = 293  # K
eta = 1E-3  # visosity Pa.s

R_hydro = 0.63  # nm

k_B = 1.38E-23
D = k_B * T / (6 * np.pi * eta * R_hydro * 1E-9)  # Diffusion

S = np.linspace(5, 1E4, 1000)

theta = gamma_sl * np.power(36*np.pi, 1/3) * np.power(v_monomer / p, 2 / 3) / (k_B * T) # Related to surface energy of aggregate
print("theta : ", theta)

n_etoile = np.power(2*theta/(3 * np.log(S)), 3)
# print("n_etoile : ", n_etoile)
# print(-4*theta**3/(27*np.log(S)*np.log(S)))
# J is in m-3 s-1 (per cubic meter per second)

pow_1_6_36_pi = np.power(36*np.pi, 1/6)
v_over_p_two_third = np.power(v_monomer/p, 2.0/3)
print("v_over_p_two_third  : ", v_over_p_two_third)
exponential_prefactor_constant_part = (V_solvent_L*1E-3) * time_step_s * pow_1_6_36_pi * D * (Ce_m3)/(pow(v_monomer*v_monomer/(p*p), 1/3) * np.sqrt(theta))
print("exponential_prefactor_constant_part : ", exponential_prefactor_constant_part)
J = (V_solvent_L * 1E-3) * time_step_s * np.power(36 * np.pi, 1 / 6) * D * (Ce * 6.02E23 / 1E3) / (np.power(v_monomer ** 2 / p ** 2, 1 / 3) * np.sqrt(theta)) * S * np.log(S) * np.exp(-4 * theta ** 3 / (27 * np.log(S) * np.log(S)))

plt.plot(S, J)
plt.xlabel("S")
plt.ylabel("J")
plt.show()