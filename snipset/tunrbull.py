import numpy as np
import matplotlib.pyplot as plt


# packing density
p = 0.7
# v : volume of the monomer (not the aggregate)
R_nm = 0.59  # radius in nm (!= R_hydro)
v = 4 / 3 * np.pi * (R_nm * 1E-9) ** 3

T = 293  # K
k_B = 1.38E-23
Na = 6.02E23

gamma_sl = np.linspace(1E-3, 7E-2, 1000)

# Stefan-Skapski-Turnbull formula
beta = 0.514
lambda_heat = np.power(v, 2/3) * gamma_sl / beta
Ce = 1/v*np.exp(-lambda_heat/(k_B * T))
# From particle to mol
Ce /= Na
# From per m^3 to liter
Ce /= 1E3

plt.plot(gamma_sl, lambda_heat/(k_B*T))
plt.xlabel("Surface tension J/m²")
plt.ylabel("lambda_heat (J)")
plt.title("Conductivité thermique")
plt.show()

plt.plot(gamma_sl, Ce)
plt.xlabel("Surface tension J/m²")
plt.ylabel("Equilibrium concentration (solubility) in mol/L")
plt.title("Turnbull equation")
plt.show()