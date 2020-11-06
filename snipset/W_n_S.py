import numpy as np
import matplotlib.pyplot as plt

# packing density
p = 0.7
# v : volume of the monomer (not the aggregate)
R = 0.59  # radius in nm (!= R_hydro)
v = 4 / 3 * np.pi * (R * 1E-9) ** 3


gamma_sl = 0.03
max_n = 100
n = np.linspace(0.1, max_n, 1000)
S = np.linspace(2, 1E6, 1000)


S_mesh, n_mesh = np.meshgrid(S, n)

T = 293  # K
eta = 1E-3  # visosity Pa.s
R_hydro = 0.63  # nm
k_B = 1.38E-23
D = k_B * T / (6 * np.pi * eta * R_hydro * 1E-9)  # Diffusion

"""
W_surface = gamma_sl * np.power(36*np.pi, 1/3)*np.power(v/p, 2/3)*np.power(n, 2/3)
W_volume = n * k_B*T*np.log(S)
W = W_surface - W_volume
"""
W_surface = gamma_sl * np.power(36*np.pi, 1/3)*np.power(v/p, 2/3)*np.power(n_mesh, 2/3)
W_volume = n_mesh * k_B*T*np.log(S_mesh)
W = W_surface - W_volume

theta = gamma_sl * np.power(36*np.pi, 1/3) * np.power(v/p, 2/3) / (k_B * T) # Related to surface energy of aggregate
# print("theta : ", theta)

n_etoile = np.power(2*theta/(3 * np.log(S)), 3)
# print("n_etoile :", n_etoile)

plt.contourf(n, S, W, 30)
plt.plot(n_etoile, S, "r")
plt.xlim(0, max_n)
plt.colorbar()
plt.xlabel("n")
plt.ylabel("S")
#plt.xscale('log')
plt.yscale('log')
plt.show()

