import numpy as np
import matplotlib.pyplot as plt

# packing density
p = 0.7
# v : volume of the monomer (not the aggregate)
R = 0.59  # radius in nm (!= R_hydro)
v = 4 / 3 * np.pi * (R * 1E-9) ** 3


gamma_sl = 0.018
# gamma = np.linspace(0.1E-3, 1E-1, 1000)
max_n = 50
n = np.linspace(0.1, max_n, 1000)

T = 293  # K
eta = 1E-3  # visosity Pa.s
R_hydro = 0.63  # nm
k_B = 1.38E-23
D = k_B * T / (6 * np.pi * eta * R_hydro * 1E-9)  # Diffusion

# Stefan-Skapski-Turnbull formula
beta = 0.514
lambda_heat = np.power(v, 2/3) * gamma_sl / beta
Ce = 1/v*np.exp(-lambda_heat/(k_B * T))
print(Ce)

S = 400

W_surface = gamma_sl * np.power(36*np.pi, 1/3)*np.power(v/p, 2/3)*np.power(n, 2/3)
W_surface /= k_B * T
W_volume = n * k_B*T*np.log(S)
W_volume /= k_B * T
W = W_surface - W_volume


theta = gamma_sl * np.power(36*np.pi, 1/3) * np.power(v/p, 2/3) / (k_B * T) # Related to surface energy of aggregate
print("theta : ", theta)

n_etoile = np.power(2*theta/(3 * np.log(S)), 3)
n_etoile_pos = np.searchsorted(n, n_etoile)
print("n_etoile :", n_etoile)
W_etoile = gamma_sl * np.power(36*np.pi, 1/3)*np.power(v/p, 2/3)*np.power(n_etoile, 2/3) - n_etoile * k_B*T*np.log(S)
W_etoile /= k_B * T

plt.plot(n, W_surface, "k--", label="surface")
plt.plot(n, -W_volume, "b-.", label="volume")
plt.plot(n, W, label="total")
plt.vlines(n_etoile, 0, np.max(W), colors="k", linestyles="dotted", alpha=0.7)
plt.hlines(0, 0, max_n, alpha=0.2)
plt.xlabel("n")
plt.ylabel(r"$W(n)/k_B T$")
plt.ylim(-200,250)
plt.xlim(0, max_n)
# plt.fill_between(n, 0, W, where=W>0 , interpolate=True, alpha=0.2)
plt.fill_between(n[:n_etoile_pos], W[:n_etoile_pos], alpha=0.2)
plt.arrow(n_etoile, -150, 0, -30, head_width=0.5, head_length=10, overhang=0)
plt.arrow(3.5, W_etoile, -2.5, 0, length_includes_head=True, head_width=5, head_length=0.8, overhang=0)
plt.plot(n_etoile, W_etoile, "ro", alpha=0.5)
plt.legend()
plt.show()

plt.plot(n, W, label="total")
plt.vlines(n_etoile, 0, np.max(W))
plt.hlines(0, 0, max_n, alpha=0.2)
plt.xlabel("n")
plt.ylabel(r"$W(n)/k_B T$")
plt.show()