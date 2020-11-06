import numpy as np
import matplotlib.pyplot  as plt



# 4/3 np.pi R_aggregate**3 = V -> volume aggregate
# v_molecule_inside_aggregate_nm3 : volume of one molecule inside the aggregate is v/p where is the v is the volume of one molecule and p the packing density
# Finally, n * v_molecule_inside_aggregate_nm3 = V where n is the number of molecule inside the aggregate


R_molecule_nm = 0.59
p = 0.7
v_molecule_nm3 = 4.0/3 * np.pi * R_molecule_nm**3
v_molecule_inside_aggregate_nm3 = v_molecule_nm3/ p

R_aggregate_nm = np.linspace(1,50, 1000)
V_aggregate_nm3 = 4.0/3 * np.pi * R_aggregate_nm**3
n = V_aggregate_nm3 / v_molecule_inside_aggregate_nm3





idx_1nm = np.searchsorted(R_aggregate_nm, 0.5)
idx_2nm = np.searchsorted(R_aggregate_nm, 1)
idx_5nm = np.searchsorted(R_aggregate_nm, 2.5)
idx_10nm = np.searchsorted(R_aggregate_nm, 5)
idx_20nm = np.searchsorted(R_aggregate_nm, 10)
idx_50nm = np.searchsorted(R_aggregate_nm, 25)


plt.plot(2*R_aggregate_nm, n)
plt.xlabel("Diameter aggregate in nm")
plt.ylabel("Nb of molecule inside the aggregate")

n_values = [1e3, 2e3, 5e3, 1e4, 2e4, 5e4, 1e5, 2e5]
idx_n_values = []
for value in n_values:
    idx_n_values.append(np.searchsorted(n, value))

for idx in idx_n_values:
    plt.vlines(2 * R_aggregate_nm[idx], 0, n[idx], alpha=0.5)
    # plt.hlines(n[idx], 0, 2*R_aggregate_nm[idx], alpha=0.5)
    plt.text(2 * R_aggregate_nm[idx]- 1, n[idx] + 1E4, str(int(2*R_aggregate_nm[idx])))

plt.show()

plt.plot(n, 2*R_aggregate_nm)
plt.ylabel("Diameter aggregate in nm")
plt.xlabel("Nb of molecule inside the aggregate")
plt.show()

