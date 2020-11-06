# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 21:09:32 2020

@author: utilisateur
"""


import numpy as np
import numba
import matplotlib.pyplot as plt

# @numba.jit(nopython=True, nogil=True)
def CTN():
    # Delta t
    pas_tps = 1E-7
    
    # 
    monitoring_time_step = 1E-5
    #FIXME
    ratio_time_step_monitoring = int(monitoring_time_step / pas_tps)


    tps_fin = 1E-2
    nb_pas_tps = int(tps_fin / pas_tps)
    tps_total = 0
    nb_monitoring_point = int(tps_fin/monitoring_time_step)
    
    C_ini = 1E-3 #mol/L
    Ce = 1E-6   # Mol/L
    C = C_ini #C(t)
    
    Na = 6.022E23
    V_L = 1 * 1E-3 # 1 mL
    n_monomer_ini = C_ini * V_L * Na
    n_monomer = n_monomer_ini # n_mono(t)
    
    # packing density
    p = 0.7
    # v : volume of the monomer (not the aggregate)
    R = 0.59  # radius in nm (!= R_hydro)
    v = 4/3 * np.pi * (R*1E-9)**3
    
    gamma_sl = 0.07 #N/m (eau/air) -> (monomer/solvant) 19.5 (±0.6), 18.5 (±0.6) and 17.0 (±0.6) cal mol-1  Å-2 in the solvent condition with a DMSO mole fraction of 0.21, 0.26 and 0.32
    gamma_sl = 5E-3

    T = 293 #K
    eta = 1E-3  # visosity Pa.s

    R_hydro = 0.63 # nm
    
    k_B = 1.38E-23
    D = k_B * T /  (6*np.pi*eta*R_hydro*1E-9) # Diffusion
    
    # Les agregats ne sont caractérisés que par une seule valeur : leur nombre d'atome, on cree donc une simple liste.
    max_nb_agregat = int(1E7)
    agregats = np.zeros(max_nb_agregat)
    nb_aggregates = np.zeros(nb_pas_tps, np.int)
    size_aggregates = np.zeros(nb_pas_tps)
    nb_agregate = 0

    #FIXME le +1
    S_data = np.zeros(nb_monitoring_point)
    J_data = np.zeros(nb_monitoring_point)
    aggregate_number = np.zeros(nb_monitoring_point)
    aggregate_mean_size = np.zeros(nb_monitoring_point)
    aggregate_std_size = np.zeros(nb_monitoring_point)
    n_etoile_data = np.zeros(nb_monitoring_point)
    concentration_monomer = np.zeros(nb_monitoring_point)
    
    num_iteration = 0
    num_monitoring = 0
    # When J is not an integer
    residual_J = 0
    J_cumul = 0

    theta = gamma_sl * np.power(36 * np.pi, 1 / 3) * np.power(v / p, 2 / 3) / (k_B * T)
    prefactor_J = (V_L * 1E-3) * pas_tps * np.power(36 * np.pi, 1 / 6) * D * (Ce * 6.02E23 / 1E3) / (
                np.power(v ** 2 / p ** 2, 1 / 3) * np.sqrt(theta))
    for num_iteration in range(nb_pas_tps):
        tps_total += pas_tps
      
        # recalculer la supersaturation
        S = C / Ce
        log_S = np.log(S)
      
        # W = gamma_sl * np.power(36*np.pi,1/3) * np.power(v/p * n,2/3) - n * k_b * T * np.log(S)
         # Related to surface energy of aggregate
        #print("theta", theta)
        n_etoile = np.power(2*theta/(3 * log_S), 3)
        #print("n_etoile", n_etoile)
        # nucleation
        # J in m-3 s-1 (per cubic meter per second
        J = prefactor_J*S*log_S * np.exp(-4*theta**3/(27*log_S**2))

        J_cumul += J
        nb_new_aggregate = int(J_cumul)
        J_cumul -= nb_new_aggregate

        nb_aggregates[num_iteration] = nb_new_aggregate
        size_aggregates = n_etoile + 1

        for n in range(nb_new_aggregate):
            # Il y a un nouveau noyau de nucleation
            # On ajoute un nouvel agregat de taille n_etoile
            agregats[nb_agregate] = n_etoile + 1
            n_monomer -= n_etoile
            nb_agregate += 1
      
        # growth
        for agregat_n in range(nb_agregate):
            # Concentration in m-3
            delta_n = 4 * np.pi * D * (Ce * 6.02E23 / 1E3) * np.power(3*agregats[agregat_n]*v/(4*np.pi*p), 1/3) * (S-1) * pas_tps
            if delta_n < 0 :
                dummy = 1
            agregats[agregat_n] += delta_n
            n_monomer -= delta_n
      
        # Actualiser la concentration (en mol/L)
        C = (n_monomer / 6.02E23) / V_L
      
        if num_iteration%ratio_time_step_monitoring == 0:
            #print("num_iteration = ", num_iteration)
            #print("num_monitoring = ", num_monitoring)
            # monitor the relevant data
            print(nb_agregate)
            S_data[num_monitoring] = S
            J_data[num_monitoring] = J
            aggregate_number[num_monitoring] = nb_agregate
            if nb_agregate > 0:
                aggregate_mean_size[num_monitoring] = np.mean(agregats[0:nb_agregate])
                aggregate_std_size[num_monitoring] = np.std(agregats[0:nb_agregate])
            else:
                aggregate_mean_size[num_monitoring] = 0
                aggregate_std_size[num_monitoring] = 0

            n_etoile_data[num_monitoring] = n_etoile
            concentration_monomer[num_monitoring] = C
            # print("Calculation %d pourcents" % int(num_monitoring*100.0 / nb_monitoring_point))
            num_monitoring += 1
      
    return agregats, S_data, J_data, aggregate_number, aggregate_mean_size, aggregate_std_size, n_etoile_data , concentration_monomer

#TODO CTN parameters
# packing density
p = 0.7
# v : volume of the monomer (not the aggregate)
R = 0.59  # radius in nm (!= R_hydro)
v = 4/3 * np.pi * (R*1E-9)**3

pas_tps = 1E-6

# 
monitoring_time_step = 1E-5
ratio_time_step_monitoring = int(monitoring_time_step / pas_tps)

tps_fin = 1E-2
nb_pas_tps = int(tps_fin / pas_tps)
tps_total = 0
nb_monitoring_point = int(tps_fin/monitoring_time_step)

# aggregate_radius 
agregats, S_data, J_data, aggregate_number, aggregate_mean_size, aggregate_std_size, n_etoile_data , concentration_monomer = CTN()    

V_aggregate = np.array(agregats) * v / p
R_aggregate_nm = np.power(3/4 * V_aggregate / np.pi, 1/3) * 1E9
R_aggregate_mean_nm = np.mean(R_aggregate_nm)
R_aggregate_std_nm = np.std(R_aggregate_nm)

time_micros = np.arange(nb_monitoring_point) * monitoring_time_step * 1E6 #in µs

plt.plot(time_micros, S_data, "ro")
plt.xlabel("time in µs")
plt.ylabel("Supersaturation")
plt.title("Evolution Supersaturation")
plt.show()

plt.plot(time_micros, J_data, "ro")
plt.xlabel("time in µs")
plt.ylabel("J Nucleation Rate")
plt.title("Evolution Nucleation rate")
plt.show()

plt.plot(time_micros, aggregate_number, "ro")
plt.xlabel("time in µs")
plt.ylabel("Number of aggregates ")
plt.title("Evolution du nombre d'aggregat")
plt.show()


# plt.errorbar(time_micros, aggregate_mean_size, yerr=aggregate_std_size)
plt.plot(time_micros, aggregate_mean_size)
plt.xlabel("time in µs")
plt.ylabel("Taille moyenne (en nombre d'atome n) des agrégats")
plt.title("Evolution taille moyenne (en nombre d'atome n) des agrégats")
plt.show()

# plt.errorbar(time_micros, R_aggregate_mean_nm, yerr=R_aggregate_std_nm)
plt.plot(time_micros, R_aggregate_mean_nm)
plt.xlabel("time in µs")
plt.ylabel("Mean Radius aggregate (nm)")
plt.show()

plt.plot(time_micros, aggregate_number, "ro")
plt.xlabel("time in µs")
plt.ylabel("Number of aggregates ")
plt.title("Evolution du nombre d'aggregat")
plt.show()

plt.plot(time_micros, n_etoile_data, "ro")
plt.xlabel("time in µs")
plt.ylabel("Number of aggregates ")
plt.title("Evolution de la taille du nucleus (n*)")
plt.show()

plt.plot(time_micros, concentration_monomer, "ro")
plt.xlabel("time in µs")
plt.ylabel("Concentration monomer")
plt.title("Evolution de la concentration en monomere")
plt.show()

  