# -*- coding: utf-8 -*-
"""
Created on Thu Oct 29 10:44:01 2020

@author: utilisateur
"""

import threading
import time
import ctypes
import os, sys
import matplotlib.pyplot as plt
import numpy as np
import configparser

class Core:

    def __init__(self, dll_name="CTN.dll"):
        self.view = None
        self.dll_name = dll_name
        self.is_using_C_dll = True
        self.monitoring_delay = 0.5
        self.thread_calculate = threading.Thread(name='Calculate', target=self.calculate)
        self.thread_monitor = threading.Thread(name='Monitor', target=self.monitor)
        self.nb_pt_monitoring = None
        self.is_stop = False
        self.is_calculating = False
        self.num_monitoring = 0
        # initialized with default values
        self.param_dict = {"time_step_ns": 1E-7,
                           "end_time_calculation_ms": 1E-3,
                           "monitoring_time_step_micros": 1E-5,
                           "v_solvent_ml": 1E-2,
                           "p_compacity": 0.7,
                           "r_nm": 0.59,
                           "gamma_sl": 0.006,
                           "temperature": 297,
                           "eta": 0.001,
                           "r_hydro_nm": 0.64,
                           "cini": 1E-3,
                           "ce_m": 1E-6,
                           "is_use_turnbull": 0
                           }
        self.config_parser_ini = None
        self.result_dict = {}
        self.result_t0_dict = {}
        self.set_param(self.param_dict)
        self.load_dll()

    def set_view(self, view):
        self.view = view

    def load_dll(self):
        dll_name = self.dll_name
        cur_dir = os.path.abspath(os.path.dirname(__file__))
        dll_file_path = os.path.normpath(os.path.join(cur_dir, dll_name))
        file_exist = os.path.isfile(dll_file_path)
        # print(dll_file_path)
        self.dll = ctypes.CDLL(dll_file_path)

        # CTN(short *is_stop, unsigned int* num_monitoring, double* S_data, double* J_data, double* aggregate_number, double* mean_size, double* std_size, double* n_star_data, double* concentration_monomer)
        self.ctn_fct = self.dll.CTN
        self.ctn_fct.argtypes = [ctypes.POINTER(ctypes.c_short), ctypes.POINTER(ctypes.c_uint),
                                 ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                 ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                 ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double),
                                 ctypes.POINTER(ctypes.c_double)]

    def set_param(self, param_dict):
        self.param_dict = param_dict
        self.monitoring_time_step_micros = self.param_dict["end_time_calculation_ms"]
        self.nb_pt_monitoring = int(self.param_dict["end_time_calculation_ms"] / self.param_dict["monitoring_time_step_micros"]*1E3)


    def allocate_memory_for_results(self):
        if self.is_using_C_dll:
            array_of_double = ctypes.c_double * self.nb_pt_monitoring
            self.S = array_of_double()
            self.J = array_of_double()
            self.nb_aggregates_array = array_of_double()
            self.mean_size = array_of_double()
            self.std_size = array_of_double()
            self.n_star = array_of_double()
            self.monomer_attachment_ns = array_of_double()

    def create_dict_result(self):
        self.result_dict = {}
        nb_pt = self.num_monitoring.value
        self.result_dict["nb_pt"] = nb_pt
        self.result_dict["t"] = np.arange(nb_pt)*self.monitoring_time_step_micros
        self.result_dict["S"] = np.array(self.S)[:nb_pt]
        self.result_dict["J"] = np.array(self.J)[:nb_pt]
        self.result_dict["nb_aggregates_array"] = np.array(self.nb_aggregates_array)[:nb_pt]
        self.result_dict["mean_size"] = np.array(self.mean_size)[:nb_pt]
        self.result_dict["std_size"] = np.array(self.std_size)[:nb_pt]
        self.result_dict["n_star"] = np.array(self.n_star)[:nb_pt]
        self.result_dict["monomer_attachment_ns"] = np.array(self.monomer_attachment_ns)[:nb_pt]

    def stop_calculation(self):
        if self.thread_calculate.is_alive():
            self.is_stop.value = 1
            self.is_calculating = False
            self.thread_calculate.join(timeout=0.5)
            self.thread_monitor.join(timeout=0.5)

    def launch_calculation(self):
        self.thread_calculate = threading.Thread(name='Calculate', target=self.calculate)
        self.thread_monitor = threading.Thread(name='Monitor', target=self.monitor)

        self.thread_calculate.start()
        self.is_calculating = True
        self.thread_monitor.start()

    def monitor_display(self):
        # print(self.num_monitoring.value)
        pass
        # can't graph something on a matplotlib figure embedded in tk from a different thread.
        # t = np.arange(self.num_monitoring) * self.monitoring_time_step_s
        # plt.plot(t, self.S)
        # plt.show()

    def monitor(self):
        while self.is_calculating:
            self.monitor_display()
            time.sleep(self.monitoring_delay)

    def import_ini_file(self, ini_file_path="CTN.ini"):
        self.param_dict = {}
        config = configparser.ConfigParser()
        config.sections()
        config.read(ini_file_path)
        for key in config["CTN"]:
            self.param_dict[key] = config["CTN"][key]

    def export_ini_file(self, ini_file_path="CTN.ini"):
        self.config_parser_ini = configparser.ConfigParser()
        self.config_parser_ini['CTN'] = self.param_dict
        with open(ini_file_path, 'w') as configfile:
            self.config_parser_ini.write(configfile)

    def plot_data(self):
        if self.view is None:
            pass
        else:
            self.view.plot_data(self.result_dict)

    def get_output_val_at_t0(self):
        S = self.param_dict["cini"] / self.param_dict["ce_m"]
        self.result_t0_dict["S"] = S

        T = self.param_dict["temperature"]
        gamma_sl = self.param_dict["gamma_sl"]
        time_step_ns = self.param_dict["time_step_ns"]
        Ce = self.param_dict["ce_m"]
        eta = self.param_dict["eta"]
        R_hydro_nm = self.param_dict["r_hydro_nm"]

        p = self.param_dict["p_compacity"]
        R_nm = self.param_dict["r_nm"]
        v = 4 / 3 * np.pi * (R_nm * 1E-9) ** 3
        V_L = self.param_dict["v_solvent_ml"] * 1E-3

        k_B = 1.38E-23
        D = k_B * T / (6 * np.pi * eta * R_hydro_nm * 1E-9)  # Diffusion

        theta = gamma_sl * np.power(36 * np.pi, 1 / 3) * np.power(v / p, 2 / 3) / (k_B * T)  # Related to surface energy of aggregate
        self.result_t0_dict["theta"] = theta

        n_star = np.power(2 * theta / (3 * np.log(S)), 3)

        J = (V_L * 1E-3) * time_step_ns * np.power(36*np.pi, 1/6) * D * (Ce * 6.02E23 / 1E3)/(np.power(v**2/p**2, 1/3) * np.sqrt(theta))*S*np.log(S) * np.exp(-4*theta**3/(27*np.log(S)*np.log(S)))

        W_surface = gamma_sl * np.power(36 * np.pi, 1 / 3) * np.power(v / p, 2 / 3) * np.power(n_star, 2 / 3)
        W_surface /= k_B * T
        W_volume = n_star * k_B * T * np.log(S)
        W_volume /= k_B * T
        Wn_star = W_surface - W_volume

        Ce_per_m3 = Ce * 6.02E23 / 1E3
        delta_n = 4 * np.pi * D * Ce_per_m3 * np.power(3 * n_star * v / (4 * np.pi * p), 1.0 / 3) * (S - 1) * 1E-9
        monomer_attachment_ns = 1 / delta_n

        self.result_t0_dict["J"] = J
        self.result_t0_dict["nstar"] = n_star
        self.result_t0_dict["Wnstar"] = Wn_star
        self.result_t0_dict["monomer_attachment_ns"] = monomer_attachment_ns

    def calculate(self):
        # Create ini file
        self.export_ini_file()
        self.set_param(self.param_dict)
        self.allocate_memory_for_results()
        self.is_calculating = True
        # Launch calculation from external dll with array and variable passed by reference.
        self.is_stop = ctypes.c_short(0)
        self.num_monitoring = ctypes.c_uint(0)
        print(self.num_monitoring)
        # CTN(short *is_stop, unsigned int* num_monitoring, double* S_data, double* J_data, double* aggregate_number, double* mean_size, double* std_size, double* n_star_data, double* concentration_monomer)
        # this function is a blocking process for this thread, the variable self.is_stop can stop remotely the main loop of the calculation
        self.ctn_fct(ctypes.byref(self.is_stop), ctypes.byref(self.num_monitoring), self.S, self.J,
                     self.nb_aggregates_array, self.mean_size, self.std_size, self.n_star,
                     self.monomer_attachment_ns)

        # self.num_monitoring = self.num_monitoring.value
        # self.create_dict_result()
        self.is_calculating = False
        # print("End calculation")


if __name__ == "__main__":
    core = Core()
    core.load_dll()
    core.launch_calculation()
