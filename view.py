import tkinter as tk
from tkinter import ttk
import matplotlib.pyplot as plt
import threading, time
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg, NavigationToolbar2Tk)
from tkinter import filedialog, messagebox, simpledialog
import os

class View():
    def __init__(self, core):
        self.core = core
        self.core.set_view(self)
        self.root = tk.Tk()
        self.fig_size_x, self.fig_size_y = 3.5,3.5
        self.fig_dpi = 100
        self.sv_input_dict, self.iv_input_dict = {}, {}
        # self.param_dict = {}
        self.sv_output_t0_dict, self.sv_output_t_dict = {}, {}
        self.create_widget()
        self.create_menu()
        self.create_string_value_dict()
        self.update_t0_value()
        self.thread_core_calculation = threading.Thread(name='Calculate Core', target=self.core_calculation)
        self.calculation_polling_delay = 500

        self.root.protocol("WM_DELETE_WINDOW",
                           self.on_quit)  # Exit when x pressed, notice that its the name of the function 'self.handler' and not a method call self.handler()

    def run(self):
        self.root.title("pyCNT")
        self.root.deiconify()
        self.root.deiconify()
        self.root.mainloop()

    def on_quit(self):
        # paramFile = open('param.ini', 'w')
        # paramFile.write(self.saveDir)
        self.root.destroy()
        self.root.quit()

    def create_widget(self):
        # Simulation
        self.input_output_frame = tk.Frame()
        self.input_output_frame.pack(side=tk.TOP, fill="both", expand=True)
        self.input_param = tk.LabelFrame(self.input_output_frame, text="Input parameters",
                                               borderwidth=4, padx=5, pady=5, font='helvetica15')
        self.input_param.pack(side=tk.LEFT, fill="both", expand=True)

        # Simulation
        self.frame_simul_param = tk.LabelFrame(self.input_param, text="Simulation",
                                               borderwidth=1, padx=5, pady=5,)
        self.frame_simul_param.pack(side=tk.LEFT, fill="both", expand=True)

        ttk.Label(self.frame_simul_param, text='Time step (ns)').grid(row=0, column=0)
        self.time_step_ns_sv = tk.StringVar(value='100')
        e = ttk.Entry(self.frame_simul_param, width=12, justify=tk.CENTER, textvariable=self.time_step_ns_sv)
        e.grid(row=0, column=1)
        e.bind('<Return>', lambda e: self.update_t0_value())

        ttk.Label(self.frame_simul_param, text='Monitor time (µs)').grid(row=1, column=0)
        self.monitor_time_micros_sv = tk.StringVar(value='1')
        e = ttk.Entry(self.frame_simul_param, width=12, justify=tk.CENTER, textvariable=self.monitor_time_micros_sv)
        e.grid(row=1, column=1)
        e.bind('<Return>', lambda e: self.update_t0_value())

        ttk.Label(self.frame_simul_param, text='End time (ms)').grid(row=2, column=0)
        self.end_time_calculation_ms = tk.StringVar(value='1')
        e = ttk.Entry(self.frame_simul_param, width=12, justify=tk.CENTER, textvariable=self.end_time_calculation_ms)
        e.grid(row=2, column=1)
        e.bind('<Return>', lambda e: self.update_t0_value())


        # Solvent
        self.frame_solvent_param = tk.LabelFrame(self.input_param, text="Solvent",
                                               borderwidth=1)
        self.frame_solvent_param.pack(side=tk.LEFT, fill="both", expand=True)

        ttk.Label(self.frame_solvent_param, text='Temperature (K)').grid(row=0, column=0)
        self.temperature_K_sv = tk.StringVar(value="293")
        e = ttk.Entry(self.frame_solvent_param, width=12, justify=tk.CENTER, textvariable=self.temperature_K_sv)
        e.grid(row=0, column=1)
        e.bind('<Return>', lambda e: self.update_t0_value())

        # TODO preset
        ttk.Label(self.frame_solvent_param, text='Viscosity (Pa.s)').grid(row=1, column=0)
        self.viscosity_sv = tk.StringVar(value="1E-3")
        e = ttk.Entry(self.frame_solvent_param, width=12, justify=tk.CENTER, textvariable=self.viscosity_sv)
        e.grid(row=1, column=1)
        e.bind('<Return>', lambda e: self.update_t0_value())

        # Molecule
        self.frame_molecule_param = tk.LabelFrame(self.input_param, text="Molecule",
                                               borderwidth=1)
        self.frame_molecule_param.pack(side=tk.LEFT, fill="both", expand=True)

        ttk.Label(self.frame_molecule_param, text='Packing density').grid(row=0, column=0)
        self.packing_density_p_sv = tk.StringVar(value="0.7")
        e = ttk.Entry(self.frame_molecule_param, width=12, justify=tk.CENTER, textvariable=self.packing_density_p_sv)
        e.grid(row=0, column=1)
        e.bind('<Return>', lambda e: self.update_t0_value())

        ttk.Label(self.frame_molecule_param, text='Radius (nm)').grid(row=1, column=0)
        self.molecule_radius_nm_sv = tk.StringVar(value="0.59")
        e = ttk.Entry(self.frame_molecule_param, width=12, justify=tk.CENTER, textvariable=self.molecule_radius_nm_sv)
        e.grid(row=1, column=1)
        e.bind('<Return>', lambda e: self.update_t0_value())

        ttk.Label(self.frame_molecule_param, text='Hydrodynamic Radius (nm)').grid(row=2, column=0)
        self.molecule_hydro_radius_nm_sv = tk.StringVar(value="0.63")
        e = ttk.Entry(self.frame_molecule_param, width=12, justify=tk.CENTER, textvariable=self.molecule_hydro_radius_nm_sv)
        e.grid(row=2, column=1)
        e.bind('<Return>', lambda e: self.update_t0_value())

        ttk.Label(self.frame_molecule_param, text='Surface tension GammaSl (N/m)').grid(row=3, column=0)
        self.molecule_tension_surface_sv = tk.StringVar(value="0.0066")
        e = ttk.Entry(self.frame_molecule_param, width=12, justify=tk.CENTER, textvariable=self.molecule_tension_surface_sv)
        e.grid(row=3, column=1)
        e.bind('<Return>', lambda e: self.update_t0_value())

        ttk.Label(self.frame_molecule_param, text='Molecule solubility Ce (mlo/L)').grid(row=4, column=0)
        self.molecule_solubility_Ce_sv = tk.StringVar(value="1E-6")
        e = ttk.Entry(self.frame_molecule_param, width=12, justify=tk.CENTER, textvariable=self.molecule_solubility_Ce_sv)
        e.grid(row=4, column=1)
        e.bind('<Return>', lambda e: self.update_t0_value())

        self.is_use_Turnbull_iv = tk.IntVar(value=0)
        tk.Checkbutton(self.frame_molecule_param, text="Use turnbull equation for Ce ?", variable=self.is_use_Turnbull_iv , command=self.update_t0_value).grid(row=5, column=0)

        # Experiment
        self.frame_exp_param = tk.LabelFrame(self.input_param, text="Experiment",
                                               borderwidth=1)
        self.frame_exp_param.pack(side=tk.LEFT, fill="both", expand=True)

        ttk.Label(self.frame_exp_param, text='Monomer concentration (mol/L)').grid(row=0, column=0)
        self.monomer_concentration_M_sv = tk.StringVar(value="1E-3")
        e = ttk.Entry(self.frame_exp_param, width=12, justify=tk.CENTER, textvariable=self.monomer_concentration_M_sv)
        e.grid(row=0, column=1)
        e.bind('<Return>', lambda e: self.update_t0_value())

        ttk.Label(self.frame_exp_param, text='Volume solution (ml)').grid(row=1, column=0)
        self.volume_solution_mL_sv = tk.StringVar(value="1")
        e = ttk.Entry(self.frame_exp_param, width=12, justify=tk.CENTER, textvariable=self.volume_solution_mL_sv)
        e.grid(row=1, column=1)
        e.bind('<Return>', lambda e: self.update_t0_value())

        # Action
        self.frame_actions = tk.LabelFrame(self.root, text="Action",
                                               borderwidth=1)
        self.frame_actions.pack(side=tk.TOP, fill="both", expand=True)

        ttk.Button(self.frame_actions, text="Recalculate", width=15, command=self.recalculate).grid(row=0,column=0)
        ttk.Button(self.frame_actions, text="STOP !", width=15, command=self.stop).grid(row=0,column=1)

        # t = 0 calculation
        self.output_val = tk.LabelFrame(self.input_output_frame, text="Outputs",
                                               borderwidth=2, padx=5, pady=5, font='helvetica15')
        self.output_val.pack(side=tk.LEFT, fill="both", expand=True)

        ttk.Label(self.output_val, text='S').grid(row=1, column=0)
        ttk.Label(self.output_val, text='C_mono').grid(row=2, column=0)
        ttk.Label(self.output_val, text='J').grid(row=3, column=0)
        ttk.Label(self.output_val, text='W(n*)').grid(row=4, column=0)
        ttk.Label(self.output_val, text='theta').grid(row=5, column=0)
        ttk.Label(self.output_val, text='n*').grid(row=6, column=0)
        ttk.Label(self.output_val, text='Attachment time (ns)').grid(row=7, column=0)

        ttk.Label(self.output_val, text='t=0').grid(row=0, column=1)
        ttk.Label(self.output_val, text='t').grid(row=0, column=2)

        self.S_t0_sv = tk.StringVar()
        self.S_t_sv = tk.StringVar()
        e = ttk.Entry(self.output_val, width=10, justify=tk.CENTER, textvariable=self.S_t0_sv)
        e.grid(row=1, column=1)
        e.configure(state='disabled')
        e = ttk.Entry(self.output_val, width=10, justify=tk.CENTER, textvariable=self.S_t_sv)
        e.grid(row=1, column=2)
        e.configure(state='disabled')
        self.J_t0_sv = tk.StringVar()
        self.J_t_sv = tk.StringVar()
        e = ttk.Entry(self.output_val, width=10, justify=tk.CENTER, textvariable=self.J_t0_sv)
        e.grid(row=3, column=1)
        e.configure(state='disabled')
        e = ttk.Entry(self.output_val, width=10, justify=tk.CENTER, textvariable=self.J_t_sv)
        e.grid(row=3, column=2)
        e.configure(state='disabled')
        self.C_t0_sv = tk.StringVar()
        self.C_t_sv = tk.StringVar()
        e = ttk.Entry(self.output_val, width=10, justify=tk.CENTER, textvariable=self.C_t0_sv)
        e.grid(row=2, column=1)
        e.configure(state='disabled')
        e = ttk.Entry(self.output_val, width=10, justify=tk.CENTER, textvariable=self.C_t_sv)
        e.grid(row=2, column=2)
        e.configure(state='disabled')
        self.Wnstar_t0_sv = tk.StringVar()
        self.Wnstar_t_sv = tk.StringVar()
        e = ttk.Entry(self.output_val, width=10, justify=tk.CENTER, textvariable=self.Wnstar_t0_sv)
        e.grid(row=4, column=1)
        e.configure(state='disabled')
        e = ttk.Entry(self.output_val, width=10, justify=tk.CENTER, textvariable=self.Wnstar_t_sv)
        e.grid(row=4, column=2)
        e.configure(state='disabled')
        self.theta_t0_sv = tk.StringVar()
        self.theta_t_sv = tk.StringVar()
        e = ttk.Entry(self.output_val, width=10, justify=tk.CENTER, textvariable=self.theta_t0_sv)
        e.grid(row=5, column=1)
        e.configure(state='disabled')
        e = ttk.Entry(self.output_val, width=10, justify=tk.CENTER, textvariable=self.theta_t_sv)
        e.grid(row=5, column=2)
        e.configure(state='disabled')
        self.nstar_t0_sv = tk.StringVar()
        self.nstar_t_sv = tk.StringVar()
        e = ttk.Entry(self.output_val, width=10, justify=tk.CENTER, textvariable=self.nstar_t0_sv)
        e.grid(row=6, column=1)
        e.configure(state='disabled')
        e = ttk.Entry(self.output_val, width=10, justify=tk.CENTER, textvariable=self.nstar_t_sv)
        e.grid(row=6, column=2)
        e.configure(state='disabled')
        self.attach_time_t0_sv = tk.StringVar()
        self.attach_time_t_sv = tk.StringVar()
        e = ttk.Entry(self.output_val, width=10, justify=tk.CENTER, textvariable=self.attach_time_t0_sv)
        e.grid(row=7, column=1)
        e.configure(state='disabled')
        e = ttk.Entry(self.output_val, width=10, justify=tk.CENTER, textvariable=self.attach_time_t_sv)
        e.grid(row=7, column=2)
        e.configure(state='disabled')







        #Graphes

        self.frame_graphs = tk.LabelFrame(self.root, text="Results graph",
                                               borderwidth=1)
        self.frame_graphs.pack(side="top", fill="both", expand=True)

        # S(t)
        self.frame_graph_S = tk.Frame(self.frame_graphs)
        self.frame_graph_S.grid(row=0, column=0)

        self.figure_S = plt.Figure(figsize=(self.fig_size_x, self.fig_size_y), dpi=self.fig_dpi)
        self.ax_S = self.figure_S.add_subplot(111)

        plt.subplots_adjust(hspace=0)
        self.figure_S.set_tight_layout(True)
        self.canvas_S = FigureCanvasTkAgg(self.figure_S, master=self.frame_graph_S)
        # self.canvas.get_tk_widget().pack(side='top', fill='both', expand=1)

        self.toolbar_S = NavigationToolbar2Tk(self.canvas_S, self.frame_graph_S)
        self.canvas_S._tkcanvas.pack(side='top', fill='both', expand=1)

        # J(T)
        self.frame_graph_J = tk.Frame(self.frame_graphs)
        self.frame_graph_J.grid(row=0, column=1)

        self.figure_J= plt.Figure(figsize=(self.fig_size_x, self.fig_size_y), dpi=self.fig_dpi)
        self.ax_J = self.figure_J.add_subplot(111)

        plt.subplots_adjust(hspace=0)
        self.figure_J.set_tight_layout(True)
        self.canvas_J = FigureCanvasTkAgg(self.figure_J, master=self.frame_graph_J)
        # self.canvas.get_tk_widget().pack(side='top', fill='both', expand=1)

        self.toolbar_J = NavigationToolbar2Tk(self.canvas_J, self.frame_graph_J)
        self.canvas_J._tkcanvas.pack(side='top', fill='both', expand=1)

        # n_star
        self.frame_graph_n_star = tk.Frame(self.frame_graphs)
        self.frame_graph_n_star.grid(row=0, column=2)

        self.figure_n_star= plt.Figure(figsize=(self.fig_size_x, self.fig_size_y), dpi=self.fig_dpi)
        self.ax_n_star = self.figure_n_star.add_subplot(111)

        plt.subplots_adjust(hspace=0)
        self.figure_n_star.set_tight_layout(True)
        self.canvas_n_star = FigureCanvasTkAgg(self.figure_n_star, master=self.frame_graph_n_star)
        # self.canvas.get_tk_widget().pack(side='top', fill='both', expand=1)

        self.toolbar_n_star = NavigationToolbar2Tk(self.canvas_n_star, self.frame_graph_n_star)
        self.canvas_n_star._tkcanvas.pack(side='top', fill='both', expand=1)


        # nb_aggregates
        self.frame_graph_n_aggregate = tk.Frame(self.frame_graphs)
        self.frame_graph_n_aggregate.grid(row=0, column=3)

        self.figure_n_aggregate= plt.Figure(figsize=(self.fig_size_x, self.fig_size_y), dpi=self.fig_dpi)
        self.ax_n_aggregate = self.figure_n_aggregate.add_subplot(111)

        plt.subplots_adjust(hspace=0)
        self.figure_n_aggregate.set_tight_layout(True)
        self.canvas_n_aggregate = FigureCanvasTkAgg(self.figure_n_aggregate, master=self.frame_graph_n_aggregate)
        # self.canvas.get_tk_widget().pack(side='top', fill='both', expand=1)

        self.toolbar_n_aggregate = NavigationToolbar2Tk(self.canvas_n_aggregate, self.frame_graph_n_aggregate)
        self.canvas_n_aggregate._tkcanvas.pack(side='top', fill='both', expand=1)

        # Mean size
        self.frame_graph_mean_size = tk.Frame(self.frame_graphs)
        self.frame_graph_mean_size.grid(row=1, column=0)

        self.figure_mean_size= plt.Figure(figsize=(self.fig_size_x, self.fig_size_y), dpi=self.fig_dpi)
        self.ax_mean_size = self.figure_mean_size.add_subplot(111)

        plt.subplots_adjust(hspace=0)
        self.figure_mean_size.set_tight_layout(True)
        self.canvas_mean_size = FigureCanvasTkAgg(self.figure_mean_size, master=self.frame_graph_mean_size)
        # self.canvas.get_tk_widget().pack(side='top', fill='both', expand=1)

        self.toolbar_mean_size = NavigationToolbar2Tk(self.canvas_mean_size, self.frame_graph_mean_size)
        self.canvas_mean_size._tkcanvas.pack(side='top', fill='both', expand=1)

        # Concentration
        self.frame_graph_concentration = tk.Frame(self.frame_graphs)
        self.frame_graph_concentration.grid(row=1, column=1)

        self.figure_monomer_attachment= plt.Figure(figsize=(self.fig_size_x, self.fig_size_y), dpi=self.fig_dpi)
        self.ax_n_monomer_attachment = self.figure_monomer_attachment.add_subplot(111)

        plt.subplots_adjust(hspace=0)
        self.figure_monomer_attachment.set_tight_layout(True)
        self.canvas_attachment_monomer = FigureCanvasTkAgg(self.figure_monomer_attachment, master=self.frame_graph_concentration)
        # self.canvas.get_tk_widget().pack(side='top', fill='both', expand=1)

        self.toolbar_monomer_attachment = NavigationToolbar2Tk(self.canvas_attachment_monomer, self.frame_graph_concentration)
        self.canvas_attachment_monomer._tkcanvas.pack(side='top', fill='both', expand=1)

    def create_menu(self):
        self.menu_system = tk.Menu(self.root)

        # FILE#############
        self.menu_file = tk.Menu(self.menu_system, tearoff=0)

        self.menu_file.add_command(label='Import ini file', underline=1, accelerator="Ctrl+n", command=self.import_ini_file)
        #self.master.bind_all("<Control-o>", self.askOpenSPC_file)

        self.menu_file.add_command(label='Export ini file', underline=1, accelerator="Ctrl+o", command=self.export_ini_file)
        #self.master.bind_all("<Control-o>", self.askOpenSPC_file)

        self.menu_system.add_cascade(label="File", menu=self.menu_file)

        self.root.config(menu=self.menu_system)

    def create_string_value_dict(self):
        self.sv_input_dict["time_step_ns"] = self.time_step_ns_sv
        self.sv_input_dict["end_time_calculation_ms"] = self.end_time_calculation_ms
        self.sv_input_dict["monitoring_time_step_micros"] = self.monitor_time_micros_sv
        self.sv_input_dict["temperature"] = self.temperature_K_sv
        self.sv_input_dict["v_solvent_ml"] = self.volume_solution_mL_sv
        self.sv_input_dict["eta"] = self.viscosity_sv
        self.sv_input_dict["p_compacity"] = self.packing_density_p_sv
        self.sv_input_dict["r_nm"] = self.molecule_radius_nm_sv
        self.sv_input_dict["r_hydro_nm"] = self.molecule_hydro_radius_nm_sv

        self.sv_input_dict["gamma_sl"] = self.molecule_tension_surface_sv
        self.sv_input_dict["cini"] = self.monomer_concentration_M_sv

        self.sv_input_dict["ce_m"] = self.molecule_solubility_Ce_sv
        self.iv_input_dict["is_use_turnbull"] = self.is_use_Turnbull_iv

        self.sv_output_t0_dict["S"] = self.S_t0_sv
        self.sv_output_t_dict["S"] = self.S_t_sv
        self.sv_output_t0_dict["J"] = self.J_t0_sv
        self.sv_output_t_dict["J"] = self.J_t_sv
        self.sv_output_t0_dict["C"] = self.C_t0_sv
        self.sv_output_t_dict["C"] = self.C_t_sv
        self.sv_output_t0_dict["nstar"] = self.nstar_t0_sv
        self.sv_output_t_dict["nstar"] = self.nstar_t0_sv
        self.sv_output_t0_dict["Wnstar"] = self.Wnstar_t0_sv
        self.sv_output_t_dict["Wnstar"] = self.Wnstar_t_sv
        self.sv_output_t0_dict["theta"] = self.theta_t0_sv
        self.sv_output_t_dict["theta"] = self.theta_t_sv
        self.sv_output_t0_dict["monomer_attachment_ns"] = self.attach_time_t0_sv
        self.sv_output_t_dict["monomer_attachment_ns"] = self.attach_time_t_sv

    def get_data_from_ui(self):
        self.core.param_dict = {}
        for key in self.sv_input_dict:
            self.core.param_dict[key] = float(self.sv_input_dict[key].get())
        for key in self.iv_input_dict:
            self.core.param_dict[key] = int(self.iv_input_dict[key].get())

    def set_data_to_ui_from_ini(self):
        for key in self.sv_input_dict:
            self.sv_input_dict[key].set(str(self.core.param_dict[key]))
        for key in self.iv_input_dict:
            self.iv_input_dict[key].set(int(self.core.param_dict[key]))

    def update_t0_value(self):
        self.get_data_from_ui()
        self.get_and_display_output_at_t0()

    def get_and_display_output_at_t0(self):
        self.core.get_output_val_at_t0()
        for key in self.core.result_t0_dict:
            self.sv_output_t0_dict[key].set("{:.2e}".format(self.core.result_t0_dict[key]))



    def plot_intermediate_data(self):
        self.core.create_dict_result()
        self.plot_data(self.core.result_dict)
        if self.core.is_calculating:
            self.root.after(self.calculation_polling_delay, self.plot_intermediate_data)

    def recalculate(self):
        self.get_data_from_ui()
        self.core.set_param(self.core.param_dict)
        # self.thread_core_calculation.start()
        self.core.launch_calculation()
        time.sleep(0.5)
        self.plot_intermediate_data()

        # self.core.plot_data()

    def core_calculation(self):
        self.core.calculate()
        while self.core.is_calculating:
            time.sleep(self.calculation_polling_delay)
        self.core.plot_data()

    def plot_data(self, result_dict):
        t = result_dict["t"]
        S = result_dict["S"]
        self.ax_S.clear()
        self.ax_S.plot(t, S)
        self.canvas_S.draw()
        self.ax_S.set_xlabel("time (µs)")
        self.ax_S.set_ylabel("Super Saturation S")


        self.ax_J.clear()
        self.ax_J.plot(t, result_dict["J"])
        self.ax_J.set_ylabel("Nucleation rate in ns/Volume")
        self.ax_J.set_xlabel("time (µs)")
        self.canvas_J.draw()

        self.ax_n_star.clear()
        self.ax_n_star.plot(t, result_dict["n_star"])
        self.ax_n_star.set_ylabel("Critical Nucleus size n*")
        self.ax_n_star.set_xlabel("time (µs)")
        self.canvas_n_star.draw()

        self.ax_n_aggregate.clear()
        self.ax_n_aggregate.plot(t, result_dict["nb_aggregates_array"])
        self.ax_n_aggregate.set_ylabel("Ng of aggregate ")
        self.ax_n_aggregate.set_xlabel("time (µs)")
        self.canvas_n_aggregate.draw()

        self.ax_mean_size.clear()
        self.ax_mean_size.plot(t, result_dict["mean_size"])
        self.ax_mean_size.set_ylabel("Mean size (in nb of molecule)")
        self.ax_mean_size.set_xlabel("time (µs)")
        self.canvas_mean_size.draw()

        self.ax_n_monomer_attachment.clear()
        self.ax_n_monomer_attachment.plot(t, result_dict["monomer_attachment_ns"])
        self.ax_n_monomer_attachment.set_ylabel("Monomer attachment time ns)")
        self.ax_n_monomer_attachment.set_xlabel("time (µs)")
        self.canvas_attachment_monomer.draw()


    def stop(self):
        self.core.stop_calculation()

    def import_ini_file(self):
        file_path = filedialog.askopenfilename(title="Open ini File")
        if file_path == None or file_path == '':
            return None
        self.core.import_ini_file(file_path)
        self.set_data_to_ui_from_ini()
        self.update_t0_value()

    def export_ini_file(self):
        file_path = filedialog.asksaveasfile(title="Export ini File")
        if file_path == None or file_path == '':
            return None
        self.core.export_ini_file(file_path.name)

