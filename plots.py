import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from constants import g, pump_eff

curve_data = pd.read_csv("pump_curves.csv")




class Plot:
    
    
    
    
    def __init__(self, model, simulation, heater, cooler, pump):
        
        
        
        
        self.model = model
        self.simulation = simulation
        self.heater = heater
        self.cooler = cooler
        self.pump = pump 




    def plot_all(self):
        
        self.heatmap()
        self.fluid_temp()
        self.mean_fluid_temp()
        self.heater_temp()
        self.cooler_temp()
        self.pump_v()
        self.pump_curve()
        self.pump_pwr_curve()
        self.pump_pwr()
        self.mass_flow_rate()
        self.fluid_v()
        self.head_loss()
        self.density()
        
        
        
        
    def heatmap(self):
    
        T_hist = np.array(self.model.T_history)
        T0 = T_hist[0]
        ΔT = T_hist - T0
        
        plt.figure(1)
        plt.imshow(T_hist, origin = "lower", cmap = "inferno", aspect = "auto")
        colour_bar = plt.colorbar()
        colour_bar.set_label("Temperature (K)")
        plt.title("Heat map of loop")
        plt.xlabel("Segment of Loop")
        plt.ylabel("Time Step")
        plt.tight_layout()
        
        plt.figure(2)
        plt.imshow(ΔT, origin = "lower", cmap = "coolwarm", aspect = "auto")
        colour_bar = plt.colorbar()
        colour_bar.set_label("Temperature (K)")
        plt.title("ΔT Heat map of loop")
        plt.xlabel("Segment of Loop")
        plt.ylabel("Time Step")
        plt.tight_layout()
        plt.show()
        
        
        
        
    def fluid_temp(self):
        
        plt.plot(self.model.fluid_T, '-r')
        plt.title("Fluid Temp Per Pipe Segment")
        plt.xlabel("Pipe Segment")
        plt.ylabel("Temperature (K)")
        plt.minorticks_on()
        plt.grid(True)
        plt.show()
        
        
        
        
    def mean_fluid_temp(self):
        
        plt.plot(self.simulation.t_arr, self.model.mean_T_history, 'r')
        plt.title("Mean Fluid Temp vs Time")
        plt.xlabel("Time (S)")
        plt.ylabel("Temperature (K)")
        plt.minorticks_on()
        plt.grid(True)
        plt.show()
        
        
        
        
    def heater_temp(self):
        
        plt.plot(self.simulation.t_arr, self.simulation.heater_core_T_arr, '-r', label = "Core Temp")
        plt.plot(self.simulation.t_arr, self.simulation.heater_ins_T_arr, '-b', label = "Insulation Temp")
        plt.plot(self.simulation.t_arr, self.simulation.heater_outer_T_arr, '-g', label = "Outer Surface Temp")
        plt.title("Heater Temp vs Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Temperature (K)")
        plt.minorticks_on()
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()
        
        
        
        
    def cooler_temp(self):
        
        plt.plot(self.simulation.t_arr, self.simulation.cooler_T_arr, '-b')
        plt.title("Cooler Temp vs Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Temperature (K)")
        plt.minorticks_on()
        plt.grid(True)
        plt.show()
        
        
        
        
    def pump_v(self):
        
        plt.plot(self.simulation.t_arr, self.simulation.pump_v_arr, '-m')
        plt.title("Pump Velocity Control vs Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Pump Velocity control")
        plt.minorticks_on()
        plt.grid(True)
        plt.show()
        
        
        
        
    def pump_curve(self):
    
        plt.plot(curve_data["Q (m^3/s)"], curve_data["H (m)"], '-r', label = "Performance Curve" )
        plt.plot(curve_data["Q (m^3/s)"], curve_data["eta (%)"], '-b', label = "Efficiency Curve" )
        plt.title("Pump Performance/Efficiency Curve")
        plt.xlabel("Volumetric Flow Rate ($m^3/s$)")
        plt.ylabel("Pump Head (m)")
        plt.minorticks_on()
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()
    
    
    
    
    def pump_pwr_curve(self):
    
        Q = curve_data["Q (m^3/s)"]
        H = curve_data["H (m)"]
        eta = curve_data["eta (%)"]/100
        
        rho = np.mean(self.model.rho_T)
        
        hydraulic = Q * (rho * g * H)
        shaft = hydraulic/eta
        elec = shaft/pump_eff
    
        plt.plot(curve_data["Q (m^3/s)"], hydraulic, '-r', label = "Hydraulic Power Curve" )
        plt.plot(curve_data["Q (m^3/s)"], shaft, '-b', label = "Shaft Power Curve" )
        plt.plot(curve_data["Q (m^3/s)"], elec, '-g', label = "Electrical Power Curve" )
        plt.title("Pump Power Curves")
        plt.xlabel("Volumetric Flow Rate ($m^3/s$)")
        plt.ylabel("Power (W)")
        plt.minorticks_on()
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()
    
    
    
    
    def pump_pwr(self):
    
        plt.plot(self.simulation.t_arr, self.simulation.hyd_pwr_arr, '-r', label = "Hydraulic Power")
        plt.plot(self.simulation.t_arr, self.simulation.shaft_pwr_arr, '-b', label = "Shaft Power")
        plt.plot(self.simulation.t_arr, self.simulation.elec_pwr_arr, '-g', label = "Electrical Power")
        plt.title("Pump Power vs Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Power (W)")
        plt.minorticks_on()
        plt.grid(True)
        plt.legend()
        plt.tight_layout()
        plt.show()
        
        
        
        
    def mass_flow_rate(self):
        
        plt.plot(self.simulation.t_arr, self.simulation.m_dot_arr, '-k')
        plt.title("Mass Flow Rate vs Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Mass Flow Rate ($kgs^{-1}$)")
        plt.minorticks_on()
        plt.grid(True)
        plt.show()
        
    
    
    
    def fluid_v(self):
    
        plt.plot(self.simulation.t_arr, self.simulation.fluid_v_arr, '-m')
        plt.title("Fluid Velocity vs Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Velocity ($ms^{-1}$)")
        plt.minorticks_on()
        plt.grid(True)
        plt.show() 
        
        
        
        
    def head_loss(self):
        
        plt.plot(self.simulation.t_arr, self.simulation.hf_arr, '-c')
        plt.title("Head Loss vs Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Head Loss (m)")
        plt.minorticks_on()
        plt.grid(True)
        plt.show()
        
        
        
        
    def density(self):
        
        plt.plot(self.simulation.t_arr, self.simulation.rho_arr, '-c')
        plt.title("Fluid Density vs Time")
        plt.xlabel("Time (s)")
        plt.ylabel("Density ($kgm^{-3}$)")
        plt.minorticks_on()
        plt.grid(True)
        plt.show()