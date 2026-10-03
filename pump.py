import numpy as np
import pandas as pd
from constants import(
    pump_steps,
    g,
    pump_eff,
    pump_max_v,)

curve_data = pd.read_csv("pump_curves.csv")




class Pump:
    
    
    
    
    def __init__(self, model):
        
        self.curve_Q = curve_data["Q (m^3/s)"]
        self.curve_H = curve_data["H (m)"]
        self.curve_eta = curve_data["eta (0-1)"]
        self.curve_P = curve_data["P (W)"]             # ELECTRICAL POWER
        
        self.H0 = 10                                   # PUMP SHUTOFF HEAD
        self.pump_k = 0.1                              # PUMP CURVE CONSTANT
        
        self.integral = 0
        self.error = 0
        
        self.pump_v_control = 0
        
        self.model = model
        
        
        
        
    def pump_pid_control(self, dt):
        
        """
        PID controller used to regulate pump velocity, allowing for a stable mass
        flow rate to be reached. This function calculates proportional, integral, 
        and derivative terms which change depending on how large the difference (error) 
        between current and the aim mass flow rate is and subsequently outputs a 
        limited pump velocity.

        Parameters
        ----------
        dt : float
            Time step duration in seconds (dt = 0.1s)

        Returns
        -------
        u_t : float
            Pump velocity control 
        """
        
        Kp = 0.0783
        Ki = 0.04
        Kd = 0.058
        
        m_dot_aim = 1.5 
        
        e_tau = m_dot_aim - self.model.m_dot
        
        self.integral += e_tau * dt
        
        derivative = (e_tau - self.error)/dt
        self.error = e_tau
        
        u_t = max(0.5, min(pump_max_v, (Kp * e_tau) + (Ki * self.integral) + (Kd * derivative)))
        
        self.pump_v_control = u_t
        
        return u_t
    
    
    
    
    def press_diff(self, m_dot, rho):
        
        """
        This function calculates the change in pressure due to the pump for the particular
        mass flow rate and density of the fluid at the time it is called. The pump curves used
        are scaled using the respective affinity laws with fluid velocity and the pump head is
        determined by interpolating the flow rate/pump head curve at the current volumetric flow rate.

        Parameters
        ----------
        m_dot : float
            Mass flow rate of fluid (kg/s)
        rho : float
            Fluid density (kg/m^3)

        Returns
        -------
        Δp : float
            Pressure change caused by pump (Pa)
        """
    
        Q = m_dot/rho
        
        Q_v = self.curve_Q * self.pump_v_control
        H_v = self.curve_H * (self.pump_v_control**2)
        
        H = np.interp(Q, Q_v, H_v)
    
        # CHANGE IN PRESSURE
        Δp = rho * g * H
        
        return Δp
    



    def pump_pwr(self, m_dot, rho): 
        
        """
        Calculates the hydraulic, shaft and electrical power of the pump for 
        the particular mass flow rate and density of the fluid at the time it is called.
        Volumetric flow rate is scaled using velocity affinity law and then multiplied by
        pressure change to calculate hydraulic power. Shaft power is calculated by dividing this 
        value by pump efficiency (eta). Electrical power is calculated by dividing shaft power by
        pump efficiency and scaling with pump affinity laws.

        Parameters
        ----------
        m_dot : float
            Mass flow rate of fluid (kg/s)
        rho : float
            Fluid density (kg/m^3)

        Returns
        -------
        hydraulic_pwr : float
            Hydraulic Power - power delivered to fluid (W)
        shaft_pwr : float
            Shaft power - mechanical power delivered to the pump shaft (W)
        elec_pwr : float
            Electrical power - power used by the motor (W)
        """
        
        Q = m_dot/rho                    # VOLUMETRIC FLOW RATE
        Δp = self.press_diff(m_dot, rho)
        
        Q_v = self.curve_Q * self.pump_v_control
        
        eta = np.interp(Q, Q_v, self.curve_eta)     # PUMP HYDRAULIC EFFICIENCY
        
        # HYDRAULIC POWER (POWER DELIVERED TO MECHANICAL SYSTEM)
        hydraulic_pwr = Q * Δp
        
        # SHAFT POWER (POWER DELIVERED TO FLUID)
        shaft_pwr = hydraulic_pwr/eta
        
        # ELECTRICAL POWER (POWER DELIVERED TO MOTOR (SCALED BY AFFFINITY LAWS))
        elec_pwr = (shaft_pwr/pump_eff) * (self.pump_v_control)**3      # 'pump_eff' REPRESENTS ELECTRICAL EFFICIENCY
        
        return hydraulic_pwr, shaft_pwr, elec_pwr
    
    
    
    
    def pwr_heating(self, m_dot, rho, dt):
        
        """
        Calculates the heat gained by the fluid due to pump losses. The pump’s
        electrical power input and hydraulic power are initially called from 'pump_pwr'. 
        The difference between these two values represents the mechanical losses in 
        the pump which are transformed into heat gain. This heat is distributed across the 
        pump segments and added to the fluid's current temperature at said segment index. 
        An exponential taper is applied to this heat gain to maintain numerical stability 
        at higher temperatures.

        Parameters
        ----------
        m_dot : float
            Mass flow rate of fluid (kg/s)
        rho : float
            Fluid density (kg/m^3)
        dt : float
            Time step duration in seconds (dt = 0.1s)

        Returns
        -------
        Q: float
            Heat generated by pump losses (W)
        """       
        
        T_aim = 303.15                  # 30°C
        mean_T = np.mean(self.model.fluid_T)
        k = 3                           # TAPER CONSTANT (ARBITRARILY TUNED)
        
        cp = self.model.cp_T
        
        hydraulic_pwr, shaft_pwr, elec_pwr = self.pump_pwr(m_dot, rho)
        
        Q = elec_pwr - hydraulic_pwr
        Q_seg = Q/len(pump_steps)
        
        pump_taper = max(0.05, min(0.3, np.exp(-(mean_T - T_aim)/k)))
        
        for i in pump_steps:
            
            self.model.fluid_T[i] += pump_taper * Q_seg/(m_dot * cp[i]) * dt
        
        return Q