import numpy as np
from constants import(
    steps,
    heater_steps,
    stef_boltz,
    
    heater_max_pwr,
    
    core_L,
    core_r,
    core_rho,
    core_cp,
    core_kappa,
    
    ins_L,
    ins_r,
    ins_rho,
    ins_cp,
    ins_kappa,
    
    outer_L,
    outer_r,
    outer_rho,
    outer_cp,
    outer_epsilon,)




class Heater:
    
    
    
    
    def __init__(self, model):
        
        self.core_T = 293.15
        self.ins_T = 293.15
        self.outer_T = 293.15
        
        self.integral = 0
        self.error = 0
        
        self.heater_pwr_control = 0
        
        self.outer_A = np.pi * ((2 * outer_r * outer_L) + 2 * outer_r**2)   # HEATER CYLINDER + END CAPS
        
        core_V = np.pi * (core_r**2) * core_L
        self.core_m = core_rho * core_V
        
        ins_V = np.pi * ((ins_r**2) - (core_r**2)) * ins_L
        self.ins_m = ins_rho * ins_V
        
        outer_V = np.pi * ((outer_r**2) - (ins_r**2)) * outer_L
        self.outer_m = outer_rho * outer_V
        
        self.model = model
        
    
        
        
    def heater_pid_control(self, dt):
        
        """
        PID controller used to regulate heater power. This function calculates 
        proportional, integral, and derivative terms which change depending on 
        how large the difference (error) between current and the aim fluid temperature
        is and subsequently outputs a limited heater power.

        Parameters
        ----------
        dt : float
            Time step duration in seconds (dt = 0.1s)

        Returns
        -------
        u_t : float
            Heater power control 
        """
        
        Kp = 0.00005
        Ki = 0.000225
        Kd = 0.0005
        
        T_aim = 303.15              # 30°C
        T = np.mean(self.model.fluid_T)       # MEAN OF FLUID IN CONTACT WITH HEATER
        
        e_tau = T_aim - T
        
        derivative = (e_tau - self.error)/dt
        self.error = e_tau
        
        u_0 = (Kp * e_tau) + (Ki * self.integral) + (Kd * derivative)
        
        # ONLY ADD TO INTEGRAL TERM WHEN POWER CONTROL IS NOT SATURATED AT 1 OR 0
        if not ((u_0 >= 1 and e_tau > 0) or (u_0 <= 0 and e_tau < 0)):
            self.integral += e_tau * dt
        
        u_t = max(0, min(1, u_0))
        
        self.heater_pwr_control = u_t
        
        return u_t
    
        
        
            
    def update_heater_T(self, dt):
        
        """
        This function updates the temperatures of the three sections included inside
        the heater component (core, insulating layer, and outer layer). The core gains heat 
        via the electrical power of the heater itself. This heat is then transferred into 
        the insulating layer and then to the outer layer via cylindrical radial conduction 
        in both instances. This heat is then eventually lost via convection and radiation from 
        the outer layer to the fluid in the system. All temperature updates are performed using 
        energy balances for each layer.

        Parameters
        ----------
        dt : float
            Time step duration in seconds (dt = 0.1s)
        """
                   
        h = np.mean(self.model.h_T[heater_steps])           # MEAN CONVECTION COEFFICIENT OF FLUID IN CONTACT WITH HEATER
        T = np.mean(self.model.fluid_T[heater_steps])       # MEAN TEMP OF FLUID IN CONTACT WITH HEATER
            
        # CORE HEATER HEAT LOSS TO INSULATING LAYER (CYLINDRICAL RADIAL CONDUCTION)
        Q_core_cond = (2 * np.pi * core_kappa * core_L * (self.core_T - self.ins_T))/np.log(ins_r/core_r)
        
        # INSULATING LAYER HEAT LOSS TO OUTER LAYER (CYLINDRICAL RADIAL CONDUCTION)
        Q_ins_cond = (2 * np.pi * ins_kappa * ins_L * (self.ins_T - self.outer_T))/np.log(outer_r/ins_r)
            
        # OUTER HEATER HEAT LOSS TO FLUID (CONVECTION)
        Q_out_conv = h * self.outer_A * (self.outer_T - T) 
        
        # OUTER HEATER HEAT LOSS TO FLUID (RADIATION)
        Q_out_rad = stef_boltz * outer_epsilon * self.outer_A * ((self.outer_T**4) - (T**4))
        
        heater_pwr = heater_max_pwr * self.heater_pwr_control
        
        self.core_T += (heater_pwr - Q_core_cond)/(self.core_m * core_cp) * dt
        self.ins_T += (Q_core_cond - Q_ins_cond)/(self.ins_m * ins_cp) * dt
        self.outer_T += (Q_ins_cond - Q_out_conv - Q_out_rad)/(self.outer_m * outer_cp) * dt
        
        

            
    def Q_heating(self, dt):
        
        """
        Calculates the heat gained by the fluid (segment based) as it travels across 
        the immersed heater component. Heat is added to the fluid via three processes. 
        Heat is added by the active powered component of the heater and is additionally 
        gained through convection and radiation between the fluid and the heater wall. 
        An index based Gaussian taper is applied to the heating of the fluid. This taper 
        represents that the heater will deliver more heat to the central region of the heater 
        and less heat toward the edges. The taper models the distribution of heating 
        effectiveness along the heater length, preventing the fluid from receiving excessive 
        and uneven heating across segments.

        Parameters
        ----------
        dt : float
            Time step duration in seconds (dt = 0.1s)

        Returns
        -------
        Q_heat_conv : Array
            Heat gain by the fluid via convection from the heater surface (W)
        Q_heat_rad : Array
            Heat gain by the fluid via radiation from the heater surface (W)
        """
        
        m_dot = self.model.m_dot
        
        h = self.model.h_T
        T = self.model.fluid_T
        cp = self.model.cp_T
        
        heater_pwr = self.heater_pwr_control * heater_max_pwr

        Q_seg = heater_pwr/len(heater_steps)
        
        seg_A = self.outer_A/len(heater_steps)
        
        Q_heat_conv = np.zeros(steps)
        Q_heat_rad = np.zeros(steps)
        Q_gain = np.zeros(steps)
        
        mu = len(heater_steps)/2
        sigma = len(heater_steps)/3
        
        for x, i in enumerate(heater_steps):
            
            # GAUSSIAN HEAT DISTRIBUTION 
            heater_taper = np.exp(-((x - mu)/sigma)**2)
            
            # CONVECTION FROM HEATER TO FLUID
            Q_heat_conv[i] = heater_taper * h[i] * seg_A * (self.outer_T - T[i])
            
            # RADIATION FROM HEATER TO FLUID
            T_rad = np.clip(self.outer_T, 293.15, 393.15)            # PREVENTS RUNAWAY TEMPS
            Q_heat_rad[i] = heater_taper * stef_boltz * outer_epsilon * seg_A * ((T_rad**4) - (T[i]**4))
            
            # TOTAL HEAT GAIN
            Q_gain[i] = Q_heat_conv[i] + Q_heat_rad[i] + Q_seg
            
        for i in heater_steps:
            
            self.model.fluid_T[i] += (Q_gain[i])/(m_dot * cp[i]) * dt
        
        return Q_heat_conv, Q_heat_rad
