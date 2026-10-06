import numpy as np
from constants import(
    steps,
    cooler_steps,
    stef_boltz,
    env_T,
    h_air,
    
    pipe_d,
    pipe_t,
     
    cooler_max_pwr,
    cooler_eff,
    cooler_L,
    cooler_rho,
    cooler_cp,
    cooler_epsilon)
    
    
    

class Cooler:
    
    
    
    
    def __init__(self, model):
        
        self.cooler_T = 293.15
        
        self.integral = 0
        self.error = 0
        
        self.cooler_control = 0
        
        self.Q_cool_conv = np.zeros(steps)

        r_inn = pipe_d/2
        r_out = r_inn + pipe_t
       
        cooler_V = np.pi * ((r_out**2) - (r_inn**2)) * cooler_L
        self.cooler_m = cooler_rho * cooler_V
        
        self.cooler_inn_A = 2 * np.pi * r_inn * cooler_L
        self.cooler_out_A = 2 * np.pi * r_out * cooler_L
        
        self.model = model
        
        
        
        
    def cooler_pid_control(self, dt):
        
        """
        PID controller used to regulate cooler power. This function calculates 
        proportional, integral, and derivative terms which change depending on 
        how large the difference (error) between current and the aim fluid temperature
        is and subsequently outputs a limited cooler power. Anti windup is implemented
        into the PID controller to ensure that if the PID saturates at 0 or 1 
        (and the error continues to push in said direction) that the integral term 
        will not build constantly. This prevents massive integral term overshoots
        which cause the PID to stay saturated for the whole simulation run.

        Parameters
        ----------
        dt : float
            Time step duration in seconds (dt = 0.1s)

        Returns
        -------
        u_t : float
            Cooler power control 
        """
         
        Kp = 0.00005
        Ki = 0.000225
        Kd = 0.0005
        
        T_aim = 299.15              # 27°C
        T = np.mean(self.model.fluid_T)       # MEAN OF FLUID IN CONTACT WITH COOLER
         
        e_tau = T - T_aim
         
        derivative = (e_tau - self.error)/dt
        self.error = e_tau
         
        u_0 = (Kp * e_tau) + (Ki * self.integral) + (Kd * derivative)
        
        # ONLY ADD TO INTEGRAL TERM WHEN POWER CONTROL IS NOT SATURATED AT 1 OR 0
        if not ((u_0 >= 1 and e_tau > 0) or (u_0 <= 0 and e_tau < 0)):
            self.integral += e_tau * dt
        
        u_t = max(0, min(1, u_0))
         
        self.cooler_control = u_t
         
        return u_t
     
        
     
        
    def Q_cooling(self, dt):
        
        """
        Calculates the heat lost by the fluid (segment based) as it travels through the cooler 
        component. Heat is removed from the fluid via two processes. Heat is removed 
        via the active powered component of the cooler and is additionally lost through convection 
        between the fluid and the cooler wall. A linear, index based taper is applied to the 
        cooling of the fluid. This taper represents that the fluid at the start of the cooler 
        will be far hotter than the fluid at the end of the cooler segments. The taper models 
        the reduction in cooling effectiveness as the fluid becomes progressively colder as it 
        passes through the length of the cooler. This prevents the fluid from being cooled excessively.
        

        Parameters
        ----------
        dt : float
            Time step duration in seconds (dt = 0.1s)

        Returns
        -------
        Q_cool_conv : Array
            Heat gain by the cooler wall via convection from fluid (W)
        """
        
        m_dot = self.model.m_dot
        
        h = self.model.h_T              
        T = self.model.fluid_T
        cp = self.model.cp_T
        
        cooler_pwr = self.cooler_control * cooler_max_pwr
        Q_seg = (cooler_pwr/len(cooler_steps))
        
        seg_A = self.cooler_inn_A/len(cooler_steps)
        
        Q_cool_conv = np.zeros(steps)
        Q_loss = np.zeros(steps)
        
        for index, i in enumerate(cooler_steps):          # SEGMENTS 100-120
        
            cooler_taper = cooler_eff * (1 - index/len(cooler_steps))     # COOLER FLUID SEGMENTS --> LESS HEAT TRANSFER
            
            # FLUID HEAT LOSS VIA COOLER (CONVECTION)
            Q_cool_conv[i] = cooler_taper * h[i] * seg_A * (T[i] - self.cooler_T)
            self.Q_cool_conv = Q_cool_conv
            
            # TOTAL HEAT LOSS (CONVECTION + COOLER POWER)
            Q_loss[i] = Q_cool_conv[i] + Q_seg
        
        for i in cooler_steps:
            
            self.model.fluid_T[i] -= (Q_loss[i])/(m_dot * cp[i]) * dt

        return Q_cool_conv

        


    def update_cooler_T(self, dt): 
        
        """
        This function updates the temperature of the cooler component. The cooler initially 
        gains heat from the fluid through the convection term called from 'Q_cooling'.
        The cooler then additionally loses heat via two processes. The electrical power of 
        the cooler itself and via convection and radiation to the environment.

        Parameters
        ----------
        dt : float
            Time step duration in seconds (dt = 0.1s)
        """
        
        # HEAT GAIN VIA FLUID
        
        Q_fluid_gain = 0
        
        for i in cooler_steps:
            Q_fluid_gain += self.Q_cool_conv[i]
        
        # PASSIVE HEAT LOSS
        
        # COOLER HEAT LOSS TO ENVIRONMENT (CONVECTION)
        Q_env_conv = h_air * self.cooler_out_A * (self.cooler_T - env_T)
    
        # COOLER HEAT LOSS TO ENVIRONMENT (RADIATION)
        T_rad = np.clip(self.cooler_T, 293.15, 393.15)
        Q_env_rad = stef_boltz * cooler_epsilon * self.cooler_out_A * ((T_rad**4) - (env_T**4))
    
        Q_loss_pass = Q_env_conv + Q_env_rad
        
        # ACTIVE HEAT LOSS THROUGH EXTERNAL COMPONENT
        Q_loss_act = self.cooler_control * cooler_max_pwr
        
        Q_total_loss = Q_loss_pass + Q_loss_act
        Q_total = Q_fluid_gain - Q_total_loss
    
        self.cooler_T += (Q_total)/(self.cooler_m * cooler_cp) * dt
