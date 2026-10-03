import numpy as np

from heater import Heater
from cooler import Cooler
from pump import Pump

from constants import(
    steps,
    cooler_steps,
    stef_boltz,
    g,
    env_T,
    pipe_L,
    pipe_d,
    pipe_t,
    pipe_rough,
    rho_pipe,
    cp_pipe,
    kappa_pipe,
    pipe_epsilon,
    h_air,)




class Model:
    
    
    
    
    def __init__(self):
        
        self.fluid_T = np.full(steps, 293.15)       # STARTS AT 20°C
        self.T_history = []
        self.mean_T_history = []
        self.m_dot = 0.5                                # MASS FLOW RATE
        self.cooler_T = env_T
        
        self.pipe_inner_T = self.fluid_T.copy()         # INNER SEGMENT STARTS AT FLUID TEMP
        self.pipe_outer_T = np.full(steps, env_T)       # OUTER SEGMENT STARTS AT ENVIRONMENT TEMP
        
        self.r_inn = pipe_d/2                           
        self.r_out = self.r_inn + pipe_t
        
        self.seg_L = pipe_L/steps                            # LENGTH OF PIPE SEGMENT
        
        self.pipe_inn_A = 2 * np.pi * self.r_inn * self.seg_L
        self.pipe_out_A = 2 * np.pi * self.r_out * self.seg_L
        
        self.pipe_cross_A = ((np.pi) * (pipe_d)**2)/4
        
        self.rho_T = np.full(steps, 9.982e2)           # DENSITY OF WATER AT 20°C
        self.eta_T = np.full(steps, 1.002e-3)          # VISCOCITY OF WATER AT 20°C                   
        self.cp_T = np.full(steps, 4.186e3)            # SPECIFIC HEAT OF WATER AT 20°C
        self.epsilon_T = np.full(steps, 0.97)          # EMISSIVITY OF WATER AT 20°C
        self.kappa_T = np.full(steps, 0.6)             # THERMAL CONDUCTIVITY OF WATER AT 20°C

        self.heater = Heater(self)
        self.cooler = Cooler(self)
        self.pump = Pump(self)
        
        
        
        
    """
    The following five functions calculate various fluid/pipe relationship properties
    using the outlined methods below:
        
        - Mass per segment: calculated using fluid density and pipe volume.
        - Fluid velocity: calculated using mass flow rate, fluid density and cross sectional pipe area
        - Friction factor: calculated using the Churchill equation which is applicable to both
          laminar and turbulent flow scenarios.
        - Head loss: calculated using the Darcy-Weisbach equation
        - friction pressure difference: calculated using Darcy-Weisbach friction pressure change
          across all pipe segments
        
    These variables are then all integrated into the thermo-mechanical simulation loop and are
    updated at every timestep.
    """
    
    
    
    
    def mass_per_seg(self):
        
        seg_L = (pipe_L/steps)
        seg_V = self.pipe_cross_A * seg_L
        
        m = self.rho_T * seg_V
        return m
    
    
    
    
    def fluid_v(self):
        
        v = np.zeros(steps)
        
        rho = np.mean(self.rho_T)
        
        for i in range(steps):
        
            v[i] = (self.m_dot)/(rho * self.pipe_cross_A)
           
        return v
    
    
    
    
    def fric_factor(self):
        
        Re = np.zeros(steps)
        f = np.zeros(steps)
        v = self.fluid_v()
        
        # REYNOLD'S NUMBER PER SEGMENT
        for i in range(steps):
            
            Re[i] = ((self.rho_T[i] * v[i] * pipe_d)/self.eta_T[i])
            
            # CHURCHILL EQUATION FOR LAMINAR/TURBULENT FLOW FRICTION FACTOR
            a = (2.457 * np.log(1/(((7/Re[i])**0.9) + (0.27 * (pipe_rough/pipe_d)))))**16
            b = (37530/Re[i])**16
            
            f[i] = 8 * (((8/Re[i])**12) + ((a + b)**-1.5))**(1/12)
        
        return f
        
    
    
            
    def head_loss(self):
        
        hf = np.zeros(steps)
        
        f = self.fric_factor()
        v = self.fluid_v()
        
        # DARCY-WEISBACH EQUATION
        for i in range(steps):
            
            hf[i] = (f[i] * self.seg_L * v[i]**2/ (pipe_d * (2 * g)))
        
        return hf
    
    
    
    
    def fric_press_diff(self):
    
        hf = self.head_loss()
        rho = np.mean(self.rho_T)       
    
        # DARCY-WEISBACH FRICTION PRESSURE CHANGE
        Δp_fric = 0.001 * np.sum(hf * rho) * g       # SCALED FOR REALISM
        return Δp_fric
    
    
    
    
    def advection(self, dt):
        
        """
        Models the circulation of fluid and the movement of fluid temperature in 
        the simulated loop.
        
        This function works by determining where a segmented section of fluid came from
        one timestep prior and then moving the previous fluid temperature downstream to allign
        with where this fluid now lies. Linear interpolation is then used to sample the previous 
        temperature distribution.
        
        A diffusion term is added in this function to model the mixing of neighbouring
        segments which each represent water of varying temperatures.
        
        A circular loop and modulo operators are used to ensure that the loop acts periodically 
        regardless of the value of the current index i.

        Parameters
        ----------
        dt : float
            Time step duration in seconds (dt = 0.1s)
        """
        
        k = 0.01 # ARBITRARILY CHOSEN DIFFUSION CONSTANT
        v = self.fluid_v()
        
        dx = v * dt
        frac_seg = dx/self.seg_L
    
        # LIMIT TO PREVENT FLUID FROM TRAVELLING >1 SEGMENT PER TIMESTEP
        frac_seg = np.clip(frac_seg, 0, 0.99)
        
        T_advec = np.full(steps, 293.15)
        x = np.arange(steps)
            
        for i in range(steps):
            
            prev_dx = i - frac_seg[i]
            
            if prev_dx < 0:
                prev_dx += steps
            
            elif prev_dx >= steps:
                prev_dx -= steps
            
            # FINDS TEMP AT LATER SEGMENT (CAN BE A FRACTION)
            T_advec[i] = np.interp(prev_dx, x, self.fluid_T) 
            
            # DIFFUSION (MIXING OF WATER OF DIFFERING TEMPS)
            i_next = (i + 1) % steps
            i_prev = (i - 1) % steps
            
            T_advec[i] += k * (self.fluid_T[i_next] - (2 * (self.fluid_T[i])) + self.fluid_T[i_prev])
            
        self.fluid_T = T_advec
            
        
        
        
    def fric_heating(self):
        
        """
        Calculates the friction based heating generated in each pipe segment which is
        caused by pressure losses and volumetric flow rate. An exponential taper is applied 
        to this value to prevent sustained high levels of friction based heating which would 
        cause numerical instabilities to arise in the system.

        Returns
        -------
        fric: Array
            Array filled with segment based friction heating power (W)
        """
        
        T_aim = 303.15                  # 30°C
        mean_T = np.mean(self.fluid_T)
        k = 3.5                           # TAPER CONSTANT (ARBITRARILY TUNED)
        
        fric = np.zeros(steps)
        
        Δp = self.rho_T * g * self.head_loss()            # CHANGE IN PRESSURE
        Q = self.pipe_cross_A * self.fluid_v()            # VOLUMETRIC FLOW RATE
        
        fric_taper = max(0.05, min(0.3, np.exp(-(mean_T - T_aim)/k)))
        
        for i in range(steps):
            
            fric[i] = fric_taper * Δp[i] * Q[i]
        
        return fric
    
    
    
    
    def fluid_phys(self, dt):
        
        """
        The main function of the simulation. Calls and integrates all functions 
        relating to the changing of the system's fluid temperature/location including: 
            
            - Advection which moves fluid around loop
            - Friction heating
            - Heat gain via heater
            - Heat loss via cooler
            - Heat loss to piping via convection
            - Heat gain via friction from pressure loss
            
        This function also updates the 'self.mean_T_history' array, allowing for
        snapshots of fluid temperature in the whole loop to be taken and plotted at every
        time interval.
        
        Note: 'pipe_phys' is called inside this function and so doesn't need to be called
              independently inside the function 'run_sim' . 

        Parameters
        ----------
        dt : float
            Time step duration in seconds (dt = 0.1s)
        """
        
        self.advection(dt)
        
        fric = self.fric_heating()
        m_seg = self.mass_per_seg()
        
        # HEAT GAIN FROM HEATER
        self.heater.Q_heating(dt)
        
        # HEAT LOSS FROM COOLER
        self.cooler.Q_cooling(dt)
        
        Q_conv_inn = self.pipe_phys(dt)[0]
        
        for i in range(steps):
            
            # HEAT GAIN FROM FRICTION
            self.fluid_T[i] += fric[i]/(m_seg[i] * self.cp_T[i]) * dt
            
            if i not in cooler_steps:
            
                # HEAT LOSS FROM CONVECTION
                self.fluid_T[i] -= Q_conv_inn[i]/(m_seg[i] * self.cp_T[i]) * dt
        
        self.mean_T_history.append(np.mean(self.fluid_T))
        
        # TRACKS FLUID TEMP OVER TIME
        self.T_history.append(self.fluid_T.copy())
        
        
        
        
    def pipe_phys(self, dt):
        
        """
        This function calculates and updates the inner and outer pipe wall temperatures for each
        segment based on four key heat transfer mechanisms and returns the according arrays
        for each respective heat transfer term for use in other functions. It does this by
        calculating the heat transfer between the fluid, inner pipe wall, outer pipe wall and
        environment through the following process:

            - Convection from the fluid to the inner pipe wall
            - Conduction through the pipe from the inner wall to the outer wall
            - Radiation from the outer pipe wall to the environment
            - Convection from the outer pipe wall to the environment

        Parameters
        ----------
        dt : float
            Time step duration in seconds (dt = 0.1s)

        Returns
        -------
        Q_conv_inn : Array
            Heat gain (inner pipe wall) from fluid via convection
        Q_cond : Array
            Heat gain (outer pipe wall) from inner pipe segment via conduction
        Q_rad : Array
            Heat loss (outer pipe wall) to environment via radiation
        Q_conv_out : Array
            Heat loss (outer pipe wall) to environment via convection
        """

        pipe_V = ((self.r_out**2) - (self.r_inn**2)) * np.pi * self.seg_L      # VOLUME OF PIPE SEGMENT
        pipe_m = rho_pipe * pipe_V                                             # MASS OF PIPE SEGMENT
        
        # CONVECTION FROM FLUID TO INNER SEGMENT
        Q_conv_inn = self.h_T * self.pipe_inn_A * (self.fluid_T - self.pipe_inner_T)
        
        # CONDUCTION FROM INNER SEGMENT TO OUTER SEGMENT
        Q_cond = (2 * np.pi * kappa_pipe * self.seg_L * (self.pipe_inner_T - self.pipe_outer_T))/(np.log(self.r_out/self.r_inn))
        
        # RADIATION FROM OUTER SEGMENT TO ENVIRONMENT
        T_rad = np.clip(self.pipe_outer_T, 293.15, 593.15)
        Q_rad = stef_boltz * pipe_epsilon * self.pipe_out_A * ((T_rad**4) - (env_T**4))
        
        # CONVECTION FROM OUTER SEGMENT TO ENVIRONMENT
        Q_conv_out = h_air * self.pipe_out_A * (self.pipe_outer_T - env_T)
        
        for i in range(steps):
            
            if i in cooler_steps:
                continue
        
            self.pipe_inner_T[i] += (Q_conv_inn[i] - Q_cond[i])/(pipe_m * cp_pipe)  * dt
            self.pipe_outer_T[i] += (Q_cond[i] - Q_rad[i] - Q_conv_out[i])/(pipe_m * cp_pipe)  * dt

        return Q_conv_inn, Q_cond, Q_rad, Q_conv_out



    
    """
    The following functions maintain and update the temperature dependent behaviour 
    of several fluid properties. Some variables are modelled using known relationships 
    and others are modelled using empirical data and approximations:
    
        - Density: modelled using the Kell equation.
        - Specific heat capacity: modelled using a polynomial approximation.
        - Viscosity: modelled using the Andrade equation.
        - Emissivity: modelled using a linear approximation.
        - Thermal conductivity: modelled using a polynomial approximation.
        - Convection coefficient: modelled using a combination of equations (outlined further below).

    All variables are updated using the function 'update_variable_T'
    and this is in turn called in the main simulation loop at every timestep.
    """




    def update_variable_T(self):
        
        self.rho_T = self.density_T()
        self.cp_T = self.specific_heat_T()
        self.eta_T = self.viscosity_T()
        self.epsilon_T = self.water_emiss_T()
        self.kappa_T = self.thermal_cond_T()
        self.h_T = self.convec_coeff_T()
        
        
        
        
    # KELL EQUATION
    def density_T(self):
        
        a = -2.8054e-10
        b = 1.0556e-7
        c = -4.6170e-5
        d = 7.9870e-3
        e = 16.945
        f = 9.9984e2
        g = 1.6890e-2
        
        num = ((a * self.fluid_T)**5) + ((b * self.fluid_T)**4) + ((c * self.fluid_T)**3) + ((d * self.fluid_T)**2) + ((e * self.fluid_T)) + f
        den = 1 + (g * self.fluid_T)
        
        rho = num/den
        return rho
        
        
        
        
    # POLYNOMIAL MODEL
    def specific_heat_T(self):
        
        a = 4.2174
        b = 3.4e-2
        c = 1.2e-5
        d = -2.8e-8
        
        cp = (d * (self.fluid_T)**3) + (c * (self.fluid_T)**2) + (b * self.fluid_T) + a
        return cp
        
    
    
    
    # EXTENDED ANDRADE EQUATION
    def viscosity_T(self):
        
        a = 1.856e-11
        b = 4209
        c = 4.527e-2
        d = -3.376e-5
        
        eta = a * np.exp((d * (self.fluid_T)**2) + (c * self.fluid_T) + (b/self.fluid_T))
        return eta

    
    

    # APPROXIMATION
    def water_emiss_T(self):
        
        ε0 = 0.97       # WATER EMISSIVITY AT 20°C
        T0 = 293.15     # 20°C
        k = -2.5e-4     # EMISSIVITY DECREASES SLIGHTLY AS TEMP INCREASES
        
        epsilon = ε0 + k * (self.fluid_T - T0)
        return epsilon
    
    
    
    # POLYNOMIAL APPROXIMATION
    def thermal_cond_T(self):
        
        kappa = 0.561 + 1.85e-3 * (self.fluid_T) -6e-6 * ((self.fluid_T)**2)
        return kappa
    
    
    
    
    def convec_coeff_T(self):
        
        """
        Calculates the temperature dependent convection coefficient
        of the liquid. The function follows the outlined path:
        
            - Reynold's number and Prandtl's number are calculated initially
              using known pipe geometry and fluid variables
              
            - For laminar flow (Re < 2300):
                
                - Graetz's number is used to identify if a thermal layer has been
                  fully developed. This determines what approximation of Nusselt's
                  number is used.
                  
            - For turbulent flow (Re ≥ 2300):
                
                - The Dittus–Boelter correlation is used to determine Nusselt's 
                  number
                  
            - The fluid's convection coefficient is then calculated using the equation
              h = (Nu x κ)/d

        Returns
        -------
        h: Array
            Convection coefficient of the coolant liquid at each pipe segment
        """
        
        Nu = np.zeros(steps)
        
        # REYNOLD'S NUMBER
        Re = (self.rho_T * self.fluid_v() * pipe_d)/self.eta_T
        
        # PRANDTL's NUMBER
        Pr = (self.eta_T * self.cp_T)/self.kappa_T
        
        for i in range(steps):
            
            L = (i + 1) * (pipe_L/steps)
            
            if Re[i] < 2300:            # LAMINAR FLOW
                    
                    # GRAETZ'S NUMBER
                    Gz = (pipe_d * Re[i] * Pr[i])/L
                    
                    if Gz < 1000:       # FULLY DEVELOPED THERMAL LAYER
                        Nu[i]= 3.66
                        
                    else:
                        Nu[i] = 3.66 + (0.0668 * Gz)/(1 + 0.04 * (Gz**(2/3)))
            
            else:                       # TURBULENT FLOW
                # DITTUS-BOELTER EQUATION
                Nu[i] = 0.023 * (Re[i]**0.8) * (Pr[i]**0.4)
        
        h = (Nu * self.kappa_T)/(pipe_d)
        return h