import numpy as np



class Simulation:
    
    
    
    
    def __init__(self, model, heater, cooler, pump):
        
        
        
        
        self.model = model
        self.heater = model.heater
        self.cooler = model.cooler
        self.pump = model.pump
        
        self.press_diff = 0
        self.flow_rate = 0
        
        self.t_count = 0
        
        # MISCELLANEOUS ARRAYS
        self.t_arr = []
        self.m_dot_arr = []
        self.fluid_v_arr = []
        self.pump_v_arr = []
        self.hf_arr = []
        self.rho_arr = []
        
        # COMPONENT TEMPS
        self.cooler_T_arr = []
        self.heater_core_T_arr = []
        self.heater_ins_T_arr = []
        self.heater_outer_T_arr = []

        # PUMP POWER
        self.hyd_pwr_arr = []
        self.shaft_pwr_arr = []
        self.elec_pwr_arr = []
        
        
        
        
    def find_m_dot(self):
        
        """
        Calculates the mass flow rate in the loop using the Newton–Raphson method.
        The function computes the value of mass flow rate that satisfies the equation
        (pressure rise - frictional pressure losses) = 0. For each iteration, the difference 
        between pump pressure and frictional pressure loss are evaluated at f(ṁ), f(ṁ + h), 
        and f(ṁ - h). This is done to estimate the derivative required for the Newton step. 
        A limit is applied to restrict how quickly the mass flow rate can fluctuate 
        between iterations to improve numerical stability. The mass flow rate of the system
        is then set equal to the final converged value of ṁ.

        Returns
        -------
        m_dot_trial : float
            Updated mass flow rate (kg/s)
        """
            
        rho = np.mean(self.model.rho_T)
        m_dot_trial = self.model.m_dot
        h_step = 1e-3

        # NEWTON-RAPHSON
        for i in range(5):
            
            self.model.m_dot = m_dot_trial
            
            # f(ṁ)
            F = self.pump.press_diff(m_dot_trial, rho) - self.model.fric_press_diff()
            
            # f(ṁ + h)
            F_plus = self.pump.press_diff(m_dot_trial + h_step, rho) - self.model.fric_press_diff()
            
            # f(ṁ - h)
            F_minus =self.pump.press_diff(m_dot_trial - h_step, rho) - self.model.fric_press_diff()
            
            # f'(ṁ)
            F_prime = (F_plus - F_minus)/(2 * h_step)
            
            if F_prime == 0:
                break

            m_dot_new = m_dot_trial - (F/F_prime)
            
            # LIMIT ON HOW FAST MASS FLOW RATE CAN FLUCTUATE
            Δ_max = 0.1 * self.model.m_dot
            Δ_m_dot = max(-Δ_max, min(Δ_max, m_dot_new - self.model.m_dot))
            
            m_dot_new = self.model.m_dot + Δ_m_dot
        
            if abs(m_dot_new - m_dot_trial) < 1e-3:
                m_dot_trial = m_dot_new
                break
            
            m_dot_trial = m_dot_new
        
        self.model.m_dot = m_dot_trial  
            
        return m_dot_trial
    



    def append_attributes(self, rho, m_dot):
        
        """
        Appends the model's main attribute's respective history arrays. This 
        includes heater and cooler inner component temperatures, time step,
        mass flow rate, fluid velocity, pump PID output, head loss, density, 
        and the hydraulic, shaft and electrical power of the pump. These are later used
        for plotting and analysis.

        Parameters
        ----------
        rho : float
            Fluid density (kg/m^3)
        m_dot : float
            Mass flow rate of fluid (kg/s)
        """

        Q = (m_dot)/rho
        
        self.heater_core_T_arr.append(self.heater.core_T)
        self.heater_ins_T_arr.append(self.heater.ins_T)
        self.heater_outer_T_arr.append(self.heater.outer_T)
        self.cooler_T_arr.append(self.cooler.cooler_T)
        
        self.t_arr.append(self.t_count)
        self.m_dot_arr.append(m_dot)
        
        self.fluid_v_arr.append(Q/self.model.pipe_cross_A)
        self.pump_v_arr.append(self.pump.pump_v_control)
        
        hf = np.sum(self.model.head_loss())
        self.hf_arr.append(hf)
        self.rho_arr.append(rho)
        
        hydraulic, shaft, elec = self.pump.pump_pwr(m_dot, rho)
        
        self.hyd_pwr_arr.append(hydraulic)
        self.shaft_pwr_arr.append(shaft)
        self.elec_pwr_arr.append(elec)
        
        
        
        
    def run_sim(self, dur, dt):
        
        """
        This function runs the simulation loop for the specified duration (dur). At the start of the
        simulation, all PID controller values for the heater and cooler are reset. 
        For each timestep, the model updates: 
            
            - Temperature dependent fluid properties
            - Calculates mass flow rate 
            - Applies PID control signals for the pump, heater, and cooler
        
        After this, the main fluid physics loop is executed, pump power calculations are made
        and component temperature updates are computed. All relevant system attributes are 
        logged by calling 'append_attributes'. The simulation time counter incrementally 
        increases until the total duration is reached, at which point the simulation ends.

        Parameters
        ----------
        dur : float
            The duration of the simulation (s)
        dt : float
            Time step duration in seconds (dt = 0.1s)
        """
        
        
        # RESET PID VALUES
        self.heater.integral = 0.0
        self.heater.error = 0.0
        self.cooler.integral = 0.0
        self.cooler.error = 0.0

        while self.t_count < dur:
            
            self.model.update_variable_T()
            
            rho = np.mean(self.model.rho_T)
            m_dot = self.find_m_dot()
            self.model.m_dot = m_dot
            
            # RUNNING ALL PID CONTROLLERS
            self.pump.pump_pid_control(dt)
            self.heater.heater_pid_control(dt)
            self.cooler.cooler_pid_control(dt)
            
            # RUNNING MAIN PHYSICS LOOP ("pipe_phys" IS ALREADY CALLED IN "fluid_phys")
            self.model.fluid_phys(dt)
            self.pump.pwr_heating(m_dot, rho, dt)
            
            # COMPONENT TEMP UPDATES
            self.heater.update_heater_T(dt)
            self.cooler.update_cooler_T(dt)
            
            self.append_attributes(rho, m_dot)
            
            self.t_count += dt
        

            

if __name__ == "__main__":      # RUNNING THE FILE AUTOMATICALLY CALLS "sim()"
    
    from model import Model
    from pump import Pump
    from heater import Heater
    from cooler import Cooler


    from plots import Plot
    
    model = Model()
    pump = Pump(model)
    heater = Heater(model)
    cooler = Cooler(model)
    
    simulation = Simulation(model, heater, cooler, pump)
    simulation.run_sim(1200, 0.1)      # SIM LASTS 1200s WITH dt = 0.1s
    
    plots = Plot(model, simulation, heater, cooler, pump)
    plots.plot_all()