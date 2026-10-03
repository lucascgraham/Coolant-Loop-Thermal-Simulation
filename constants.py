import numpy as np

steps = 100
pump_steps = range(0, 2)
heater_steps = range(25, 27)
cooler_steps = range(75, 85)

stef_boltz = 5.67e-8    # STEFAN-BOLTZMANN CONSTANT
g = 9.81                # GRAVITATIONAL ACCELERATION
env_T = 293.15          # 20°C

h_air = 7.5             # CONVECTION COEFFICIENT OF AIR AT 20°C




                        ### PIPE/FLUID SPECIFIC CONSTANTS ###
                        
pipe_L = 20             # TOTAL PIPE LENGTH
pipe_d = 5e-2           # PIPE (INNER) DIAMETER
pipe_t = 6e-3           # PIPE THICKNESS

pipe_rough = 2e-6       # PIPE ROUGHNESS

rho_pipe = 7.85e3       # DENSITY OF STEEL PIPING AT 20°C
cp_pipe = 490           # SPECIFIC HEAT OF STEEL PIPING AT 20°C
kappa_pipe = 40         # THERMAL CONDUCTIVITY OF STEEL PIPING AT 20°C
pipe_epsilon = 0.1      # EMISSIVITY OF STEEL PIPING AT 20°C

rho_fluid = 1e3         # DENSITY OF WATER AT 20°C
fluid_m = rho_fluid * np.pi * ((pipe_d/2)**2) * pipe_L       # TOTAL MASS OF FLUID




                        ### HEATER SPECIFIC CONSTANTS ###

heater_max_pwr = 3e3    # 3kW IMMERSION HEATER
heater_eff = 0.95

outer_r = 0.3 * (pipe_d/2)
ins_r = outer_r * (2/3)
core_r = outer_r * (1/3)

outer_L = (len(heater_steps)/steps) * pipe_L
ins_L = outer_L * (2/3)
core_L = outer_L * (1/3)
   
outer_rho = 7.75e3      # DENSITY OF STAINLESS STEEL
outer_cp = 0.460e3      # SPECIFIC HEAT OF STAINLESS STEEL AT 20°C     
outer_epsilon = 0.16    # EMISSIVITY OF STAINLESS STEEL AT 20°C

ins_rho = 3.58e3        # DENSITY OF MAGNESIUM OXIDE
ins_cp = 0.983e3        # SPECIFIC HEAT OF MAGNESIUM OXIDE AT 20°C
ins_kappa = 50          # THERMAL CONDUCTIVITY OF MAGNESIUM OXIDE AT 20°C

core_rho = 8.4e3        # DENSITY OF NICHROME
core_cp = 0.44e3        # SPECIFIC HEAT OF NICHROME AT 20°C
core_kappa = 11         # THERMAL CONDUCTIVITY OF NICHROME AT 20°C
      



                        ### COOLER SPECIFIC CONSTANTS ###
             
k_cool = 5
cooler_max_pwr = 1e3  # 1kW POWERED COOLER
cooler_eff = 0.95
cooler_L = (len(cooler_steps)/steps) * pipe_L
cooler_rho = 8.96e3     # DENSITY OF COPPER AT 20°C
cooler_cp = 0.385e3     # SPECIFIC HEAT OF COPPER AT 20°C
cooler_kappa = 400      # THERMAL CONDUCTIVITY OF COPPER AT 20°C
cooler_epsilon = 0.07   # EMISSIVITY OF COPPER AT 20°C




                        ### PUMP SPECIFIC CONSTANTS ###
               
pump_eff = 0.9
pump_max_v = 3