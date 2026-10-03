Coolant Loop Thermal Simulation





A segmented thermal-mechanical segmented coolant loop simulation modelling mass flow rate, system temperature distribution, pump behaviour, heat transfer and temperature based fluid variables. Pump curves, advection, frictional heating, conduction, convection, radiation, pump affinity laws and PID controlled heater/cooler components are also incorporated into the modelling to enhance the physical realism of the simulation





Features:



* Segmented loop with advection based temperature flow between neighbouring pipe sections
* Temperature dynamics modelling through convection, conduction and radiation heat transfer equations in loop components
* PID controlled pump, heater and cooler component behaviour
* Temperature based behaviour modelling of fluid variables
* Pump power and efficiency curves are implemented and plotted for analysis of pump performance
* Mass flow rate is estimated using a Newton-Raphson iteration method
* Automatic plotting of coolant loop characteristics (Mass flow rate, component temperatures, fluid velocity/temperature, pump curves, head loss etc)





Structure:



simulation.py - Where coolant model is simulated from and mass flow rate is estimated



model.py - Main physical loop is modelled (fluid/pipe physics)



pump.py - Pump power calculated and PID controller modelling



heater.py - Heat transfer to coolant fluid from heater modelled and heater temperature updated. Heater PID controller modelled also



cooler.py - Heat removal from coolant fluid via cooler modelled and cooler temperature updated. Cooler PID controller modelled also



constants.py - Where all physical constants and parameters are stored and imported from



plots.py - Imported into simulation.py to automatically plot modelled behaviour (mass flow rate, temperatures, pump curves, head loss, etc)





Installation:



This model requires the use of all packages listed in "requirements.txt". To install all required packages,

run the following command:



"pip install -r requirements.txt"



This will automatically install all necessary packages (NumPy, matplotlib and pandas).





Running Simulation:



To run the coolant loop simulation, execute the "simulation.py" file. The full model is ran by the 'if \_\_name\_\_ == "\_\_main\_\_":' block, which runs the simulation loop and automatically plots all relevant data.





Model Overview:



Fluid Physics - The main body of the simulation that governs the temperature of each fluid segment by modelling advection, heat gain from the heater, heat loss via the cooler, frictional heating due to flow resistance, and convection based heat loss to the surrounding piping. A copy of the array which holds fluid segment temperature is stored at the end of the function to be used in plotting and data analysis.



Pipe Physics - Exchanges heat with coolant fluid through convection and transfers heat to the innermost pipe layer. The heat of the inner piping is then transferred to the outer pipe layer via heat conduction and finally this heat is lost to the environment via radiation and convection. The inner and outer pipe temperatures are then updated each timestep, allowing the piping to influence the coolant temperature distribution along the loop.



Major Components - Designated heater, cooler and pump components are integrated into the overall coolant loop model. The heater and cooler components transfer and remove heat respectively to and from the coolant fluid. The internal temperature of these components, which influence the coolant fluid temperature are also modelled and updated. The pump component interpolates hydraulic, shaft and electrical power curves and transfers heat to the coolant fluid through power dissipation. All three components have integrated PID controllers which stabilise their outputs. The PID controller of the heater and cooler stabilise the temperature of the coolant fluid and the pump PID maintains a stable mass flow rate. Working together, these components create a regulation system that pushes the simulation toward thermal/hydraulic stability.



Newton-Raphson Method - The selected numerical method used to estimate the mass flow rate of the system. The mass flow rate is estimated ultimately by finding the point in a function of mass flow rate where the pump pressure rise and frictional pressure drop equalise. The solving loop used repeatedly evaluates this difference and its numerical derivative, updating the mass flow rate estimate at every iteration until the pressure difference falls below a certain tolerance level, at which point the loop breaks and a final value of mass flow rate is reached.



Temperature based variable behaviour - Several key coolant fluid properties can vary significantly with temperature, and these variations influence the loop behaviour. In this model, the temperature based nature of these properties has been has been implemented. Some properties are governed by well defined physical laws (density and specific heat), and these are applied directly to the model. For properties governed by more complex and empirical behaviour (convection coefficient), appropriate assumptions, correlations and data have all been used. All properties have been modelled so that they follow these behaviours and subsequently influence the surrounding physics used in the model realistically.





Assumptions Made:



* Water used as coolant is incompressible
* Kell equation assumes no air pockets are present in water
* The loop is treated essentially in one dimension
* Pipe geometry is assumed to be uniform along the whole path
* Each segment is assumed to have uniform temperature and fluid properties.
* Head loss is calculated via friction only and head loss due to pipe geometry and other components is not modelled
* Heater components are uniform along their length, instead of segmented sections which may have different temperatures
* It is assumed that the coolant fluid remains a liquid at all times and no mass is lost to boiling
* Properties of heater and coolant components (thermal conductivity, specific heat etc) are treated as constant (the physical changes are minor)





Future Additions:



* Heater/cooler geometry to be modelled so they act as realistic components.
* Add a component such as a reactor core which the coolant loop could influence and regulate
* Add head loss due to pipe geometry and other components such as valves and changing pipe materials

