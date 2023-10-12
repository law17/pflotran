# Mineral reaction code for a batch reaction model (i.e., transient state)
# TST mineral reaction rate equation: dc/dt = KA(1 - C/ceq) = KA/Ceq * (Ceq - C)
# Assume 1.0 m3 bulk volume 

import numpy as np
import matplotlib.pylab as plt

# Rate constants, concentrations, and area

growth_rate_constant = 3e-8      # mol/m2/s
Area = 10                        # m2/m3 of bulk
C = 0.15                         # Aqueous ionic species concentration, mol/m3
Ceq = 1.934e-3                   # equilibrium concnetration, mol/m3

# Time discretization 
final_time = 86400               # Total simulation time, seconds
time_step = 864                  # Time step, seconds
num_of_steps = int(final_time/time_step)
time = np.linspace(0, final_time, num_of_steps)


# Print options
output_time = 872.7272727272727            # Print results after every this number of time in secods
output_steps = int(final_time/output_time) # Total number of times to print results
time_hr = 0.0
seconds_to_hr = 0.000277778
prtingting = 1
period = 1

# Other variables 
remained_amount_of_C = 0.0       # Conectration of C that remained during the goechemical reaction for the time step
reacted_amount_of_C = 0.0        # Conectration of C that reacted during the goechemical reaction for the time step
reaction_rate = 0.0              # The reaction rate obtained during a specified time 
moles = 0.0                      # mMoles of the minerals reacted. Thus, moles of minerals produced or dissolved
k1 = 0.0
k2 = 0.0
k3 = 0.0
k4 = 0.0
    
# Initial conditions:  C = C at t = 0   

def calculate_rate(C):
   rate = (growth_rate_constant * Area)/Ceq * (Ceq - C)
   return rate
          

def rk_method(C):  
      k1 = time_step * calculate_rate(C)
      k2 = time_step * calculate_rate(C + 0.5 * k1)
      k3 = time_step * calculate_rate(C + 0.5 * k2)
      k4 = time_step * calculate_rate(C + k3)
     
      C += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0 
      return C      
          
          
# Now call "rk_method" and "calculate_rate" fucntions

for it, time_val in enumerate(time): 
  if time_val > 0.0:
    C = rk_method(C)
    reaction_rate = calculate_rate(C)
    moles += reaction_rate * time_step * (-1)  # Since precipitation rate is negative
    
    print(reaction_rate)
    
    # Now, output time, rate, and moles at every this number of steps   
    if period == prtingting:  
       time_hr = time_val * seconds_to_hr 
       np.savetxt('output.dat', np.c_[time_hr, reaction_rate, moles], delimiter = "  ")
       output_time = output_time * 11
       print(moles)
       prtingting += 1
    
    period += 1
        

     
         




