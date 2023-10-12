import numpy as np
import matplotlib.pylab as plt 

# Testing values
A_0 = 100
phi_initial = 0.0

## constants
Domain_length = 2 # meters
NBL = 200  # number of grid blocks
assumed_unreactive_mineral_fraction = 0.4
porosity_initial = 0.6
poro_minimum = 0.05
bulk_volume = Domain_length/NBL

k_g = 3e-8 # mol m^2 /s
S_m = 30 # kg / m^3
rho_m = 2071.87
V_m = 29e-6  # molar_volume [m3/mol]
k_n = 1/bulk_volume
TC = 25  # C
TK = TC + 273.15 # K 
gas_constant = 8.3144
pi_value = np.pi
NA = 6.022e23  # avogadros_constant
surface_energy = 100e-3  # milli Joule/m^2
hetero_phi = 0.28
shape_factor = (16 * pi_value)/3

# Solute concentrations
Ceq = 1.934e-3   # mol/m3 
C_inj = 0.25 

# space discretization 
nodes = NBL + 1  # (N)
x = np.linspace(0, Domain_length, nodes)
dx = x[1] - x[0]

# Time discretization  
dt = 10.0 # seconds
final_time = 172800 # seconds
num_steps = int(final_time / dt) + 1
time = np.linspace(0, final_time, num_steps)

##  calculate  nucleation rate without (In(C/Ceq))**2. This will be included later
top_value = shape_factor * V_m**2 * NA * surface_energy**3 * hetero_phi
bottom_value = (gas_constant * TK)**3
beta = top_value/bottom_value

fig, ax = plt.subplots(figsize = (10,8), nrows = 2)

# state variables
poro_curr = porosity_initial * np.ones_like(x)
poro_next = np.zeros_like(x) 
phi_curr = np.zeros_like(x)  # volume fraction of the secondary mineral, which is zero at the initial condition.
phi_next = np.zeros_like(x)
velocity = 1e-4 * np.ones_like(x)   # 2.54e-4 ; 5e-4 ; 0.01; 0.03

C_intial = Ceq
C_curr = C_intial * np.ones_like(x)
C_next = C_intial * np.ones_like(x)
C = C_intial
#Area = np.zeros_like(x) 
Area = A_0 * np.zeros_like(x) 
mole_value = 0.0 * np.zeros_like(x)   # mMoles of the minerals reacted. Thus, moles of minerals produced or dissolved
moles = 0.0


# Runge Kutta Variables
k1 = 0.0
k2 = 0.0
k3 = 0.0
k4 = 0.0
growth_rate = 0.0
nucleation_rate = 0.0
Concentration = 0.0
rate = 0.0


# Calculate reaction rate using Runge Kutta method
def calculate_growth_rate(C):
   rate = (k_g * Area)/Ceq * (Ceq - C)
   return rate

def calculate_nucleation_rate(C):
   rate = k_n * np.exp(-beta/np.log(C_curr/Ceq)**2)
   if np.isnan(rate):
      rate = 0.0 
   return rate   
   
          

def rk_method(C):  
      k1 = dt * calculate_growth_rate(C)
      k2 = dt * calculate_growth_rate(C + 0.5 * k1)
      k3 = dt * calculate_growth_rate(C+ 0.5 * k2)
      k4 = dt * calculate_growth_rate(C + k3)
     
      C += (k1 + 2 * k2 + 2 * k3 + k4) / 6.0 
      return C      
          

# Start flow and transport calculations

for it, time_val in enumerate(time):   
        # Surafce area calculations
    for ix, val in enumerate(x):
        #Area[ix] = S_m * rho_m * phi_curr[ix]
        Area[ix] = A_0  #* (phi_curr[ix]/phi_initial)**0.6667

        # Reactive Transport calculations 
    for ix, val in enumerate(x):
        if ix == 0:
            C_curr[ix] = C_inj   # inlet BC
        elif ix == nodes - 1:
            C_curr[ix] = C_curr[ix - 1] # outlet BC 
        else:    
            RHS_1 = - velocity[ix] * (C_curr[ix] - C_curr[ix - 1])/dx
            C_next[ix] =  (poro_curr[ix]*C_curr[ix] + dt * RHS_1)/poro_curr[ix]
            C = C_next[ix]
            
        # Now call calculate rate using Runge Kutta_method
            Concentration = rk_method(C)
            growth_rate = calculate_growth_rate(Concentration)
            #nucleation_rate = calculate_nucleation_rate(Concentration)
            moles += (growth_rate) * dt * (-1)  # Since the growth rate is negative
            mole_value[ix] = moles         
            RHS = RHS_1 + growth_rate #+ nucleation_rate
            
            phi_next[ix] = phi_curr[ix] + (RHS * dt * V_m)   # volume fraction             
            C_next[ix] =  (poro_curr[ix]*C_curr[ix] + dt * RHS)/poro_curr[ix]
            
              
            # Volume fraction calculations

    poro_curr = porosity_initial - phi_next 

    if it == 0:
        ax[0].plot(x,C_curr)
        ax[1].plot(x,phi_curr)
        
    if time_val % 8640 == 0:
        print(time_val, RHS_1, C_curr[10])
        np.savetxt('transport.dat', np.c_[x, poro_curr, phi_next, Area], delimiter = "  ")
        ax[0].plot(x,C_curr)
        ax[1].plot(x,phi_curr)

    phi_curr = phi_next
    C_curr = C_next

plt.savefig("solution.png") 
