#!/usr/bin/python3

print("This script calculates the surface impedance and propagation constant of an absorber based on the Komatsu Model.")

import numpy as np
import matplotlib.pyplot as plt

# Enviromental conditions
speed_of_sound = 343.2
air_density = 1.2041 

# Material A properties
porosity = 0.975
toruosity = 1.05
airflow_resistivity = 5000
viscious_char_length = 150e-6
thermal_char_length = 300e-6

# Domain
omega = np.linspace(1,2e4,80000)

prop_const = (omega / speed_of_sound ) * ( 0.0069  * np.power( 2 - np.log10( omega / ( 2 * np.pi * airflow_resistivity ) ), 4.1 ) + 1j * ( 1 + 0.0004 * np.power( 2 - np.log10( omega / ( 2 * np.pi * airflow_resistivity ) ), 6.2 ) ) )
surface_impedance = ( air_density * speed_of_sound ) * ( 1 + 0.00027 * np.power( 2 - np.log10( omega / ( 2 * np.pi * airflow_resistivity ) ), 6.2 ) - 1j * 0.0047 * np.power( 2 - np.log10( omega / ( 2 * np.pi * airflow_resistivity ) ), 4.1 ) ) 
#plt.plot(omega,prop_const)
#plt.title('Propagation constant')
#plt.xlabel('omega')
#plt.show()

plt.semilogx(omega,surface_impedance)
plt.grid()
plt.show()
