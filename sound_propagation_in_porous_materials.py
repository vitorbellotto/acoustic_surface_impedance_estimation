#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp
import scipy.integrate as integ

class Absorber:

    def __init__( self, airflow_resistivity, porosity, tortuosity, visc_char_length, thrm_char_leng ):
        self.airflow_resistivity = airflow_resistivity
        self.porosity = porosity
        self.tortuosity = tortuosity
        self.visc_char_leng = visc_char_length
        self.thrm_char_leng = thrm_char_leng
        self.air_density = 1.2041
        self.speed_of_sound = 343.2
        self.gamma = 1.4
        self.air_viscosity = 18.4e-6
        self.atom_pres = 101325
        self.prandtl = 0.7179

    def cal_fluid(self, omega):
        self.char_dens = self.surface_impedance * self.prop_const * 1j / omega
        self.bulk_modu = - omega * 1j * self.surface_impedance / self.prop_const

    def komatsu(self, omega ):
        self.surface_impedance = ( self.air_density * self.speed_of_sound ) * ( 1 + 0.00027 * np.power( 2 - np.log10( omega / ( 2 * np.pi * self.airflow_resistivity ) ), 6.2 ) - 1j * 0.0047 * np.power( 2 - np.log10( omega / ( 2 * np.pi * self.airflow_resistivity ) ), 4.1 ) ) 
        self.prop_const = ( omega / self.speed_of_sound ) * ( 0.0069  * np.power( 2 - np.log10( omega / ( 2 * np.pi * self.airflow_resistivity ) ), 4.1 ) + 1j * ( 1 + 0.0004 * np.power( 2 - np.log10( omega / ( 2 * np.pi * self.airflow_resistivity ) ), 6.2 ) ) )
        self.cal_fluid( omega ) 

    def zwikker_kosten(self, omega ):
        Zo = 1
        ko = 1
        self.surface_impedance = ( Zo / self.porosity ) * np.sqrt( ( self.tortuosity / self.gamma ) * ( 1 - 1j * ( self.airflow_resistivity * self.porosity ) / ( omega * self.air_density * self.tortuosity ) ) )
        self.prop_const = 1j * ko * np.sqrt( self.tortuosity * self.gamma * ( 1 - 1j * ( self.airflow_resistivity * self.porosity ) / ( omega * self.air_density * self.tortuosity ) ) )

    def johnson_champoux_allard( self, omega ):
        self.char_dens = ( self.tortuosity * self.air_density / self.tortuosity ) * ( 1 + self.airflow_resistivity * self.tortuosity * np.sqrt( 1 + 1j * ( 4 * self.tortuosity ** 2 * self.air_viscosity * self.air_density * omega ) / ( self.airflow_resistivity ** 2 * self.visc_char_leng **2 * self.porosity )  ) / ( 1j * omega * self.air_density * self.tortuosity )) 
        self.bulk_modu =  ( self.gamma * self.atm_pressure / self.porosity ) / ( self.gamma - ( self.gamma - 1 ) * np.power( 1 - 1j * ( 8 * self.air_viscosity / ( self.thrm_char_leng ** 2 * self.prandtl * self.air_density * omega ) ) * np.sqrt( 1 + 1j * (( self.thrm_char_leng ** 2 * self.prandtl * self.air_density * omega ) / ( 16 * self.air_viscosity ) ) ), - 1 ) )

    def plot_surface_impedance(self, omega):
        plt.semilogx(omega,self.surface_impedance)
        plt.grid()
        plt.show()

    def plot_propagation_constant(self, omega):
        plt.semilogx(omega,self.prop_const)
        plt.grid()
        plt.show()


# Domain
omega = np.linspace(1,2e4,80000)
x = np.linspace(0.01,40,1000)

def third_term(k, z_source, z_receiver, air_wavenumber, absorber_wavenumber, r, absorber_density, air_density, layer_thickness ):
    air_nu = np.sqrt( np.power( k , 2 ) - np.power( air_wavenumber, 2 ) )
    absorber_nu = np.sqrt( np.power( k , 2 ) - np.power( absorber_wavenumber, 2 ) )
    return np.exp(- air_nu  * ( z_source + z_receiver ) ) * ( (2*absorber_density) / ( absorber_density * air_nu + air_density * absorber_nu * np.tanh(absorber_nu * layer_thickness ) ) ) * sp.j0( k * r ) * k 

def velocity_potencial(r_source, r_receiver, angular_frequency ):
    air_wavenumber = 2 * np.pi * angular_frequency / speed_of_sound
    absorber_wavenumber = 2 * np.pi * angular_frequency / ( 2 * speed_of_sound )
    r_refl = np.copy(r_source)
    r_refl[2] =  - r_refl[2]
    r = np.sqrt( np.power( r_source[0] - r_receiver[0], 2) + np.power( r_source[1] - r_receiver[1], 2) )
    Rinc = np.linalg.norm( r_source - r_receiver )
    Rrefl = np.linalg.norm( r_refl - r_receiver )
    Q = 1000
    inputs = ( r_source[2], r_receiver[2], air_wavenumber, absorber_wavenumber , r, absorber_density,) 
    #return ( Q/ ( 4 * np.pi ) ) * ( np.exp( - 1j * air_wavenumber * Rinc ) / Rinc )
    return ( Q/ ( 4 * np.pi ) ) * ( np.exp( - 1j * air_wavenumber * Rinc ) / Rinc - np.exp( - 1j * air_wavenumber * Rrefl ) / Rrefl  ) 

speed_of_sound = 343.2
angular_frequency = np.linspace(1,20000,80000)
r_source = np.array([0,0,0.6])
r_receiver = np.array([0,0,0.05])

plt.plot( angular_frequency , velocity_potencial( r_source, r_receiver, angular_frequency ) )
plt.show()












# Allard paper - The acoustic sound field above a porous layer and the estimation of the acoustic surface impedance from free-field measurements
#integral_bessel = integ.quadrature(velocity_potencial,1,2e4)
# pressure
# p_tot = air_density * omega * velocity_potencial * 1j



material_A = Absorber( 5e3, 0.975, 1.05, 150e-6, 300e-6 )
#material_A.komatsu( omega )
#material_A.zwikker_kosten( omega )
#material_A.plot_propagation_constant(omega)
