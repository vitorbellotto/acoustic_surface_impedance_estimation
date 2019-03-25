#!/usr/bin/python3

import numpy as np
import matplotlib.pyplot as plt
import scipy.special as sp

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
plt.plot(x,sp.j0(x))
plt.plot(x,sp.j1(x))
plt.show()

material_A = Absorber( 5e3, 0.975, 1.05, 150e-6, 300e-6 )
#material_A.komatsu( omega )
#material_A.zwikker_kosten( omega )
#material_A.plot_propagation_constant(omega)
