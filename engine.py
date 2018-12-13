# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 14:28:58 2018

@author: mandreini
"""

class Engine(object):
    
    def __init__(self, P_sl, name='not given', mass=None, sfc=None, size=None):
        '''
        Create an engine object
        To be changed as i.e. Max does the engine work, this is a basic structure
        
        Parameters
        ----------
        P_sl : float
            Power available at sea level
        weight : float
            Weight of the engine, in kg
        sfc : float
            Specific fuel consumption of the engine, in N/s
        '''
        
        self.P_sl = P_sl
        self.mass = mass
        self.sfc = sfc
        self.size = size
        
    @staticmethod
    def get_rho(h):
        return 1.225*288.16/(15+271.16-1.983*h/304.8)*(1-(0.001981*h/0.3048)/288.16)**5.256
    
    @staticmethod
    def get_T(h):
        return 298.13 + -0.0065*h
        
    def determine_power_at_altitude(self, h):
        '''
        Function to determine power provided by the engine at given altitude
        Power available decreases 1:1 with density decrease
        Power available increases 0.1% per degree K temperature decrease
        
        Parameters
        ----------
        altitude : array_like
            (An array of) operating altitude(s)
        
        Returns
        -------
        array_like :
            power provided at given altitude(s)
        '''
        rhoa = self.get_rho(h)
        rho0 = self.get_rho(0)
        Ta = self.get_T(h)
        T0 = self.get_T(0)
        
        return self.P_sl*(rhoa/rho0)*(1+(T0-Ta)*0.001)

lbshrhp2nsw = lambda x: 0.453592*0.277778*9.81/0.745699872*1000 * x   # lbs/hr/hp to N/s/W... I think 
   
HoneywellT35 = Engine(1175, name='AH-1J SeaCobra', mass=312, sfc=lbshrhp2nsw(0.65), size=(120/9, 24.5))  # https://aerospace.honeywell.com/en/products/engines/t53-turboshaft-engine, 
# https://aerospace.honeywell.com/en/~/media/aerospace/files/brochures/n61-1733-000-000-t53-v6-bro.pdf, https://en.wikipedia.org/wiki/Lycoming_T53
KlimovVK2500 = Engine(1985, name='Kamov Ka-50', mass=300, sfc=0.22*3.6, size=(205.5, 72.8))  # https://en.wikipedia.org/wiki/Klimov_VK-2500

print(HoneywellT35.determine_power_at_altitude(9000))

