# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 09:26:19 2018

@author: mandreini
"""
import abc
import numpy as np
import matplotlib.pyplot as plt


class SRheli(object):
    
    __metaclass__ = abc
    configuration = 'single-rotor'
    iscoax = False
    zR = 0
    
    def __init__(self, name, Nr, Nb, W, h, R, c, f, P_e, eta_m, Omega, CDS, CDp_bar, k=1.15, *args, **kwargs):
        '''
        Initialise a single rotor helicopter object
        
        Parameters
        ----------
        Nr : int
            Number of rotors
        Nb : int
            Number of blades
        W : int
            Weight of the rotorcraft [N]
        h : float
            Operating altitude of the helicopter
        R : float
            Radius of the blades
        c : float
            Chord length of the blade
        f : float
            Equivalent wetted area of the helicopter. Recommend 0.93 for small helicopter, and 4.65 for large helicopter        
        P_e : float
            Engine power
        eta_m : float
            Mechanical efficiency of the engine, 0<eta_m<1
        Omega : float
            Angular velocity of the rotor in rad/s          
        CDS : float
            Equivalent flat plate area of the helicopter. Values range from 1 (aerodynamicly efficient) to 6 (aerodynamically inefficient)        
        CDp_bar : array_like, float
            (Average) drag coefficient of the rotor blade. Constant when blade is constant        
        '''
        self.name0, self.Nr0, self.Nb0, self.W0, self.h0, self.R0, self.c0, self.f0, self.P_e0, self.eta_m0, self.Omega0, self.CDS0, self.CDp_bar0, self.k0 = (
            name, Nr, Nb, W, h, R, c, f, P_e, eta_m, Omega, CDS, CDp_bar, k)  # record the given values
        
        self.name = name
        self.Nr = Nr
        self.Nb = Nb
        self.W = W
        self.h = h
        self.R = R
        self.c = c
        self.f = f
        self.P_e = P_e
        self.eta_m = eta_m
        self.Omega = Omega
        self.CDS = CDS
        self.CDp_bar = CDp_bar
        self.k = k
        
    @staticmethod 
    def convert(Nr, Nb, R, c, Omega):
        '''
        Converts the given single-rotor data into the equivalent (solidity) coaxial rotor data
        
        Parameters
        ----------
        Nr : int
            Number of rotors
        Nb : int
            Number of blades (per rotor)
        R : float
            Blade radius
        c : float
            Blade chord length
        Omega : float
            Angular velocity of the rotor
            
        Returns
        -------
        Nr, Nb, R, c, Omega equivalents for a coaxial rotor
        
        '''
        return 2, Nb//2, R/np.sqrt(2), c/np.sqrt(2), Omega*np.sqrt(2)
        
    def reset_values(self):
        '''Reset all the parameters of the HAMRAC to their initial ones'''
        self.name, self.Nr, self.Nb, self.W, self.h, self.R, self.c, self.f, self.P_e, self.eta_m, self.Omega, self.CDS, self.CDp_bar, self.k, self.alpha = (
        self.name0, self.Nr0, self.Nb0, self.W0, self.h0, self.R0, self.c0, self.f0, self.P_e0, self.eta_m0, self.Omega0, self.CDS0, self.CDp_bar0, self.k0, self.alpha0)
        
    def calc_params(self):
        '''Calculcate design parameters. Allows for them to change as new variables are introduced'''
        self.sigma = self.Nb * self.c / (np.pi*self.R)
        self.v_tip = self.Omega*self.R
        self.rho = 1.225*288.16/(15+271.16-1.983*self.h/304.8)*(1-(0.001981*self.h/0.3048)/288.16)**5.256
        self.vi_hov = (np.sqrt(self.W/(2*self.rho*np.pi*self.R**2)))
        self.P_a = self.P_e*self.eta_m

    def get_vihov(self, R):
        return np.sqrt(self.W/(2*self.rho*np.pi*R**2))
        
    def setR(self, Rs):
        self.R = Rs

    def setOmega(self, Omegas):
        self.Omega = Omegas

    def setalpha(self, alphas):
        self.alpha = alphas    
    
    def setW(self, Ws):
        self.W = Ws
        
    def seth(self, hs):
        self.h = hs
        
    def setc(self, cs):
        self.c = cs
        
    def setCDp_bar(self, CDp_bars):
        self.CDp_bar = CDp_bars
        
    def setCDS(self, CDSs):
        self.CDS = CDSs
    
    def plot_val(self, x, y, title=None, xlabel=None, ylabel=None, get_min=False, xlim=None, ylim=None):
        '''Plot x against y'''
        plt.figure()
        plt.grid()
        
        if not xlim: xlim = (0, max(x))
        if not ylim: ylim = (0, max(y))
        plt.xlim(xlim)    
        plt.ylim(ylim)
        
        if get_min:
            xmin, ymin = self.get_min(x, y)
            plt.axhline(ymin, 0, xmin/xlim[1])#, linestlye='dotted')
            plt.axvline(xmin, 0, ymin/ylim[1])#, linestyle='dotted')
        
        if title: plt.title(title)
        if xlabel: plt.xlabel(xlabel)
        if ylabel: plt.ylabel(ylabel)
        plt.plot(x, y)
        plt.show()
    
    @staticmethod
    def get_min(x, y):
        ymin = min(y)
        xmin = x[np.where(y==ymin)]
        return xmin, ymin
    
    def gen_pcurve(self, v, Pto, Pi, Ppd, Ppar, figtitle, fname):
        '''
        Generate the power curve
        
        Parameters
        ----------
        v : array_like
            airspeed of the rotorcraft
        Pto : array_like
            Total power (required) of the rotorcraft
        Pi : array_like
            Induced power of the rotorcraft
        Ppd : array_like
            Profile drag power of the rotorcraft
        Ppar : array_like
            Parasite drag power of the rotorcraft
        figtitle : string
            What to name the graph
        fname : string
            What to title the file to save, None to not save the figure
            
        Returns
        -------
        None
        '''
        Pto, Pi, Ppd, Ppar = Pto/1000, Pi/1000, Ppd/1000, Ppar/1000
        plt.figure()
        xscale=1.5
        yscale=1.1
        plt.xlim(v[0], v[-1]*xscale)
        plt.ylim(0, max(Pto)*yscale)
        plt.grid()
        plt.xlabel('Airspeed V [m/s]')
        plt.ylabel('Power P [kW]')
        
        plt.title(figtitle)
        plt.plot(v, Pto, label='Pto')
        plt.plot(v, Pi, label='Pi')
            
        plt.plot(v, Ppd, label='Ppd')
        plt.plot(v, Ppar, label='Ppar')
        
        v_A, v_B, v_C, Pmin, Psfr = self.determineV_chars(v, Pto)
        
        plt.axhline(self.P_a, 0, v_C/v.max()/xscale, label='Pa: %i kW' % int(self.P_a), color='blue', linestyle='dotted')
        plt.axvline(v_C, 0, self.P_a/Pto.max()/yscale, label='vmax: %i m/s' % int(v_C), color='blue', linestyle='dashed')
        plt.plot((0,v_B), (0, Psfr), color='purple')
        plt.axvline(v_B, 0, Psfr/Pto.max()/yscale, label='vsfr: %i m/s' % int(v_B), color='purple', linestyle='dashed')
        plt.axvline(v_A, 0, Pmin/Pto.max()/yscale, label='vpmin: %i m/s' % int(v_A), color='brown', linestyle='dashed')
        plt.axhline(Pmin, 0, v_A/v.max()/xscale, label='Pmin: %i kW' % int(Pmin), color='brown', linestyle='dotted')
        plt.legend(loc='lower right')
        if fname: plt.savefig(fname)
        plt.show()

    def determineV_chars(self, v, P_to):
        '''
        Calculates relevant velocity values of a helicopter from its power requirements
        
        Parameters
        ----------
        P_to : array_like
            Power curve of the helicopter [W]
        v : array_like
            Airspeeds associated with the power curve [m/s]
        P_e : float
            Engine power of the helicopter [W]
        eta_m : float
            Mechanical efficiency
        W : float
            Weight of the helicopter [N]
            
        Returns
        -------
        v_A : float
            Minimum power/fuel flow/maximum climb/minimum descent airspeed [m/s]
        v_B : float
            Maximum specific flight range airspeed, airspeed for minimum descent angle [m/s]
        v_C : float 
            Maximum airspeed, 0 rate of climb airspeed [m/s]
        p_min : float
            Power required at v_A (minimum power) [W]
        P_sfr : float
            Power required at v_A (power required for maximum specific flight range) [W]
        '''
        # transmission losses
        P_a = self.eta_m * self.P_e
        
        # minimum power airspeed
        p_min = P_to.min(axis=0) 
        v_A = v[np.where(P_to==p_min)]
        
        # minimum fuel flow airspeed (maximum range airspeed)
        # method: determine slope at all points on graph, determine slope of line between origin and said point, then find which one is the closest
#        dydx = [(P_to[i+1]-P_to[i-1])/(v[i+1]-v[i-1]) for i in range(1, len(P_to)-1)]
#        m = [P_to[j]/v[j] for j in range(1, len(P_to)-1)]
#        difference0 = [abs((dydx[k]-m[k])/m[k]) for k in range(len(m))]
#        closest_fit = min(difference)
#        closest_loc = difference.index(closest_fit)
#        p_sfr = P_to[closest_loc]
#        v_B = v[closest_loc]
#        
        slope_power_to = np.gradient(P_to, axis=0)  # identical to dydx, but also includes first and last
        # plt.plot(v[:,0], P_to[:,0]/max(P_to[:,0]), v[:,0], sl1[:,0]/max(sl1[:,0]))
        slope_origin_to_point = P_to/v  # identical to m, but also includes first and last
#        P_to_I, P_to_II = np.meshgrid(P_to, P_to)
#        v_I, v_II = np.meshgrid(v, v)
#        dydx = (P_to_II-P_to_I)/(v_II-v_I)
#        m = P_to/v
        difference = np.abs((slope_power_to - slope_origin_to_point)/slope_origin_to_point)
        closest_fit = difference.min(axis=0)
        closest_loc = np.where(difference==closest_fit)   
        v_B = v[closest_loc]
        p_sfr = P_to[closest_loc]
        
        # maximum airspeed
        equal_power = np.where(P_to>P_a*1000)
        if len(equal_power) > 0:
            v_C = v[np.where(P_to>P_a)[0][0]]  # does not work when hover is critical
        else:
            v_C = max(v)
        
        return v_A, v_B, v_C, p_min, p_sfr
        
    def determineP_to(self, v_air, P_to0=None):
        ''' 
        Determine the total power required of a helicopter in level, forward flights
        
        Parameters
        ----------
        v_inf : array_like
            The airspeed(s) of the helicopter in m/s 
        P_to0 : array_like
            The power required for level flight. Used to determine climbing flight
            
        Returns
        -------
        v_air : array_like
            Airspeed values of the heli used
        v_cl : array_like
            Climb speed values of the heli calculated. 0 in horizontal flight (P_to0=None)
            
        P_to : array_like
            The total power required for the helicopter at airspeed
        P_i : array_like
            The induced power for the helicopter at airspeed
        P_pd : array_like
            The profile drag power for the helicopter at airspeed
        P_par : array_like
            The parasite drag power for the helicopter at airspeed
        '''
        self.calc_params()
        
        Pexcess = 0 if P_to0 is None else self.P_a*1000 - P_to0
        v_cl = 2*Pexcess/self.W
        v_inf = np.sqrt(v_air**2 + v_cl**2)
        
        myu = v_air/(self.Omega*self.R)
        P_par = self.CDS*0.5*self.rho*v_inf**3
        alpha_tpp = (Pexcess+P_par)/self.W/v_inf if P_to0 is not None else 0
        
        vi_cr = self.vi_hov**2 / np.sqrt((v_inf*np.cos(alpha_tpp))**2+(v_inf*np.sin(alpha_tpp)+self.vi_hov)**2)  # from cunha lecture 6
        P_i = self.W*vi_cr
        P_pd = self.sigma*self.CDp_bar/8 * self.rho * ((self.Omega*self.R)**3) * np.pi*self.R**2 * (1+3*myu**2)  # to be given from Aero department, note: highly sensity to blade radius (R^5)!
        
        P_to = P_i + P_pd + P_par # marilena
        
        return (v_air, v_cl), (P_to, P_i, P_pd, P_par)
    
    def powerCurve(self, vs=np.arange(0, 105, 1), P_to0=None, figtitle='Power requirements for a helicopter in level forward flight', fname=None):
        '''
        Generate a power curve for various airspeeds
        
        Parameters
        ----------
        vs : array_like
            The airspeed of the rotorcraft
            
        P_to0 : array_like
            The total_power required at level flight. Used to determine climbing flight
            
        figtitle : string
            The title of the figure on the plot
            
        fname : string
            The filename of the figure to save. Incude extension (i.e. 'name.png')
            
        Returns
        -------
        array_like
            The total power required for the helicopter at given airspeeds
        '''
        
        v_vals, power_vals = self.determineP_to(vs, P_to0)
        self.gen_pcurve(vs, *power_vals, figtitle, fname)
        return v_vals, power_vals[0]
            
    def idealPhovafoR(self, Rs=np.linspace(1,11)):
        '''
        This function will generate a graph P_hov(R)
        
        Parameters
        ----------
        Rs: array_like
            An array of radii of the rotorblades
            
        Returns
        -------
        vs : array_like
            The airspeed values used
        Rs : array_like
            The radii used in the procedure
        P_hovs : array_like
            The power at hover values calculated
        '''
        self.setR(Rs)
        self.calc_params()
        vs, P_hovs = self.determineP_to(v_air=0)
        
        return vs, (Rs, P_hovs)
    
    def idealPhovafoc(self, cs=np.linspace(0.1,0.6)):
        '''
        This function will generate a graph P_hov(c)
        
        Parameters
        ----------
        cs : array_like
            An array of blade chords
            
        Returns
        -------
        vs : array_like
            The airspeed values used
        cs :  array_like
            The blade chord lengths used
        Pi_hovs : array_like
            The power at hover values calculated
        '''
        self.setc(cs)
        self.calc_params()
        vs, P_hovs = self.determineP_to(v_air=0)
        
        return vs, (cs, P_hovs)
    
    def idealPhovafoOmega(self, Omegas=np.linspace(15, 25)):
        '''
        This function will generate a graph P_hov(Omega)
        
        Parameters
        ----------
        Omegas : array_like
            An array of angular velocities
            
        Returns
        -------
        vs : array_like
            The airspeed values used
        Omegas :  array_like
            The angular velocities used
        Pi_hovs : array_like
            The power at hover values calculated
        '''
        self.setOmega(Omegas)
        self.calc_params()
        vs, P_hovs = self.determineP_to(v_air=0)
        
        return vs, (Omegas, P_hovs)
    
    def v_charsafoR(self, vs, Rs=np.linspace(1,10,50)):
        '''
        This function will generate a graph v_A(R).
        
        Parameters
        ----------
        Rs : array_like
            An array of rotor blade radii
            
        Returns
        -------
        vs : array_like
            The airspeed values used
        Rs :  array_like
            The radii used
        v_As : array_like
            The maximum airspeed v_C values calculated
        '''
        Rs, vs = np.meshgrid(Rs, vs)
        self.setR(Rs)
        self.calc_params()
        vs, Ps = self.determineP_to(v_air=vs)
        vchars = self.determineV_chars(v=vs[0], P_to=Ps[0])
        return vs, (Rs, vchars)
    

class Coaxheli(SRheli):
    
    configuration = 'coaxial-rotor'
    iscoax = True
    zR = 0

    def __init__(self, *args, **kwargs):
        super(Coaxheli, self).__init__(*args, **kwargs)
        
        if 'zR' in kwargs.keys():
            self.set_zR(kwargs['zR'])
        
        if 'convert' in kwargs.keys() and kwargs['convert']:
            self.convert()

    def set_zR(self, zRs):
        self.zR = zRs
        
    def convert(self): 
        '''
        Converts the given single-rotor data into the equivalent (solidity) coaxial rotor data
        '''
        self.Nr = 2 
        self.Nb = self.Nb//2
        self.R = self.R/np.sqrt(2)
        self.c = self.c/np.sqrt(2)
        self.Omega = self.Omega*np.sqrt(2)
    
    def reconvert(self):
        '''
        Takes the (current) coaxial rotor parameters, and returns them into their single-rotor equivalent sizes
        '''        
        SR_Nr = 1
        SR_Nb = self.Nb*2
        SR_R = self.R*np.sqrt(2)
        SR_c = self.c*np.sqrt(2)
        SR_Omega = self.Omega/np.sqrt(2)
        
        return SR_Nr, SR_Nb, SR_R, SR_c, SR_Omega


