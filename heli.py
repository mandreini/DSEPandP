# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 09:26:19 2018

@author: mandreini
"""
import numpy as np
import matplotlib.pyplot as plt


class heli(object):
    
    
    def __init__(self, Nr, Nb, W, h, R, c, f, P_e, eta_m, Omega, CDS, CDp_bar, k=1.15):
        '''
        Initialise a helicopter object
        
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
        self.Nr0, self.Nb0, self.W0, self.h0, self.R0, self.c0, self.f0, self.P_e0, self.eta_m0, self.Omega0, self.CDS0, self.CDp_bar0, self.k0 = (
            Nr, Nb, W, h, R, c, f, P_e, eta_m, Omega, CDS, CDp_bar, k)  # record the given values
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
        
    def _reset_values(self):
        '''Reset all the parameters of the HAMRAC to their initial ones'''
        self.Nr, self.Nb, self.W, self.h, self.R, self.c, self.f, self.P_e, self.eta_m, self.Omega, self.CDS, self.CDp_bar, self.k, self.alpha = (
        self.Nr0, self.Nb0, self.W0, self.h0, self.R0, self.c0, self.f0, self.P_e0, self.eta_m0, self.Omega0, self.CDS0, self.CDp_bar0, self.k0, self.alpha0)
        
    def _calc_params(self):
        '''Calculcate design parameters. Allows for them to change as new variables are introduced'''
        self.sigma = self.Nb * self.c / (np.pi*self.R)
        self.v_tip = self.Omega*self.R
        self.rho = 1.225*288.16/(15+271.16-1.983*self.h/304.8)*(1-(0.001981*self.h/0.3048)/288.16)**5.256
        self.vi_hov = (np.sqrt(self.W/(2*self.rho*np.pi*self.R**2)))
        self.P_a = self.P_e*self.eta_m

    def _get_vihov(self, R):
        return np.sqrt(self.W/(2*self.rho*np.pi*R**2))
    
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
        xmin = x[np.where(y==ymin)][0]
        return xmin, ymin
    
    def gen_pcurve(self, v, Pto, Pi, Ppd, Ppar, figtitle, fname):
        '''
        Plot a figure
        '''
        Pto, Pi, Ppd, Ppar = Pto/1000, Pi/1000, Ppd/1000, Ppar/1000
        plt.figure()
        xscale=1.3
        yscale=1.1
        plt.xlim(v[0], v[-1]*xscale)
        plt.ylim(0, max(Pto)*yscale)
        plt.grid()
        plt.xlabel('airspeed V [m/s]')
        plt.ylabel('power P [kW]')
        
        plt.title(figtitle)
        plt.plot(v, Pto, label='Pto')
        plt.plot(v, Pi, label='Pi')
        fname = 'power_req_HAMRAC.png'
            
        plt.plot(v, Ppd, label='Ppd')
        plt.plot(v, Ppar, label='Ppar')
        
        v_A, v_B, v_C, Pmin, Psfr = self.determineV_chars(v, Pto)
        
        plt.axhline(self.P_a, 0, v_C/max(v)/xscale, label='Pa', color='blue', linestyle='dotted')
        plt.axvline(v_C, 0, self.P_a/max(Pto)/yscale, label='vmax', color='blue', linestyle='dashed')
        plt.plot((0,v_B), (0, Psfr), color='purple')
        plt.axvline(v_B, 0, Psfr/max(Pto)/yscale, label='vsfr', color='purple', linestyle='dashed')
        plt.axvline(v_A, 0, Pmin/max(Pto)/yscale, label='vpmin', color='brown', linestyle='dashed')
        plt.axhline(Pmin, 0, v_A/max(v)/xscale, color='brown', linestyle='dotted')
        plt.legend(loc='lower right')
        plt.savefig(fname)
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
        p_min = min(P_to)
        v_A = v[np.where(P_to==p_min)]
        
        # minimum fuel flow airspeed (maximum range airspeed)
        # method: determine slope at all points on graph, determine slope of line between origin and said point, then find which one is the closest
        dydx = [(P_to[i+1]-P_to[i-1])/(v[i+1]-v[i-1]) for i in range(1, len(P_to)-1)]
        m = [P_to[j]/v[j] for j in range(1, len(P_to)-1)]
        difference = [abs((dydx[k]-m[k])/m[k]) for k in range(len(m))]
        closest_fit = min(difference)
        closest_loc = difference.index(closest_fit)
        p_sfr = P_to[closest_loc]
        v_B = v[closest_loc]
        
        # maximum airspeed
        v_C = v[np.where(P_to>P_a)[0][0]]  # does not work when hover is critical
        
        # climb rate
        v_cl = (P_a-P_to)*2/self.W
        v_cl = np.where(v_cl<0, 0, v_cl)  # cunha 16
        
        return v_A, v_B, v_C, p_min, p_sfr
    
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
        
    def determineP_to(self, v_air, P_to0=None):
        ''' 
        Determine the total power required of a helicopter in level, forward flight
        
        Parameters
        ----------
        v_inf : array_like
            The airspeed(s) of the helicopter in m/s 
            
        Returns
        -------
        P_to : array_like
            The total power required for the helicopter at airspeed
        P_i : array_like
            The induced power for the helicopter at airspeed
        P_pd : array_like
            The profile drag power for the helicopter at airspeed
        P_par : array_like
            The parasite drag power for the helicopter at airspeed
        '''
        self._calc_params()
        
        Pexcess = 0 if P_to0 is None else self.P_a - P_to0
        v_cl = 2*Pexcess/self.W
        v_inf = np.sqrt(v_air**2 + v_cl**2)
        
        myu = v_inf/(self.Omega*self.R)
        P_par = self.CDS*0.5*self.rho*v_inf**3
        alpha_tpp = (Pexcess+P_par)/self.W/v_inf
        
        vi_cr = self.vi_hov**2 / np.sqrt((v_inf*np.cos(alpha_tpp))**2+(v_inf*np.sin(alpha_tpp)+self.vi_hov)**2)  # from cunha lecture 6
        P_i = self.W*vi_cr
        P_pd = self.sigma*self.CDp_bar/8 * self.rho * ((self.Omega*self.R)**3) * np.pi*self.R**2 * (1+3*myu**2)  # to be given from Aero department, note: highly sensity to blade radius (R^5)!
        
        P_to = P_i + P_pd + P_par # marilena
        
        return P_to, P_i, P_pd, P_par
    
    def powerCurve(self, vs=np.arange(0, 105, 1), P_to0=0, figtitle='Power requirements for a helicopter in level forward flight', fname='power_req_HAMRAC.png'):
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
        
        power_vals = self.determineP_to(vs, P_to0)
        self.gen_pcurve(vs, *power_vals, figtitle, fname)
        return power_vals[0]
            
    def idealPhovafoR(self, Rs=np.linspace(1,11)):
        '''
        This function will generate a graph P_hov(R)
        
        Parameters
        ----------
        Rs: array_like
            An array of radii of the rotorblades
            
        Returns
        -------
        Rs : array_like
            The radii used in the procedure
        Pi_hovs : array_like
            The induced power values calculated
        '''
        self.setR(Rs)
        self._calc_params()
        Pi_hovs = self.determineP_to(v_air=0)
        
        return Rs, Pi_hovs
            
HAMRAC = heli(1, 8, 2657*9.81, 3500, 8.49, 0.457, 3, 632, 0.95, 21.21, 2, 0.012, 1.15)
airspeed = np.arange(1, 80, 1)
level_power = HAMRAC.powerCurve(airspeed, figtitle='Power requirements for a rotorcraft in level, horizontal, forward flight', fname='PowerCurveHAMRACLevel.png')
climbing_power = HAMRAC.powerCurve(airspeed, level_power, 'Power requirements for a rotorcraft in non-level, climbing flight', fname='PowerCurveHAMRACClimb.png')
#Pcurve2 = HAMRAC.powerCurve(airspeed, P_to0=Pcurve[1])
#HAMRAC.gen_pcurve(airspeed, *Pcurve2[1:], figtitle='Power requirements for a rotorcraft in non-level, climbing, forward flight', fname='PowerCurveHAMRACClimb.png')  # works!!!!!!
#Rdata = HAMRAC.idealPhovafoR()
#HAMRAC.plot_val(Rdata[0], Rdata[1][0]/1000, title='Optimizing power required to hover wrt blade radius', xlabel='Blade radius R [m]', ylabel='P_hov [kW]', get_min=True)
