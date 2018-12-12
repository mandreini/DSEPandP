# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 09:26:19 2018

@author: mandreini
"""
import numpy as np
import matplotlib.pyplot as plt


class heli(object):
    
    
    def __init__(self, name, Nr, Nb, W, h, R, c, f, P_e, eta_m, Omega, CDS, CDp_bar, k=1.15):
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
        
    def _reset_values(self):
        '''Reset all the parameters of the HAMRAC to their initial ones'''
        self.name, self.Nr, self.Nb, self.W, self.h, self.R, self.c, self.f, self.P_e, self.eta_m, self.Omega, self.CDS, self.CDp_bar, self.k, self.alpha = (
        self.name0, self.Nr0, self.Nb0, self.W0, self.h0, self.R0, self.c0, self.f0, self.P_e0, self.eta_m0, self.Omega0, self.CDS0, self.CDp_bar0, self.k0, self.alpha0)
        
    def _calc_params(self):
        '''Calculcate design parameters. Allows for them to change as new variables are introduced'''
        self.sigma = self.Nb * self.c / (np.pi*self.R)
        self.v_tip = self.Omega*self.R
        self.rho = 1.225*288.16/(15+271.16-1.983*self.h/304.8)*(1-(0.001981*self.h/0.3048)/288.16)**5.256
        self.vi_hov = (np.sqrt(self.W/(2*self.rho*np.pi*self.R**2)))
        self.P_a = self.P_e*self.eta_m

    def _get_vihov(self, R):
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
        xscale=1.3
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
        
        plt.axhline(self.P_a, 0, v_C/max(v)/xscale, label='Pa', color='blue', linestyle='dotted')
        plt.axvline(v_C, 0, self.P_a/max(Pto)/yscale, label='vmax', color='blue', linestyle='dashed')
        plt.plot((0,v_B), (0, Psfr), color='purple')
        plt.axvline(v_B, 0, Psfr/max(Pto)/yscale, label='vsfr', color='purple', linestyle='dashed')
        plt.axvline(v_A, 0, Pmin/max(Pto)/yscale, label='vpmin: %i m/s' % int(v_A), color='brown', linestyle='dashed')
        plt.axhline(Pmin, 0, v_A/max(v)/xscale, color='brown', linestyle='dotted')
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
        equal_power = np.where(P_to>P_a)
        if len(equal_power) > 0:
            v_C = v[np.where(P_to>P_a)[0][0]]  # does not work when hover is critical
        else:
            v_C = max(v)
        
        return v_A, v_B, v_C, p_min, p_sfr
        
    def determineP_to(self, v_air, P_to0=None):
        ''' 
        Determine the total power required of a helicopter in level, forward flight
        
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
        self._calc_params()
        
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
        self._calc_params()
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
        self._calc_params()
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
        self._calc_params()
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
        self._calc_params()
        vs, Ps = self.determineP_to(v_air=vs)
        vchars = self.determineV_chars(v=vs[0], P_to=Ps[0])
        # fix the return for the plot vals
        return vs, (Rs, vchars)
    
            
HAMRAC = heli('HAMRAC', 1, 8, 2657*9.81, 8950, 8.49, 0.457, 3, 450, 0.95, 21.21, 2, 0.012, 1.15)
marilena = heli("Marilena's Example", 1, 4, 90e3, 0, 9.5, 0.457, 3, 2000, 0.95, 21, 3, 0.01, 1.15)
airspeed = np.arange(1, 105, 1)
airspeed2 = np.arange(1, 120, 1)
#mp = marilena.powerCurve(airspeed2)

#level_speed, level_power = HAMRAC.powerCurve(airspeed, P_to0=None, figtitle='Power requirements for a rotorcraft in level, horizontal, forward flight', fname='PowerCurveHAMRACLevel.png')
#level_vchars = HAMRAC.determineV_chars(airspeed, level_power/1000)
#climb_speed, climbing_power = HAMRAC.powerCurve(airspeed, level_power, 'Power requirements for a rotorcraft in non-level, climbing flight', fname='PowerCurveHAMRACClimb.png')
#climb_vchars = HAMRAC.determineV_chars(airspeed, level_power/1000)
#
#Rvs, Rdata = HAMRAC.idealPhovafoR()
#HAMRAC.plot_val(Rdata[0], Rdata[1][0]/1000, title='Optimising power required to hover wrt blade radius', xlabel='Blade radius R [m]', ylabel='P_hov [kW]', get_min=True)
#
#cvs, cdata = HAMRAC.idealPhovafoc()
#HAMRAC.plot_val(cdata[0], cdata[1][0]/1000, title='Optimising power required to hover wrt blade chord', xlabel='Blade chord c [m]', ylabel='P_hov [kW]', get_min=True)
#
#Omegavs, Omegadata = HAMRAC.idealPhovafoOmega()
#HAMRAC.plot_val(Omegadata[0], Omegadata[1][0]/1000, title='Optimising power required to hover wrt angular velocity', xlabel='Angular velocity Omega [rad/s]', ylabel='P_hov [kW]', get_min=True)

v_Avs, v_Adata = HAMRAC.v_charsafoR(airspeed)#np.linspace(1,100,10), Rs=np.arange(5,10,1))
v_Avals = v_Adata[1][::-1]
HAMRAC.plot_val(v_Adata[0][0], v_Avals, xlabel='Radius of rotor R (m)', ylabel='v_pmin; minimum power airspeed', title='Comparison between minimum power airspeed and rotor diameter')

#h2 = HAMRAC
#h2.setR(1)
#h2p1 = h2.powerCurve()
#h2.setR(6)
#h2p2 = h2.powerCurve()
#h2.setR(7)
#h2p3 = h2.powerCurve()
#h2.setR(8)
#h2p4 = h2.powerCurve()
#h2.setR(9)
#h2p5 = h2.powerCurve()