# -*- coding: utf-8 -*-
"""
Created on Thu Dec  6 09:26:19 2018

@author: mandreini
"""
import abc
import numpy as np
import matplotlib.pyplot as plt
import engine
import copy

class SRheli(object):
    
    __metaclass__ = abc.ABCMeta
    configuration = 'single-rotor'
    iscoax = False
    zR = 0
    gamma = 1.4
    v_cr = 140/1.9
    vmax = v_cr*1.15 
    vne = vmax*1.1

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
        
        self.P_a = P_e*eta_m
        self.engine_choice = None
        
    def __str__(self):
        return (self.name, self.Nr, self.Nb, self.W, self.h, self.R, self.c, 
                self.f, self.P_e, self.eta_m, self.Omega, self.CDS, self.CDp_bar, 
                self.k, self.P_a, self.engine_choice)
                
    def __repr__(self):
        return 'Name: %s\nNr: %i\nNb: %i\nh: %f\nR: %f\nc: %f\nOmega: %f\nP_e: %f\nP_a(h): %f\nEngine: %s' % (
        self.name, self.Nr/2, self.Nb*2, self.h, self.R, self.c, self.Omega, self.P_e, self.P_a, self.engine_choice.name)
                
#    @staticmethod 
#    def convert(Nr, Nb, R, c, Omega):
#        '''
#        Converts the given single-rotor data into the equivalent (solidity) coaxial rotor data
#        
#        Parameters
#        ----------
#        Nr : int
#            Number of rotors
#        Nb : int
#            Number of blades (per rotor)
#        R : float
#            Blade radius
#        c : float
#            Blade chord length
#        Omega : float
#            Angular velocity of the rotor
#            
#        Returns
#        -------
#        Nr, Nb, R, c, Omega equivalents for a coaxial rotor
#        
#        '''
#        return 2, Nb/2, R/np.sqrt(2), c/np.sqrt(2), Omega*np.sqrt(2)
        
    def reset_values(self):
        '''Reset all the parameters of the HAMRAC to their initial ones'''
        self.name, self.Nr, self.Nb, self.W, self.h, self.R, self.c, self.f, self.P_e, self.eta_m, self.Omega, self.CDS, self.CDp_bar, self.k, self.alpha = (
        self.name0, self.Nr0, self.Nb0, self.W0, self.h0, self.R0, self.c0, self.f0, self.P_e0, self.eta_m0, self.Omega0, self.CDS0, self.CDp_bar0, self.k0, self.alpha0)
        
    @staticmethod
    def get_rho(h):
        return 1.225*288.16/(15+271.16-1.983*h/304.8)*(1-(0.001981*h/0.3048)/288.16)**5.256

    @staticmethod
    def get_T(h):
        return (15+271.16-1.983*h/304.8)
        
    def get_a(self):
        return np.sqrt(self.gamma * 287 * self.get_T(self.h))
        
    def get_bl(self, v_air=v_cr):
        v_tip = self.v_tip + v_air
        return self.W/np.pi/self.R**2 / self.rho / v_tip / v_tip / self.sigma
        
    @staticmethod
    def get_vtip(R, Omega, v):
        return (Omega*R+v) / v

    def calc_params(self, v_air=v_cr):
        '''Calculcate design parameters. Allows for them to change as new variables are introduced'''
        self.sigma = self.Nb * self.c / (np.pi*self.R)
        self.v_tip = self.Omega*self.R
        self.rho = self.get_rho(self.h)
        self.vi_hov = (np.sqrt(self.W/(2*self.rho*np.pi*self.R**2)))
        self.A = self.R/self.c
        self.cbar = self.c/self.R
        self.v_tip = self.v_tip + v_air
        self.bl = self.get_bl(self.v_tip)
        self.myu_c = self.v_tip/v_air

    def get_vihov(self, R):
        return np.sqrt(self.W/(2*self.rho*np.pi*R**2))
        
    def setname(self, name):
        self.name = name
        
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

    def setEngine(self, engine):
        self.engine_choice = engine
        self.P_a = engine.determine_power_at_altitude(self.h)
        self.P_e = engine.P_sl
        
    def setTwinEngine(self, engine):
        self.engine_choice = engine
        self.P_a = 2*engine.determine_power_at_altitude(self.h)
        self.P_e = 2*engine.P_sl

    def plot_val(self, x, y, title=None, xlabel=None, ylabel=None, get_min=False, xlim=None, ylim=None, fname=None):
        '''Plot x against y'''
        plt.figure()
        plt.grid()
                    
        if not xlim: xlim = (0, max(x))
        if not ylim: ylim = (0, max(y))
        plt.xlim(xlim)    
        plt.ylim(ylim)
        
        if get_min:
            xmin, ymin = self.get_min(x, y)
            plt.axhline(ymin, 0, xmin/xlim[1], label='y= %s' % str(ymin), linestyle='dotted')
            plt.axvline(xmin, 0, ymin/ylim[1], label='x= %s' % str(xmin), linestyle='dotted')
        
        if title: plt.title(title)
        if xlabel: plt.xlabel(xlabel)
        if ylabel: plt.ylabel(ylabel)
        plt.plot(x, y)
        
        plt.legend(loc='upper right')
        if fname: plt.savefig(fname)
#        plt.show()

    def plot_vals(self, x, ys, labels=None, title=None, xlabel=None, ylabel=None, get_min=False, xlim=None, ylim=None, fname=None):
        '''Plot x against ys'''
        plt.figure()
        plt.grid()
                    
        for ind, y in enumerate(ys):
            if xlim: plt.xlim(xlim)    
            if ylim: plt.ylim(ylim)
            
            if get_min:
                xmin, ymin = self.get_min(x, y)
                plt.axhline(ymin, 0, xmin/xlim[1], label='y%i= %s' % (ind, str(ymin)), linestyle='dotted')
                plt.axvline(xmin, 0, ymin/ylim[1], label='x%i= %s' % (ind, str(xmin)), linestyle='dotted')
            
            if not labels:
                label='Curve%i' % ind
            else:
                label=labels[ind]

            plt.plot(x, y, label=label)
        
        if title: plt.title(title)
        if xlabel: plt.xlabel(xlabel)
        if ylabel: plt.ylabel(ylabel)
        plt.legend(loc='upper right')
        if fname: plt.savefig(fname)
#        plt.show()
        
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
        plt.plot((0,v_B), (0, Psfr), color='purple', linestyle='dotted')
        plt.axvline(v_B, 0, Psfr/Pto.max()/yscale, label='vsfr: %i m/s' % int(v_B), color='purple', linestyle='dashed')
        plt.axvline(v_A, 0, Pmin/Pto.max()/yscale, label='vpmin: %i m/s' % int(v_A), color='brown', linestyle='dashed')
        plt.axhline(Pmin, 0, v_A/v.max()/xscale, label='Pmin: %i kW' % int(Pmin), color='brown', linestyle='dotted')
        plt.legend(loc='lower right')
        plt.tight_layout()
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
        # minimum power airspeed
        p_min = P_to.min(axis=0) 
        v_A = v[np.where(P_to==p_min)]
        
        # minimum fuel flow airspeed (maximum range airspeed)
        # method: determine slope at all points on graph, determine slope of line between origin and said point, then find which one is the closest
        slope_power_to = np.gradient(P_to, axis=0)  # identical to dydx, but also includes first and last
        slope_origin_to_point = P_to/v  # identical to m, but also includes first and last
        difference = np.abs((slope_power_to - slope_origin_to_point)/slope_origin_to_point)
        closest_fit = difference.min(axis=0)
        closest_loc = np.where(difference==closest_fit)   
        
        v_B = v[closest_loc]
        p_sfr = P_to[closest_loc]
        
        # maximum airspeed
        equal_power = np.where(P_to>self.P_a)[0]
        if len(equal_power) > 0:
            v_C = v[equal_power][0]  # does not work when hover is critical, probably does not work at all with 2D
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
        self.calc_params()
        
        # velocity of the air seen by the heli
        Pexcess = 0 if P_to0 is None else self.P_a*1000 - P_to0
        v_cl_uncor = Pexcess/self.W  # currently based on energy analysis from navy test pilot manual
#        v_cl = v_cl_uncor/self.vi_hov + (1./(v_cl_uncor+self.vi_hov))  # from momentum equation in naval test pilot manual, still not sold on this https://apps.dtic.mil/dtic/tr/fulltext/u2/a061671.pdf
        v_cl = v_cl_uncor        
        v_inf = np.sqrt(v_air**2 + v_cl**2)
        
        # angle of the air seen by the rotor
        myu = v_air/(self.Omega*self.R)
        P_par = self.CDS*0.5*self.rho*v_inf**3
        alpha_tpp = -(Pexcess+P_par)/self.W/v_inf if P_to0 is not None else 0  # cunha 6
        
        vi_cr = self.vi_hov**2 / np.sqrt((v_inf*np.cos(alpha_tpp))**2+(v_inf*np.sin(alpha_tpp)+self.vi_hov)**2)  # from cunha lecture 6; vi_cr is induced velocity in cruise
        P_i = self.W*vi_cr
        P_pd = self.sigma*self.CDp_bar/8 * self.rho * ((self.Omega*self.R)**3) * np.pi*self.R**2 * (1+3*myu**2)  # to be given from Aero department, note: highly sensity to blade radius (R^5)!
        
        P_to = P_i + P_pd + P_par
        
        return (v_air, v_cl), (P_to, P_i, P_pd, P_par)
    
    def powerCurve(self, vs=np.arange(0, 105, 1), P_to0=None, figtitle='Power requirements for a rotorcraft in level forward flight', fname=None, select_engine=True):
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
        if select_engine:
            self.select_engine(power_vals[0])
        self.gen_pcurve(vs, *power_vals, figtitle=figtitle, fname=fname)
        return v_vals, power_vals
    
    def select_engine(self, P_to):
        '''
        Using the EngineSpider object in engine.py, iterate through all engine possibilities
        and select the "best" one. The best one is (currently) determined by whichever 
        engine provides the lowest power at altitude that meets hover power requirements.
        This is being updated to comply with FAR regulations for twin-engine rotorcraft
        Regulation 1: Take-off at 6400m
        Regulation 2: Vertical climb of 0.76m/s at 6705m
        
        Parameters
        ----------
        P_to : array_like
            An array of power values in kW, with P_to[0] being power at hover
        
        Returns
        -------
        Engine_best : engine.Engine
            The engine with the available power closest to the hover requirements 
            that meets them
        Engine_alternatives : list; engine.Engine
            All engines that meet the hover power requirements
        '''
        far_regulation_TO = copy.copy(self)
        far_regulation_TO.setname('far-regulation-1, TO')
        far_regulation_TO.seth(6400)
        farTO_v, farTO_power = far_regulation_TO.determineP_to(np.linspace(1, 100))#, P_to0=P_to[0])        
        
        engine_best, engine_alternatives = engine.engineoptions.select_engine(farTO_power[0][0], far_regulation_TO.h)
        self.setEngine(engine_best)

        far_regulation_vcl = copy.copy(self)
        far_regulation_vcl.setname('far-regulation-2, vcl')
        far_regulation_vcl.seth(6705)
        farvcl_v, farvcl_power = far_regulation_vcl.determineP_to(np.linspace(1, 100))
        farvcl_v_cl, farvcl_power_cl = far_regulation_vcl.determineP_to(np.linspace(1, 100), P_to0=farvcl_power[0])
                    
        vcl_req = 0.76  # m/s
        if farvcl_v_cl[1][0] < vcl_req:
#            print('2nd FAR req')
            vclreq_power = vcl_req / (far_regulation_vcl.vi_hov - vcl_req*far_regulation_vcl.vi_hov)*far_regulation_vcl.W/9.81 + farvcl_power[0][0]
            engine_best, engine_alternatives = engine.engineoptions.select_engine(vclreq_power, far_regulation_vcl.h)
            
        self.setTwinEngine(engine_best)
        self.engine_alternatives = engine_alternatives 
        return engine_best, engine_alternatives
            
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
        self.calc_params(v_air=0)
        vs, P_hovs = self.determineP_to(v_air=0)
        
        return vs, (Rs, P_hovs)
    
    def idealPhovafoc(self, cs=np.linspace(0.1, 0.6)):
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
        self.calc_params(v_air=0)
        vs, P_hovs = self.determineP_to(v_air=0)
        
        return vs, (cs, P_hovs)
    
    def idealPhovafoOmega(self, Omegas=np.linspace(15, 40)):
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
        self.calc_params(v_air=0)
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
        
    def v_charsafoc(self, vs, cs=np.linspace(0.1, 0.6, 50)):
        '''
        same as v_charsafoR, but wrt chord c
        '''
        cs, vs = np.meshgrid(cs, vs)
        self.setc(cs)
        self.calc_params()
        vs, Ps = self.determineP_to(v_air=vs)
        vchars = self.determineV_chars(v=vs[0], P_to=Ps[0])
        return vs, (cs, vchars)

    def v_charsafoOmega(self, vs, Omegas=np.linspace(15, 40, 50)):
        '''
        same as v_charsafoR, but wrt angular velocity Omega
        '''
        Omegas, vs = np.meshgrid(Omegas, vs)
        self.setOmega(Omegas)
        self.calc_params()
        vs, Ps = self.determineP_to(v_air=vs)
        vchars = self.determineV_chars(v=vs[0], P_to=Ps[0])
        return vs, (Omegas, vchars)
        
    def Mtip_ne_constraint(self):
        #vtip < 0.92 mach @ vne @ 8950m
        vtip_ne = self.v_tip + self.vne
        Mtip_ne = vtip_ne/self.get_a()
        return Mtip_ne < 0.92
        
    def Mtip_cr_constraint(self):
        #Mtip @ cruise @ 8950m < 0.85
        vtip_cr = self.v_tip + self.v_cr
        Mtip_cr = vtip_cr/self.get_a()
        return Mtip_cr < 0.85
    
    def myu_constraint(self):
        #myu < 0.45 in cruise
        myu = self.v_tip / self.v_cr
        return myu < 0.45
        
    def bl_constraint(self):
        #bl < 0.12
        v_tip = self.v_cr+self.v_tip
        s = 0.6
        bl = s*self.W*9.81 / (self.c*self.Nb*self.R*v_tip**2)
        return bl < 0.12

class Coaxheli(SRheli):
    
    configuration = 'coaxial-rotor'
    iscoax = True
    zD = 0

    def __init__(self, *args, **kwargs):
        super(Coaxheli, self).__init__(*args, **kwargs)
        
        if 'zD' in kwargs.keys():
            self.set_zD(kwargs['zD'])
        
        if 'convert' in kwargs.keys() and kwargs['convert']:
            self.convert()

    def set_zD(self, zDs):
        self.zD = zDs
        
    def convert(self): 
        '''
        Converts the given single-rotor data into the equivalent (solidity) coaxial rotor data
        '''
#        self.Nr = 2 
#        self.Nb = self.Nb/2
        self.R = self.R*np.sqrt(2)
        self.c = self.c*np.sqrt(2)
        self.Omega = self.Omega/np.sqrt(2)
    
    def reconvert(self):
        '''
        Takes the (current) coaxial rotor parameters, and returns them into their single-rotor equivalent sizes
        '''        
        SR_Nr = 1
        SR_Nb = self.Nb*2
        SR_R = self.R/np.sqrt(2)
        SR_c = self.c/np.sqrt(2)
        SR_Omega = self.Omega*np.sqrt(2)
        
        return SR_Nr, SR_Nb, SR_R, SR_c, SR_Omega
