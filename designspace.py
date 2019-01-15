# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 11:35:24 2018

@author: Matt
"""
from __future__ import division 
import numpy as np
import matplotlib.pyplot as plt

class DesignSpace(object):
    
    def __init__(self, leg, h, v_cr):
        '''
        Initialise and object to check the design space for the HAMRAC.
        Values given are fixed for the current mission.
        Statistical boundaries and constraints are given by Leon.
        
        Parameters
        ----------
        leg : str
            The name of the leg for the design space to analyse
        h : int
            The altitude of the
        '''
        self.leg = leg
        self.v_cr = v_cr
        self.vmax = v_cr*1.15 
        self.vne = self.vmax*1.1
        self.Nb = 6
        self.W = 2487*9.81
        self.h = h

        self.statistical_boundaries = {
            'c': [0.19, 0.45],
            'R': [3.5, 7],
            'Omega': [25, 42],
            'Nb': [3, 4],
            'cbar': [0.04, 0.08],  # c/R
            'A': [12.5, 25],  # R/c
            'sigma': [0.03*2, 0.12*2],  # Nb*c/pi/R
            'bl': [0.07, 0.1],  #  s*Mto*g/sigma/c/Nb/R/vtip
            'vtip': [156, 193],  # Omega * R
            'myu': [0.25*np.sqrt(2), 0.45],  # vair / vtip
        }       

    def Mtip_cr_constraint(self, R, c, Omega):
        '''Tip Mach at cruise'''
        vtip_cr = Omega*R + self.v_cr
        Mtip_cr = vtip_cr/self.get_a(self.h)
        return Mtip_cr
       
    def Mtip_ne_constraint(self, R, c, Omega):
        '''Tip Mach at never exceed speed'''
        vtip_ne = Omega*R + self.vne
        Mtip_ne = vtip_ne/self.get_a(self.h)
        return Mtip_ne
    
    def myu_constraint(self, R, c, Omega):
        '''Advance ratio, v_cr/v_tip'''
        myu = self.v_cr / (self.v_cr+Omega*R) 
        return myu 
        
    def bl_constraint(self, R, c, Omega):
        '''blade loading'''
#        v_tip = Omega*R +self.v_cr
#        s = 0.6
#        sigma = self.Nb * c / np.pi / R
        Cts = self.W / (self.get_rho(self.h) * self.Nb * c * R * (Omega*R)**2)
#        Cts = s*self.W / (sigma * c * self.Nb * R * v_tip**2)
        return Cts
        
    def det_constraints(self, R, c, Omega):
        ''' Mtip @ cruise @ 8950m < 0.85 '''
        ''' vtip < 0.92 mach @ vne @ 8950m '''
        ''' myu < 0.45 in cruise '''
        ''' bl < 0.12 '''
        Mtip_cr = self.Mtip_cr_constraint(R, c, Omega)
        Mtip_ne = self.Mtip_ne_constraint(R, c, Omega)
        myu = self.myu_constraint(R, c, Omega)
        bl = self.bl_constraint(R, c, Omega)
        return (Mtip_cr, Mtip_ne, myu, bl), (Mtip_cr < 0.85, Mtip_ne < 0.92, myu < 0.45, bl < 0.12)

    def get_a(self, h):
        return np.sqrt(1.4 * 287 * self.get_T(h))
        
    @staticmethod
    def get_rho(h):
        return 1.225*288.16/(15+271.16-1.983*h/304.8)*(1-(0.001981*h/0.3048)/288.16)**5.256

    @staticmethod
    def get_T(h):
            return (15+271.16-1.983*h/304.8)
          
    def get_Rafovtip_Omega(self):
        '''
        Determine values for R as a function of v_tip and Omega
        Fix vtip to minimum (vtip1), and determine R1, R2 for Omegamin, Omegamax
        Fix vtip to maximum (vtip2), and determine R3, R4 for Omegamin, Omegamax
        '''
        vtip1, vtip2 = self.statistical_boundaries['vtip']
        Omega1, Omega2 = self.statistical_boundaries['Omega']
       
        R1 = vtip1/Omega1
        R2 = vtip1/Omega2
        R3 = vtip2/Omega2
        R4 = vtip2/Omega1
        
        return R1, R2, R3, R4
        
    def get_Omegaafovtip_R(self):
        '''
        Determine values for Omega as a function of v_tip and R
        Fix vtip to minimum (vtip1), and determine Omega1, Omega2 for Rmin, Rmax
        Fix vtip to maximum (vtip2), and determine Omega3, Omega4 for Rmin, Rmax
        '''
        vtip1, vtip2 = self.statistical_boundaries['vtip']
        R1, R2 = self.statistical_boundaries['R']
        
        Omega1 = vtip1/R1
        Omega2 = vtip1/R2
        Omega3 = vtip2/R2
        Omega4 = vtip2/R1
        
        return Omega1, Omega2, Omega3, Omega4
        
    def get_Rafomyu_Omega(self):
        ''' same as get_Rafovtip_Omega but for advance ratio myu '''
        myu1, myu2 = self.statistical_boundaries['myu']
        Omega1, Omega2 = self.statistical_boundaries['Omega']
        
        R1 = self.v_cr / Omega1 / myu1
        R2 = self.v_cr / Omega1 / myu2
        R3 = self.v_cr / Omega2 / myu2
        R4 = self.v_cr / Omega2 / myu1
        
        return R1, R2, R3, R4
        
    def get_Omegaafomyu_R(self):
        ''' same as get_Omegaafovtip_R but for advance ratio myu '''
        myu1, myu2 = self.statistical_boundaries['myu']        
        R1, R2 = self.statistical_boundaries['R']
        
        Omega1 = self.v_cr / R1 / myu1
        Omega2 = self.v_cr / R1 / myu2
        Omega3 = self.v_cr / R2 / myu2
        Omega4 = self.v_cr / R2 / myu1
        
        return Omega1, Omega2, Omega3, Omega4
        
    def get_Rafosigma_c(self):
        '''
        Determine values for R as a function of sigma and c
        Fix sigma to minimum (sigma), and determine R1, R2 for cmin, cmax
        Fix sigma to maximum (sigma), and determine R3, R4 for cmin, cmax
        '''
        c1, c2 = self.statistical_boundaries['c']
        sigma1, sigma2 = self.statistical_boundaries['sigma']
        
        R1 = self.Nb*c1/np.pi/sigma1
        R2 = self.Nb*c1/np.pi/sigma2
        R3 = self.Nb*c2/np.pi/sigma2
        R4 = self.Nb*c2/np.pi/sigma1
        
        return R1, R2, R3, R4

    def get_cafosigma_R(self):
        '''
        Determine values for c as a function of sigma and R
        Fix sigma to minimum (sigma), and determine c1, c2 for Rmin, Rmax
        Fix sigma to maximum (sigma), and determine c3, c4 for Rmin, Rmax
        '''
        R1, R2 = self.statistical_boundaries['R']
        sigma1, sigma2 = self.statistical_boundaries['sigma']
        
        c1 = R1*np.pi*sigma1/self.Nb
        c2 = R1*np.pi*sigma2/self.Nb
        c3 = R2*np.pi*sigma2/self.Nb
        c4 = R2*np.pi*sigma1/self.Nb

        return c1, c2, c3, c4
        
    def get_RafoA_c(self):
        ''' same as get_Rafosigma_c but for aspect ratio A '''
        c1, c2 = self.statistical_boundaries['c']
        A1, A2 = self.statistical_boundaries['A']
        
        R1 = A1*c1
        R2 = A2*c1
        R3 = A2*c2
        R4 = A1*c2
        
        return R1, R2, R3, R4
        
    def get_cafoA_R(self):
        ''' same as get_cafosigma_R but for aspect ratio A '''
        R1, R2 = self.statistical_boundaries['R']
        A1, A2 = self.statistical_boundaries['A']
        
        c1 = R1/A1
        c2 = R1/A2
        c3 = R2/A2
        c4 = R2/A1
        
        return c1, c2, c3, c4 
        
    def get_OmegaafoCts_Rc(self, R, c):
        # take R-c design space, use to determine Omega = f(Ct/sigma, R, c)
        Cts1, Cts2 = self.statistical_boundaries['bl']
        
        Omega1 = np.sqrt(self.W/(R**3*self.get_rho(self.h)*2*self.Nb*c*Cts1)) 
        Omega2 = np.sqrt(self.W/(R**3*self.get_rho(self.h)*2*self.Nb*c*Cts2))
        
        return Omega1, Omega2
        
    def get_cafoCts_ROmega(self, R, Omega):
        # take R-Omega design space, use to determine Omega = f(Ct/sigma, R, Omega)
        Cts1, Cts2 = self.statistical_boundaries['bl']
        
        c1 = self.W/ self.get_rho(self.h)/self.Nb/R/Cts1/(Omega*R)**2               
        c2 = self.W/ self.get_rho(self.h)/self.Nb/R/Cts2/(Omega*R)**2        
        
        return c1, c2
        
    def determine_borders(self, Rv, Rm, Omegav, Omegam, Rs, RA, cs, cA):
        '''
        Determine the lines that will determine the borders of the design space.
        
        Parameters
        ----------
        Rv, Rm, Omegav, Omegam, Rs, RA, cs, cA: array_like
            A 1x4 array with the one axis coordinates of the box constraining the 
            design space. The other axis is the min/max allowable in the 
            statistical boundaries. 
            suffixes: v -> vtip, m -> myu, s -> sigma, A -> A, 1 -> minimum, 2 -> maximum
        
        
        Returns
        -------
        ((xsOv1, xsOv2, xsRv1, xsRv2, xsOm1, xsOm2, xsRm1, xsRm2), 
         (xsRs1, xsRs2, xsRA1, xsRA2, xscs1, xscs2, xscA1, xscA2),) : array_like
            The arrays constructed are 50 length between the minimum and maximum values
            obtained from each constraint (get_XXXafoX_XX methods).
        '''        
        R1, R2 = self.statistical_boundaries['R']
        Omega1, Omega2 = self.statistical_boundaries['Omega']
        
        # R-Omega relations
        xsOv1 = np.linspace(Omegav[0], Omegav[1], 50)
        xsOv2 = np.linspace(Omegav[3], Omegav[2], 50)

        xsOm1 = np.linspace(Omegam[1], Omegam[2])
        xsOm2 = np.linspace(Omegam[0], Omegam[3])
        
        xsRv1 = np.linspace(Rv[0], Rv[1], 50)
        xsRv2 = np.linspace(Rv[3], Rv[2], 50)
        
        xsRm1 = np.linspace(Rm[1], Rm[2], 50)
        xsRm2 = np.linspace(Rm[0], Rm[3], 50)

        # R-c relations
        xsRs1 = np.linspace(Rs[1], Rs[2], 50) # works
        xsRs2 = np.linspace(Rs[0], Rs[3], 50) # works
        
        xsRA1 = np.linspace(RA[0], RA[3], 50) # works
        xsRA2 = np.linspace(RA[1], RA[2], 50) # works
        
        xscs1 = np.linspace(cs[0], cs[3], 50) # works
        xscs2 = np.linspace(cs[1], cs[2], 50) # works

        xscA1 = np.linspace(cA[1], cA[0], 50) # works
        xscA2 = np.linspace(cA[2], cA[3], 50) # works

        return ((xsOv1, xsOv2, xsRv1, xsRv2, xsOm1, xsOm2, xsRm1, xsRm2), 
#                (Ovb1, Ovb2, Rvb1, Rvb2, Omb1, Omb2, Rmb1, Rmb2),
                (xsRs1, xsRs2, xsRA1, xsRA2, xscs1, xscs2, xscA1, xscA2),)
        
    def compare_statistics(self):
        Rmin, Rmax = self.statistical_boundaries['R']
        Omegamin, Omegamax = self.statistical_boundaries['Omega']
        cmin, cmax = self.statistical_boundaries['c']
        
        # R/Omega plot
        Rv1, Rv2, Rv3, Rv4 = self.get_Rafovtip_Omega()
        Omegav1, Omegav2, Omegav3, Omegav4 = self.get_Omegaafovtip_R()
        Rm1, Rm2, Rm3, Rm4 = self.get_Rafomyu_Omega()
        Omegam1, Omegam2, Omegam3, Omegam4 = self.get_Omegaafomyu_R()
        
        # R/c plot
        Rs1, Rs2, Rs3, Rs4 = self.get_Rafosigma_c()
        cs1, cs2, cs3, cs4 = self.get_cafosigma_R()
        RA1, RA2, RA3, RA4 = self.get_RafoA_c()
        cA1, cA2, cA3, cA4 = self.get_cafoA_R()
                
        return (((Rv1, Rv2, Rv3, Rv4), (Omegav1, Omegav2, Omegav3, Omegav4)),
                ((Rm1, Rm2, Rm3, Rm4), (Omegam1, Omegam2, Omegam3, Omegam4)),
                ((Rs1, Rs2, Rs3, Rs4), (cs1, cs2, cs3, cs4)),
                ((RA1, RA2, RA3, RA4), (cA1, cA2, cA3, cA4)),
                )
        
    def get_available_space(self, coordsX, coordsY, coordsZ, xsRO, xsRc):
        # rip python usage
        # some of the lines are still backwards, (currently) they do not constrain the design space at all
        '''
        Determine the design space of the HAMRAC.
        Checks each coordinate such that it is above the minimum value constraint under that value,
        and below the maximum value constraint under that value.
        R-Omega plot is x-y; R-c plot is x-z. y-z space is assumed and does not further constrain the space 
        Sadly, not applicable for real estate         
        
        Example: R=5, Omega=40. 
            Check that point (5,40) is to the right of line xsRv1 and to the left of line xsRv2
            xsRv1 is a series of R values, and is indexed by the col value, as the col value 
            checks where the xsRv1 line is in Y. This gives the R value (in x) at the desired Y index
            xsOv1 is a series of Omega values, and is indexed by the row value, as the row value
            checks where the xsOv1 line is in X. This gives the Omega value (in y) at the desired X index
            Roughly speaking, check the X coordinate for changes in Y (O[row]), and
            check the Y coordinate for changes in X (R[col]).            
            For 1 suffixes, R, Omega, c must be larger
            For 2 suffixes, R, Omega, c must be smaller
            
        Parameters
        ----------
        coordsX : array_like
            An array of the coordinates between Rmin and Rmax
        coordsY : array_like
            An array of the coordinates between Omegamin and Omegamax
        coordsZ : array_like
            An array of the coordinates between cmin and cmax.
        xsRO : array_like
            A series of arrays corresponding to the constraint coordinates 
            determined from self.determine_borders for the R-Omega plot
        xsRc : array_like
            A series of arrays corresponding to the constraint coordinates
            determined from self.determine_borders for the R-c plot
        '''        
        xsOv1, xsOv2, xsRv1, xsRv2, xsOm1, xsOm2, xsRm1, xsRm2 = xsRO
        xsRs1, xsRs2, xsRA1, xsRA2, xscs1, xscs2, xscA1, xscA2 = xsRc
        R1, R2 = self.statistical_boundaries['R']
        O1, O2 = self.statistical_boundaries['Omega']
        c1, c2 = self.statistical_boundaries['c']
        feasible_cells_RO = []
        feasible_cells_Rc = []
 
        
        for row in range(len(coordsX)):
            for col in range(len(coordsY)):
                R, O, c = coordsX[row], coordsY[col], coordsZ[col]
                
                if R > xsRv1[col] and R > xsRm1[col] and O > xsOv1[row] and O > xsOm1[row]:
                    if R < xsRv2[col]  and R < xsRm2[col] and O < xsOv2[row] and O < xsOm2[row]:
                        feasible_cells_RO.append([row, col])
                        
                if R > xsRs1[col] and R > xsRA1[col] and c > xscs1[row] and c > xscA1[row]:
                    if R < xsRs2[col] and R < xsRA2[col] and c < xscs2[row] and c < xscA2[row]:
                        feasible_cells_Rc.append([row, col])

        
        return feasible_cells_RO, feasible_cells_Rc
        
    def plot_boundaries(self, Rvs, Rms, Ovs, Oms, Rss, RAs, css, cAs, idealR=None, idealc=None, idealOmega=None):
        '''
        Determine the borders based on given parameters, then create a plot of them, then
        determine the available design space and create a scattergraph of them
        
        Parameters
        ----------
            Rvs, Rms, Ovs, Oms, Rss, RAs, css, cAs : array_like
                A series of arrays to form the boxes for the constraints. Each array
                is 4 in length, and represents the one axis coordinates of the box
                The other axis coordinates are set to min/max of the 'other' parameter
                Other ->: R-Omega; R-c
        '''        
        Rmin, Rmax = self.statistical_boundaries['R']
        Omegamin, Omegamax = self.statistical_boundaries['Omega']
        cmin, cmax = self.statistical_boundaries['c']
        coordsX = np.linspace(Rmin, Rmax)
        coordsY = np.linspace(Omegamin, Omegamax)
        coordsZ = np.linspace(cmin, cmax)
        
        ((Rv1, Rv2, Rv3, Rv4), (Omegav1, Omegav2, Omegav3, Omegav4),
         (Rm1, Rm2, Rm3, Rm4), (Omegam1, Omegam2, Omegam3, Omegam4),
         (Rs1, Rs2, Rs3, Rs4), (cs1, cs2, cs3, cs4),
         (RA1, RA2, RA3, RA4), (cA1, cA2, cA3, cA4)) = (Rvs, Ovs, Rms, Oms, Rss, css, RAs, cAs)
         
        xsRO, xsRc = self.determine_borders(Rvs, Rms, Ovs, Oms, Rss, RAs, css, cAs)
        xsOv1, xsOv2, xsRv1, xsRv2, xsOm1, xsOm2, xsRm1, xsRm2 = xsRO
        xsRs1, xsRs2, xsRA1, xsRA2, xscs1, xscs2, xscA1, xscA2 = xsRc
        
        dspaceRO, dspaceRc= self.get_available_space(coordsX, coordsY, coordsZ, xsRO, xsRc)
        dspaceRO = np.array(dspaceRO)
        dspaceRc = np.array(dspaceRc)
        
        Ropt = coordsX[dspaceRO[:,0]]
        Oopt = coordsY[dspaceRO[:,1]]
        R2opt = coordsX[dspaceRc[:,0]]
        copt = coordsZ[dspaceRc[:,1]]
#        
        
#        Cts1 = self.get_OmegaafoCts_Rc(blmin[:,0], blmin[:,1])
        Cts2 = self.get_cafoCts_ROmega(Ropt, Oopt)
        
        Rset = set(Ropt)
        Rset = sorted(Rset)  
        Rlist = list(Ropt)
        Cts21 = list(Cts2[0])
        Cts22 = list(Cts2[1])
        Cts2upper = []
        Cts2lower = []
            
        for Rval in Rset:
            # list goes the wrong way, need to find last instance of Rval, but cannot list.reverse()
            Rindleft = Rlist.index(Rval)
            Cts2upper.append(Cts22[Rindleft])

        fig = plt.figure()
#        axRO = fig.add_subplot(111, projection='3d')
        axRO = fig.add_subplot(111)
        plt.title('Design Space at %s for the HAMRAC R-Omega parameters' % self.leg)
        plt.grid(which='major')
        plt.xlabel('R [m]')
        plt.ylabel('Omega [rad/s]')
        plt.xlim(3, 9)
#        plt.ylim(22, 50)
        
        axRO.plot((Rmin, Rmax, Rmax, Rmin, Rmin), (Omegamin, Omegamin, Omegamax, Omegamax, Omegamin), label='R, Omega constraints')
        axRO.plot((Rv1, Rv2, Rv3, Rv4, Rv1), (Omegamin, Omegamax, Omegamax, Omegamin, Omegamin), label='vtip, R constraints')  
        axRO.plot((Rmin, Rmax, Rmax, Rmin, Rmin), (Omegav1, Omegav2, Omegav3, Omegav4, Omegav1), label='vtip, Omega constraints')        
        axRO.plot((Rm1, Rm2, Rm3, Rm4, Rm1), (Omegamin, Omegamin, Omegamax, Omegamax, Omegamin), label='myu, R constraints')  
        axRO.plot((Rmin, Rmin, Rmax, Rmax, Rmin), (Omegam1, Omegam2, Omegam3, Omegam4, Omegam1), label='myu, Omega constraints')        
        
        plt.scatter(Ropt, Oopt, c='#fff9f6', marker='o', s=2)  
        if idealR and idealOmega:
            plt.scatter(idealR, idealOmega, marker='x', c='y', s=100)
#        plt.scatter(R2opt, Cts1[0])
#        plt.scatter(R2opt, Cts1[1])

        plt.legend()
        plt.savefig('mdo/desspaceRO_%s.png' % self.leg)
        plt.show()
#        
        fig = plt.figure()
        axRc = fig.add_subplot(111)#, projection='3d')
        plt.title('Design Space at %s for the HAMRAC R-c parameters' % self.leg)
        plt.grid(which='major')
        plt.xlabel('R [m]')
        plt.ylabel('c [m]')
        plt.xlim(3, 9)
#        plt.ylim(0.1, 0.5)

        axRc.plot((Rmin, Rmax, Rmax, Rmin, Rmin), (cmin, cmin, cmax, cmax, cmin))
        axRc.plot((Rs1, Rs2, Rs3, Rs4, Rs1), (cmin, cmin, cmax, cmax, cmin), label='sigma, R constraints')
        axRc.plot((Rmin, Rmin, Rmax, Rmax, Rmin), (cs1, cs2, cs3, cs4, cs1), label='sigma, c constraints')        
        axRc.plot((RA1, RA2, RA3, RA4, RA1), (cmin, cmax, cmax, cmin, cmin), label='A, R constraints')
        axRc.plot((Rmin, Rmin, Rmax, Rmax, Rmin), (cA1, cA2, cA3, cA4, cA1), label='A, c constraints')        
        
        plt.scatter(R2opt, copt, c=[(1, 1, 1, 0) for i in range(len(copt))], s=2)  
        if idealR and idealc:
            plt.scatter(idealR, idealc, marker='x', c='y', s=200)
#        plt.scatter(blmin[:,0], Cts2[0])
        plt.plot(Ropt, Cts2[1])
        plt.plot(Ropt, Cts2[0])
        plt.plot(Rset, Cts2upper)
        
        plt.legend(loc='lower right')
        plt.savefig('mdo/desspaceRc_%s.png' % self.leg)
        plt.show()
        
rescue = DesignSpace('rescue', 8950, 50) # 120 kts
far1 = DesignSpace('far1', 6400, 72)
far2 = DesignSpace('far2', 6705, 72)
takeoff1 = DesignSpace('takeoff intl airport', 1402, 72)
cruise1 = DesignSpace('cruise1', 3930, 80) # 150 kts
climb1 = DesignSpace('climb1', 6200, 40)
takeoff2 = DesignSpace('takeoff refuelling', 3780, 72)


data = rescue.compare_statistics()
far1space = far1.compare_statistics()
far2space = far2.compare_statistics()
takeoff1space = takeoff1.compare_statistics()
cruise1space = cruise1.compare_statistics()
climb1space = climb1.compare_statistics()
takeoff2space = takeoff2.compare_statistics()

R = 5.5
c = 0.36
Omega = 28
Omega2 = 35

#rescue.plot_boundaries(data[0][0], data[1][0], data[0][1], data[1][1], data[2][0], data[3][0], data[2][1], data[3][1], idealR=R, idealc=c, idealOmega=Omega)
#far1.plot_boundaries(far1space[0][0], far1space[1][0], far1space[0][1], far1space[1][1], far1space[2][0], far1space[3][0], far1space[2][1], far1space[3][1], idealR=idealR, idealc=idealc, idealOmega=idealOmega)
#far2.plot_boundaries(far2space[0][0], far2space[1][0], far2space[0][1], far2space[1][1], far2space[2][0], far2space[3][0], far2space[2][1], far2space[3][1], idealR=idealR, idealc=idealc, idealOmega=idealOmega)
#takeoff1.plot_boundaries(takeoff1space[0][0], takeoff1space[1][0], takeoff1space[0][1], takeoff1space[1][1], takeoff1space[2][0], takeoff1space[3][0], takeoff1space[2][1], takeoff1space[3][1], idealR=idealR, idealc=idealc, idealOmega=idealOmega)
#takeoff2.plot_boundaries(takeoff2space[0][0], takeoff2space[1][0], takeoff2space[0][1], takeoff2space[1][1], takeoff2space[2][0], takeoff2space[3][0], takeoff2space[2][1], takeoff2space[3][1], idealR=idealR, idealc=idealc, idealOmega=idealOmega)

# rotor 2: at takeoff+cruise+refuel
#f12bounds, f12validity = far1.det_constraints(*ideals2)
#f22bounds, f22validity = far2.det_constraints(*ideals2)
to2bounds, to2validity = takeoff2.det_constraints(R, c, Omega)
c1bounds, c1validity = cruise1.det_constraints(R, c, Omega)
to1bounds, to1validity = takeoff1.det_constraints(R, c, Omega)

#blade 2: at refuel+FAR+rescue
rbounds, rvalidity = rescue.det_constraints(R, c, Omega2)
f2bounds, f2validity = far2.det_constraints(R, c, Omega2)
f1bounds, f1validity = far1.det_constraints(R, c, Omega2)
to22bounds, to22validity = takeoff2.det_constraints(R, c, Omega2)
#c1bounds, c1validity = cruise1.det_constraints(*ideals2)

   
vals1 = [('takeoff1', takeoff1.h, to1bounds, to1validity),
         ('cruise1', cruise1.h, c1bounds, c1validity), 
         ('takeoff2', takeoff2.h, to2bounds, to2validity), ]

vals2 = [('takeoff2', takeoff2.h, to22bounds, to22validity),
         ('far1', far1.h, f1bounds, f1validity), 
         ('far2', far2.h, f2bounds, f2validity),
         ('rescue', rescue.h, rbounds, rvalidity), ]
     

f = open('mdo/designspaceconstraintsresults.txt', 'w')

paramstr1 = '\nvals1: R=%f, c=%f, Omega=%f\n' % (R, c, Omega)
print(paramstr1)
f.write(paramstr1)

for v in vals1:
    output = '''%s (h=%i):\n----------------
    Mtip_cr: %f - %s
    Mtip_ne: %f - %s
    myu: %f - %s
    bl: %f - %s\n''' % (v[0], v[1], v[2][0], v[3][0], v[2][1], v[3][1], v[2][2], v[3][2], v[2][3], v[3][3])
    print(output)
    f.write(output)

paramstr2 = '\nvals2: R=%f, c=%f, Omega=%f\n' % (R, c, Omega2)
print(paramstr2)
f.write(paramstr2)

for v in vals2:
    output = '''%s (h=%i):\n----------------
    Mtip_cr: %f - %s
    Mtip_ne: %f - %s
    myu: %f - %s
    bl: %f - %s\n\n''' % (v[0], v[1], v[2][0], v[3][0], v[2][1], v[3][1], v[2][2], v[3][2], v[2][3], v[3][3])
    print(output)
    f.write(output)


f.close() 