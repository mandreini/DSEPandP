# -*- coding: utf-8 -*-
"""
Created on Mon Dec 17 11:35:24 2018

@author: Matt
"""
from __future__ import division 
import numpy as np
import matplotlib.pyplot as plt

class DesignSpace(object):
    
    def __init__(self):
        self.v_cr = 140/1.9
        self.Nb = 3

        self.statistical_boundaries = {
            'c': [0.19, 0.45],
            'R': [3.5, 7],
            'Omega': [25, 42],
            'Nb': [3, 4],
            'cbar': [0.04, 0.08],  # c/R
            'A': [12.5, 25],  # R/c
            'sigma': [0.03, 0.12],  # Nb*c/pi/R
            'bl': [0.07, 0.1],  #  W/pi/R/R/rho/vtip/vtip/sigma 
            'vtip': [180, 220],  # Omega * R + vair
            'myu': [0.25, 0.4],  # vair / vtip
        }        
        
    def calca(self,):
        return np.sqrt(self.gamma, self.R, self.get_T(self.h))
        
    
    def get_Rafovtip_Omega(self):
        vtip1, vtip2 = self.statistical_boundaries['vtip']
        Omega1, Omega2 = self.statistical_boundaries['Omega']
       
        R1 = vtip1/Omega1
        R2 = vtip1/Omega2
        R3 = vtip2/Omega2
        R4 = vtip2/Omega1
        
        return R1, R2, R3, R4
        
    def get_Omegaafovtip_R(self):
        vtip1, vtip2 = self.statistical_boundaries['vtip']
        R1, R2 = self.statistical_boundaries['R']
        
        Omega1 = vtip1/R1
        Omega2 = vtip1/R2
        Omega3 = vtip2/R2
        Omega4 = vtip2/R1
        
        return Omega1, Omega2, Omega3, Omega4
        
    def get_Rafomyu_Omega(self):
        myu1, myu2 = self.statistical_boundaries['myu']
        Omega1, Omega2 = self.statistical_boundaries['Omega']
        
        R1 = self.v_cr / Omega1 / myu1
        R2 = self.v_cr / Omega1 / myu2
        R3 = self.v_cr / Omega2 / myu2
        R4 = self.v_cr / Omega2 / myu1
        
        return R1, R2, R3, R4
        
    def get_Omegaafomyu_R(self):
        myu1, myu2 = self.statistical_boundaries['myu']        
        R1, R2 = self.statistical_boundaries['R']
        
        Omega1 = self.v_cr / R1 / myu1
        Omega2 = self.v_cr / R1 / myu2
        Omega3 = self.v_cr / R2 / myu2
        Omega4 = self.v_cr / R2 / myu1
        
        return Omega1, Omega2, Omega3, Omega4
        
    def get_Rafosigma_c(self):
        c1, c2 = self.statistical_boundaries['c']
        sigma1, sigma2 = self.statistical_boundaries['sigma']
        
        R1 = self.Nb*c1/np.pi/sigma1
        R2 = self.Nb*c1/np.pi/sigma2
        R3 = self.Nb*c2/np.pi/sigma2
        R4 = self.Nb*c2/np.pi/sigma1
        
        return R1, R2, R3, R4

    def get_cafosigma_R(self):
        R1, R2 = self.statistical_boundaries['R']
        sigma1, sigma2 = self.statistical_boundaries['sigma']
        
        c1 = R1*np.pi*sigma1/self.Nb
        c2 = R1*np.pi*sigma2/self.Nb
        c3 = R2*np.pi*sigma2/self.Nb
        c4 = R2*np.pi*sigma1/self.Nb

        return c1, c2, c3, c4
        
    def get_RafoA_c(self):
        c1, c2 = self.statistical_boundaries['c']
        A1, A2 = self.statistical_boundaries['A']
        
        R1 = A1*c1
        R2 = A1*c2
        R3 = A2*c2
        R4 = A2*c1
        
        return R1, R2, R3, R4
        
    def get_cafoA_R(self):
        R1, R2 = self.statistical_boundaries['R']
        A1, A2 = self.statistical_boundaries['A']
        
        c1 = R1/A1
        c2 = R1/A2
        c3 = R2/A2
        c4 = R2/A1
        
        return c1, c2, c3, c4 
        
    def determine_borders(self, Rv, Rm, Omegav, Omegam, Rs, RA, cs, cA):
        # v -> vtip, m -> myu 
        R1, R2 = self.statistical_boundaries['R']
        Omega1, Omega2 = self.statistical_boundaries['Omega']
        xs = np.linspace(R1, R2)
        
        # R-Omega relations
        xsOv1 = np.linspace(Omegav[1], Omegav[2], 50)
        xsOv2 = np.linspace(Omegav[3], Omegav[0], 50)

        xsOm1 = np.linspace(Omegam[1], Omegam[2])
        xsOm2 = np.linspace(Omegam[3], Omegam[0])
        
        xsRv1 = np.linspace(Rv[1], Rv[2], 50)
        xsRv2 = np.linspace(Rv[0], Rv[3], 50)
        
        xsRm1 = np.linspace(Rm[1], Rm[2], 50)
        xsRm2 = np.linspace(Rm[3], Rm[0], 50)

        # R-c relations
        xsRs1 = np.linspace(Rs[0], Rs[1], 50)
        xsRs2 = np.linspace(Rs[2], Rs[3], 50)
        
        xsRA1 = np.linspace(RA[0], RA[1], 50)
        xsRA2 = np.linspace(RA[2], RA[3], 50)
        
        xscs1 = np.linspace(cs[0], cs[1], 50)
        xscs2 = np.linspace(cs[2], cs[3], 50)

        xscA1 = np.linspace(cA[1], cA[0], 50)
        xscA2 = np.linspace(cA[3], cA[2], 50)

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
        
#        fig = plt.figure()
##        axRO = fig.add_subplot(111, projection='3d')
#        axRO = fig.add_subplot(111)
#        plt.title("Design Space for the HAMRAC")
#        plt.grid(which='major')
#        plt.xlabel('R [m]')
#        plt.ylabel('Omega [rad/s]')
##        plt.xlim(3, 8.5)
##        plt.ylim(22, 50)
#        
#        axRO.plot((Rmin, Rmax, Rmax, Rmin, Rmin), (Omegamin, Omegamin, Omegamax, Omegamax, Omegamin), label='R, Omega constraints')
#        axRO.plot((Rv1, Rv2, Rv3, Rv4, Rv1), (Omegamin, Omegamin, Omegamax, Omegamax, Omegamin), label='vtip, R constraints')  
#        axRO.plot((Rmin, Rmin, Rmax, Rmax, Rmin), (Omegav1, Omegav2, Omegav3, Omegav4, Omegav1), label='vtip, Omega constraints')        
#        axRO.plot((Rm1, Rm2, Rm3, Rm4, Rm1), (Omegamin, Omegamin, Omegamax, Omegamax, Omegamin), label='myu, R constraints')  
#        axRO.plot((Rmin, Rmin, Rmax, Rmax, Rmin), (Omegam1, Omegam2, Omegam3, Omegam4, Omegam1), label='myu, Omega constraints')        
##        
#        plt.legend()
#        plt.show()
#        
#        fig = plt.figure()
#        axRc = fig.add_subplot(111)#, projection='3d')
#        plt.title("Design Space for the HAMRAC")
#        plt.grid(which='major')
#        plt.xlabel('R [m]')
#        plt.ylabel('c [m]')
##        plt.xlim(3, 10)
##        plt.ylim(0.1, 0.5)
#
#        axRc.plot((Rmin, Rmax, Rmax, Rmin, Rmin), (cmin, cmin, cmax, cmax, cmin))
#        axRc.plot((Rs1, Rs2, Rs3, Rs4, Rs1), (cmin, cmax, cmax, cmin, cmin), label='sigma, R constraints')
#        axRc.plot((Rmin, Rmax, Rmax, Rmin, Rmin), (cs1, cs2, cs3, cs4, cs1), label='sigma, c constraints')        
#        axRc.plot((RA1, RA2, RA3, RA4, RA1), (cmin, cmax, cmax, cmin, cmin), label='A, R constraints')
#        axRc.plot((Rmin, Rmax, Rmax, Rmin, Rmin), (cA1, cA2, cA3, cA4, cA1), label='A, c constraints')        
#        
#        # find borders of the (current) design space with a (series of) np.where calls -> where above/below lines
#        # use these borders as inputs into the blade loading limit twice, once for each graph (x no of limits)
#        plt.legend(loc='lower right')
#        plt.show()
        
        return (((Rv1, Rv2, Rv3, Rv4), (Omegav1, Omegav2, Omegav3, Omegav4)),
                ((Rm1, Rm2, Rm3, Rm4), (Omegam1, Omegam2, Omegam3, Omegam4)),
                ((Rs1, Rs2, Rs3, Rs4), (cs1, cs2, cs3, cs4)),
                ((RA1, RA2, RA3, RA4), (cA1, cA2, cA3, cA4)),
                )
        
    def get_available_space(self, coordsX, coordsY, coordsZ, xsRO, xsRc):
        # rip python usage
        # some of the lines are still backwards, (currently) they do not constrain the design space at all
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
                
                if (R > xsRv1[col] and R < xsRv2[col] and R > xsRm1[col] and R < xsRm2[col]):
                    if O > xsOv1[row] and O < xsOv2[row] and O > xsOm1[row] and O < xsOm2[row]:
                        feasible_cells_RO.append([row, col])
                        
                if (R > xsRs1[col] and R < xsRs2[col] and R > xsRA1[col] and R < xsRA2[col]):
                    if c > xscs1[row] and c < xscs2[row] and c > xscA1[row] and c < xscA2[row]:
                        feasible_cells_Rc.append([row, col])

        
        return feasible_cells_RO, feasible_cells_Rc
        
    def plot_boundaries(self, Rvs, Rms, Ovs, Oms, Rss, RAs, css, cAs):
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
         
        xsRO, xsRc = HAMRspace.determine_borders(Rvs, Rms, Ovs, Oms, Rss, RAs, css, cAs)
        xsOv1, xsOv2, xsRv1, xsRv2, xsOm1, xsOm2, xsRm1, xsRm2 = xsRO
        xsRs1, xsRs2, xsRA1, xsRA2, xscs1, xscs2, xscA1, xscA2 = xsRc
        
        dspaceRO, dspaceRc = HAMRspace.get_available_space(coordsX, coordsY, coordsZ, xsRO, xsRc)
        dspaceRO = np.array(dspaceRO)
        dspaceRc = np.array(dspaceRc)
        Ropt = coordsX[dspaceRO[:,0]]
        Oopt = coordsY[dspaceRO[:,1]]
        R2opt = coordsX[dspaceRc[:,0]]
        copt = coordsZ[dspaceRc[:,1]]
        
        fig = plt.figure()
#        axRO = fig.add_subplot(111, projection='3d')
        axRO = fig.add_subplot(111)
        plt.title("Design Space for the HAMRAC R-Omega parameters")
        plt.grid(which='major')
        plt.xlabel('R [m]')
        plt.ylabel('Omega [rad/s]')
#        plt.xlim(3, 8.5)
#        plt.ylim(22, 50)
        
        axRO.plot((Rmin, Rmax, Rmax, Rmin, Rmin), (Omegamin, Omegamin, Omegamax, Omegamax, Omegamin), label='R, Omega constraints')
        axRO.plot((Rv1, Rv2, Rv3, Rv4, Rv1), (Omegamin, Omegamin, Omegamax, Omegamax, Omegamin), label='vtip, R constraints')  
        axRO.plot((Rmin, Rmin, Rmax, Rmax, Rmin), (Omegav1, Omegav2, Omegav3, Omegav4, Omegav1), label='vtip, Omega constraints')        
        axRO.plot((Rm1, Rm2, Rm3, Rm4, Rm1), (Omegamin, Omegamin, Omegamax, Omegamax, Omegamin), label='myu, R constraints')  
        axRO.plot((Rmin, Rmin, Rmax, Rmax, Rmin), (Omegam1, Omegam2, Omegam3, Omegam4, Omegam1), label='myu, Omega constraints')        
        
        plt.scatter(Ropt, Oopt, s=2)

        plt.legend()
        plt.savefig('desspaceRO.png')
        plt.show()
#        
        fig = plt.figure()
        axRc = fig.add_subplot(111)#, projection='3d')
        plt.title("Design Space for the HAMRAC")
        plt.grid(which='major')
        plt.xlabel('R [m]')
        plt.ylabel('c [m]')
#        plt.xlim(3, 10)
#        plt.ylim(0.1, 0.5)

        axRc.plot((Rmin, Rmax, Rmax, Rmin, Rmin), (cmin, cmin, cmax, cmax, cmin))
        axRc.plot((Rs1, Rs2, Rs3, Rs4, Rs1), (cmin, cmax, cmax, cmin, cmin), label='sigma, R constraints')
        axRc.plot((Rmin, Rmax, Rmax, Rmin, Rmin), (cs1, cs2, cs3, cs4, cs1), label='sigma, c constraints')        
        axRc.plot((RA1, RA2, RA3, RA4, RA1), (cmin, cmax, cmax, cmin, cmin), label='A, R constraints')
        axRc.plot((Rmin, Rmax, Rmax, Rmin, Rmin), (cA1, cA2, cA3, cA4, cA1), label='A, c constraints')        
        
        plt.scatter(R2opt, copt, s=2)        
        plt.legend(loc='lower right')
        plt.savefig('desspaceRc.png')
        plt.show()
        
                    
HAMRspace = DesignSpace()
data = HAMRspace.compare_statistics()
HAMRspace.plot_boundaries(data[0][0], data[1][0], data[0][1], data[1][1], data[2][0], data[3][0], data[2][1], data[3][1])
#R1, R2 = HAMRspace.statistical_boundaries['R']
#Omega1, Omega2 = HAMRspace.statistical_boundaries['Omega']
#c1, c2 = HAMRspace.statistical_boundaries['c']
#coordsX = np.linspace(R1, R2)
#coordsY = np.linspace(Omega1, Omega2)
#coordsZ = np.linspace(c1, c2)
#xsRO, xsRc = HAMRspace.determine_borders(data[0][0], data[1][0], data[0][1], data[1][1], data[2][0], data[3][0], data[2][1], data[3][1])
#xsOv1, xsOv2, xsRv1, xsRv2, xsOm1, xsOm2, xsRm1, xsRm2 = xsRO
#xsRs1, xsRs2, xsRA1, xsRA2, xscs1, xscs2, xscA1, xscA2 = xsRc
##Ovb1, Ovb2, Rvb1, Rvb2, Omb1, Omb2, Rmb1, Rmb2 = RVB
##plt.plot(RVB[2], Os, RVB[3], Os, Rs, RVB[0], Rs, RVB[1])
#
#dspaceRO, dspaceRc = HAMRspace.get_available_space(coordsX, coordsY, coordsZ, xsRO, xsRc)
#dspaceRO = np.array(dspaceRO)
#dspaceRc = np.array(dspaceRc)
#Ropt = coordsX[dspaceRO[:,0]]
#Oopt = coordsY[dspaceRO[:,1]]
#
#plt.scatter(Ropt, Oopt)
#
#plt.figure()
#R2opt = coordsX[dspaceRc[:,0]]
#copt = coordsZ[dspaceRc[:,1]]
#
#plt.scatter(R2opt, copt)

