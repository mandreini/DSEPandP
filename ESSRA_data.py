5# -*- coding: utf-8 -*-
"""
Created on Wed Dec 12 14:37:18 2018

@author: mandreini
"""

import numpy as np
import heli
import copy

MTOW = 2779
R = 4.46
c = 0.227
Omega = 41

#R = 6.2
#c = 0.27
#Omega = 29

heliinputs = 'HAMRAC', 2, 3, 1.05*MTOW*9.81, 8950, R, c, 3, 650, 0.95, Omega, 2, 0.012, 1.15
# 5% added to weight due to the downdraft on the fuselage, must be changed if the actual weight of the rotorcraft becomes relevant
airspeed = np.arange(1, 120, 1, dtype=float)

#far = heli.Coaxheli('FAR-regulation-check', 1, 6, MTOW*9.81, 6400, 0.457, 3, 650, 0.95, 21.21, 2, 0.012, 1.15, convert=True)
#SR_HAMRAC = heli.SRheli(*heliinputs)
#marilena = heli.SRheli("Marilena's Example", 1, 4, 90e3, 0, 9.5, 0.457, 3, 2500, 0.95, 21, 3, 0.01, 1.15)

#mp = marilena.powerCurve(airspeed, select_engine=False)
#mpcl = marilena.powerCurve(airspeed, P_to0=mp[1], select_engine=False)

# Single Rotor Sample Power curve
#level_speed, level_power = SR_HAMRAC.powerCurve(airspeed, P_to0=None, figtitle='Power requirements for a rotorcraft in level, horizontal, forward flight', fname='powercurves/PowerCurveHAMRACLevel.png')
#level_vchars = SR_HAMRAC.determineV_chars(airspeed, level_power/1000)
#climb_speed, climb_power = SR_HAMRAC.powerCurve(airspeed, level_power, 'Power requirements for a rotorcraft in non-level, climbing flight', fname='powercurves/PowerCurveHAMRACClimb.png', select_engine=False)
#climb_vchars = SR_HAMRAC.determineV_chars(airspeed, level_power/1000)

# Coaxial Rotor Power Curve
coax_HAMRAC = heli.Coaxheli(*heliinputs, zD=0.05, convert=True)
clvl_v, clvl_p = coax_HAMRAC.powerCurve(airspeed, P_to0=None, figtitle='Power requirements for a coaxial rotorcraft in \nlevel, horizontal, forward flight', fname='powercurves/PowerCurveCoaxLevel.png', select_engine=True)
#ccl_v, ccl_p = coax_HAMRAC.powerCurve(airspeed, clvl_p, 'Power requirements for a coaxial rotorcraft in \nnon-level, climbing flight', fname='powercurves/PowerCurveCoaxClimb.png', select_engine=False)

# Rotor size difference between single-rotor and coaxial-rotor configurations
#print('Parameter | single-rotor | coaxial-rotor\nNr |', SR_HAMRAC.Nr, '|', Coax_HAMRAC.Nr, '\nNb |', SR_HAMRAC.Nb, '|', SR_HAMRAC.Nb, '\nR |', SR_HAMRAC.R, '|', Coax_HAMRAC.R, '\nc |', SR_HAMRAC.c, '|', Coax_HAMRAC.c, '\nOmega |', SR_HAMRAC.Omega, '|', Coax_HAMRAC.Omega)

# Parameter optimisation
#Rvs, Rdata = coax_HAMRAC.idealPhovafoR()
#coax_HAMRAC.plot_val(Rdata[0], Rdata[1][0]/1000, title='Optimising power required to hover wrt blade radius', xlabel='Blade radius R [m]', ylabel='P_hov [kW]', get_min=True, fname='mdo/P_hovafoR.png')
##
#cvs, cdata = coax_HAMRAC.idealPhovafoc()
#coax_HAMRAC.plot_val(cdata[0], cdata[1][0]/1000, title='Optimising power required to hover wrt blade chord', xlabel='Blade chord c [m]', ylabel='P_hov [kW]', get_min=True, fname='mdo/P_hovafoc.png')
##
#Omegavs, Omegadata = coax_HAMRAC.idealPhovafoOmega(np.linspace(1, 40))
#coax_HAMRAC.plot_val(Omegadata[0], Omegadata[1][0]/1000, title='Optimising power required to hover wrt angular velocity', xlabel='Angular velocity Omega [rad/s]', ylabel='P_hov [kW]', get_min=True, fname='mdo/P_hovafoOmega.png')

# 'ideal' rotorcraft wrt hover power
i_Phov_coax = copy.copy(coax_HAMRAC)
#i_Phov_coax.setR(7.73)
#i_Phov_coax.setc(0.4265)
#i_Phov_coax.setOmega(29.795)
i_Phov_coax.setR(5.55)
i_Phov_coax.setc(0.304)
i_Phov_coax.setOmega(33)
i_Phov_coax.convert()
#i_Phov_coax.Nb = 3  # copying coax_HAMRAC copies Nb and Nr, but Nr is then divided by 2 in convert
i_Phov_vals = i_Phov_coax.reconvert()

HAMRAC_vals = coax_HAMRAC.reconvert()
HAMRAC_vals = 1, 6, R, c, Omega

#validity = i_Phov_coax.Mtip_ne_constraint(), i_Phov_coax.Mtip_cr_constraint(), i_Phov_coax.myu_constraint(), i_Phov_coax.bl_constraint()
ilvl_v, ilvl_p = i_Phov_coax.powerCurve(airspeed, P_to0=None, figtitle="Power requirements for the 'ideal' coaxial rotorcraft in \nlevel, horizontal, forward flight \n'ideal' with respect to minimum hover power", fname='powercurves/PowerCurveidealLevel.png', select_engine=True)
#icl_v, icl_p = i_Phov_coax.powerCurve(airspeed, P_to0=ilvl_p, figtitle="Power requirements for the 'ideal' coaxial rotorcraft in \nnon-level, climbing, forward flight \n'ideal' with respect to minimum hover power", fname='powercurves/PowerCurveidealClimb.png', select_engine=False)
#
#print('P_hov_HAMRAC: %f' % clvl_p[0])
#print('P_hov_ideal: %f' % ilvl_p[0])
#print('Parameter | HAMRAC | ideal')
##print('Nr        | %i      | %i' % (1, 6))
##print('Nb        | %i      | %i' % (1, 6))
#print('R         | %s | %s' % (str(HAMRAC_vals[2])[:6], str(i_Phov_vals[2])[:5]))
#print('c         | %s | %s' % (str(HAMRAC_vals[3])[:6], str(i_Phov_vals[3])[:5]))
#print('Omega     | %s | %s' % (str(HAMRAC_vals[4])[:6], str(i_Phov_vals[4])[:5]))

# comparison between old and new
#i_Phov_coax.plot_vals(airspeed, (ilvl_p/1000, clvl_p/1000), labels=('ideal curve', 'current curve'), title='Comparison between current and current ideal rotor design', xlabel='airspeed [m/s]', ylabel='P_hov [kW]', get_min=False, fname='PowerCurveComparison.png')


better_HAMRAC = copy.copy(coax_HAMRAC)
# SR rotor
better_HAMRAC.setR(6.5)
better_HAMRAC.setc(0.301)
better_HAMRAC.setOmega(31) 
better_HAMRAC.convert()

blvl_v, blvl_p = better_HAMRAC.powerCurve(airspeed)
better_vals = better_HAMRAC.reconvert()

print('P_hov_HAMRAC: %f' % clvl_p[0])
print('P_hov_ideal: %f' % ilvl_p[0])
print('P_hov_better: %f' % blvl_p[0])
print('Parameter | HAMRAC | ideal | better')
#print('Nr        | %i      | %i' % (1, 6))
#print('Nb        | %i      | %i' % (1, 6))
print('R         | %s | %s | %s' % (str(HAMRAC_vals[2])[:6], str(i_Phov_vals[2])[:5], str(better_vals[2])[:5]))
print('c         | %s | %s | %s' % (str(HAMRAC_vals[3])[:6], str(i_Phov_vals[3])[:5], str(better_vals[3])[:5]))
print('Omega     | %s | %s | %s' % (str(HAMRAC_vals[4])[:6], str(i_Phov_vals[4])[:5], str(better_vals[4])[:5]))

# Characteristic airspeed parameter optimisation (2D parameter optimisation)
# R
#v_charRvs, v_charRdata = coax_HAMRAC.v_charsafoR(airspeed)
#v_charRs, v_charRspeeds = v_charRdata
#v_charRs = v_charRs[0]  # 'unmesh'
#v_ARs, v_BRs, v_CRs, p_minRs, p_sfrRs = v_charRspeeds
#v_ARs, v_BRs, v_CRs = v_ARs[::-1], v_BRs[::-1], v_CRs[::-1]
#
## c
#v_charcvs, v_charcdata = coax_HAMRAC.v_charsafoc(airspeed)
#v_charcs, v_charcspeeds = v_charcdata
#v_charcs = v_charcs[0]  
#v_Acs, v_Bcs, v_Ccs, p_minRs, p_sfrRs = v_charcspeeds
#v_Acs, v_Bcs, v_Ccs = v_Acs[::-1], v_Bcs[::-1], v_Ccs[::-1]
#
## Omega
#v_charOmegavs, v_charOmegadata = coax_HAMRAC.v_charsafoOmega(airspeed)
#v_charOmegas, v_charOmegaspeeds = v_charOmegadata
#v_charOmegas = v_charOmegas[0]  
#v_AOmegas, v_BOmegas, v_COmegas, p_minOmegas, p_sfrOmegas = v_charOmegaspeeds
#v_AOmegas, v_BOmegas, v_COmegas = v_AOmegas[::-1], v_BOmegas[::-1], v_COmegas[::-1]

#coax_HAMRAC.plot_val(v_charRs, v_ARs, xlabel='Radius of rotor R [m]', ylabel='v_pmin; minimum power airspeed [m/s]', title='Comparison between minimum power airspeed and rotor diameter', fname='mdo/vpminR.png')
#coax_HAMRAC.plot_val(v_charRs, v_BRs, xlabel='Radius of rotor R [m]', ylabel='v_sfr; specific flight range airspeed [m/s]', title='Comparison between specific flight range airspeed and rotor diameter', fname='mdo/vsfrR.png')
#
#coax_HAMRAC.plot_val(v_charcs, v_Acs, xlabel='Length of chord c [m]', ylabel='v_pmin; minimum power airspeed [m/s]', title='Comparison between minimum power airspeed and chord length', fname='mdo/vpminc.png')
#coax_HAMRAC.plot_val(v_charcs, v_Bcs, xlabel='Length of chord c [m]', ylabel='v_sfr; specific flight range airspeed [m/s]', title='Comparison between specific flight range airspeed and chord length', fname='mdo/vsfrc.png')
#
#coax_HAMRAC.plot_val(v_charOmegas, v_AOmegas, xlabel='Angular velocity Omega [rad/s]', ylabel='v_pmin; minimum power airspeed [m/s]', title='Comparison between minimum power airspeed and angular velocity', fname='mdo/vpminOmega.png')
#coax_HAMRAC.plot_val(v_charOmegas, v_BOmegas, xlabel='Angular velocity Omega [rad/s]', ylabel='v_sfr; specific flight range airspeed [m/s]', title='Comparison between specific flight range airspeed and angular velocity', fname='mdo/vsfrOmega.png')
# assume that v_C has similar comparison to v_B as they both are largely dependent on parasite drag

# Preliminary data power curve, Leon, 8.2.4
#m1pavel3 = heli.Coaxheli('m1pavel3', 2, 3, 2657*9.81, 8950, 4.46, 0.304, 3, 650, 0.95, 41.8, 2, 0.012, 1.15)
#m2krenik3 = heli.Coaxheli('m2krenik3', 2, 3, 2657*9.81, 8950, 4.46, 0.304, 3, 650, 0.95, 41.8, 2, 0.012, 1.15)
#m3campfens3 = heli.Coaxheli('m3krenik3', 2, 3, 2657*9.81, 8950, 4.46, 0.293, 3, 650, 0.95, 48.202, 2, 0.012, 1.15)
#
#m1pavel4 = heli.Coaxheli('m1pavel4', 2, 4, 2657*9.81, 8950, 4.46, 0.227, 3, 650, 0.95, 41.8, 2, 0.012, 1.15)
#m2krenik4 = heli.Coaxheli('m2krenik3', 2, 4, 2657*9.81, 8950, 4.46, 0.227, 3, 650, 0.95, 41.8, 2, 0.012, 1.15)
#m3campfens4 = heli.Coaxheli('m3krenik3', 2, 4, 2657*9.81, 8950, 4.46, 0.220, 3, 650, 0.95, 48.202, 2, 0.012, 1.15)
##
#m13_lvl_v, m13_lvl_p = m1pavel3.powerCurve(vs=airspeed, P_to0=None, figtitle='Power requirements for %s in level forward flight' % m1pavel3.name, fname='leonvals/pcurvem13_lvl.png')
#m23_lvl_v, m23_lvl_p = m2krenik3.powerCurve(vs=airspeed, P_to0=None, figtitle='Power requirements for %s in level forward flight' % m2krenik3.name, fname='leonvals/pcurvem23_lvl.png')
#m33_lvl_v, m33_lvl_p = m3campfens3.powerCurve(vs=airspeed, P_to0=None, figtitle='Power requirements for %s in level forward flight' % m3campfens3.name, fname='leonvals/pcurvem33_lvl.png')
#m14_lvl_v, m14_lvl_p = m1pavel4.powerCurve(vs=airspeed, P_to0=None, figtitle='Power requirements for %s in level forward flight' % m1pavel4.name, fname='leonvals/pcurvem14_lvl.png')
#m24_lvl_v, m24_lvl_p = m2krenik4.powerCurve(vs=airspeed, P_to0=None, figtitle='Power requirements for %s in level forward flight' % m2krenik4.name, fname='leonvals/pcurvem24_lvl.png')
#m34_lvl_v, m34_lvl_p = m3campfens4.powerCurve(vs=airspeed, P_to0=None, figtitle='Power requirements for %s in level forward flight' % m3campfens4.name, fname='leonvals/pcurvem34_lvl.png')
#
#m13_cl_v, m13_vl_p = m1pavel3.powerCurve(vs=airspeed, P_to0=m13_lvl_p, figtitle='Power requirements for %s in non-level, climbing flight' % m1pavel3.name, fname='leonvals/pcurvem13_cl.png')
#m23_cl_v, m23_cl_p = m2krenik3.powerCurve(vs=airspeed, P_to0=m23_lvl_p, figtitle='Power requirements for %s in non-level, climbing flight' % m2krenik3.name, fname='leonvals/pcurvem23_cl.png')
#m33_cl_v, m33_cl_p = m3campfens3.powerCurve(vs=airspeed, P_to0=m33_lvl_p, figtitle='Power requirements for %s in non-level, climbing flight' % m3campfens3.name, fname='leonvals/pcurvem33_cl.png')
#m14_cl_v, m14_cl_p = m1pavel4.powerCurve(vs=airspeed, P_to0=m14_lvl_p, figtitle='Power requirements for %s in non-level, climbing flight' % m1pavel4.name, fname='leonvals/pcurvem14_cl.png')
#m24_cl_v, m24_cl_p = m2krenik4.powerCurve(vs=airspeed, P_to0=m24_lvl_p, figtitle='Power requirements for %s in non-level, climbing flight' % m2krenik4.name, fname='leonvals/pcurvem24_cl.png')
#m34_cl_v, m34_cl_p = m3campfens4.powerCurve(vs=airspeed, P_to0=m34_lvl_p, figtitle='Power requirements for %s in non-level, climbing flight' % m3campfens4.name, fname='leonvals/pcurvem34_cl.png')

#data = SR_HAMRAC.compare_statistics()

