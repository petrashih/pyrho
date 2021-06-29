# -*- coding: utf-8 -*-
"""
Evaluate equilibrium time correlation function
"""
import numpy as np
from scipy import integrate
from pyrho.lib import const, utils

def switch_to_zero(x,x0,x1):
    '''Switching function to go from 1 at x0 to 0 at x1'''
    # Cubic spline with zero derivative on both ends
    y = (x-x0)/(x1-x0)
    if x>x0:
        switch = 1 - 3*y**2 + 2*y**3
    else:
        switch = 1
    return switch  

class EquilibriumTCF(object):
    
    def __init__(self, position, dynamics, picture='None'):
        # dynamics instance must have a dynamics.propagate() method
        # position operator in DVR basis
        self.position = position
        # must be HEOM
        self.dynamics = dynamics
        #
        self.picture = picture
        # not sure whether we need hbar
        # const.hbar = 5308.8
        const.hbar = 1.0
    
    def EquilibrateDynamics(self, rho_0, t_equil, dt):
        # initialize ADOs
        rho_hierarchy0 = self.dynamics.initialize_from_rdm(rho_0)
        # propagate rho (ADOs)
        time1s, rho_hierarchys = self.dynamics.propagate_full(rho_hierarchy0, 0.0, t_equil, dt, self.picture)
        # the last set of ADOs, time, all RDM
        rhos_rdm = self.dynamics.reduce_to_rdm(rho_hierarchys)
        # take care of picture transform
        if self.picture == 'schroedinger':
        	rhos_site = rhos_rdm
        	return rho_hierarchys[-1], time1s, rhos_site
        elif self.picture == 'interaction':
        	rhos_int = rhos_rdm
        	rhos_site = self.dynamics.rhos_int2site(self.dynamics.ham.ham_sys, rhos_int, time1s)
        	return rho_hierarchys[-1], time1s, rhos_site, rhos_int        
    
    def EvaluateTCF(self, rho_hierarchy_init, T_init, T_final, dT, e_min, e_max, de):
        # operate on position operator
        qrho_hierarchy_init = self.dynamics.act_from_left(self.position, rho_hierarchy_init)
        # propagate qrho (ADOs)
        time2s, qrho_hierarchys = self.dynamics.propagate_full(qrho_hierarchy_init, T_init, T_final, dT, 'schroedinger')
        # take RDM at each time                                         ## needs to be modified
        qrho_RDMs = self.dynamics.reduce_to_rdm(qrho_hierarchys)
        # take care of picture transform
#        if self.picture == 'schroedinger':
#        	qrhos_site = qrho_RDMs
#        elif self.picture == 'interaction':
#        	qrhos_int = qrho_RDMs
#        	qrhos_site = self.dynamics.rhos_int2site(self.dynamics.ham.ham_sys, qrhos_int, time2s)
        # calculate C(t) and C(w) by fourier transform of C_damp(t)
        Ct = np.zeros(len(time2s), dtype=complex)
        # Ct_damp = np.zeros(len(time2s), dtype=complex)
        # T_damp = 0.
        
        energies = np.arange(e_min, e_max, de)
        omegas = energies/const.hbar
        Cw = np.zeros(len(omegas))
        # Cw_damp = np.zeros(len(omegas))
        expi = np.exp(1j*(np.outer(omegas,time2s)))        
        for T, Time in enumerate(time2s):
            weight = dT - 0.5*dT*(T==0 or T==len(time2s)-1)
            Ct[T] = np.trace(np.dot(self.position, qrho_RDMs[T]))
            Cw += weight*(expi[:,T]*Ct[T]).real
            # if T > T_damp:
            #     Ct_damp[T] = Ct[T] * switch_to_zero(Time, T_damp, T_final)
            # else:
            #     Ct_damp[T] = Ct[T]
            # Cw_damp += weight*(expi[:,T]*Ct_damp[T]).real

        # # calculate D_imp by (Eq. 6)
        # Dt = 1/1j * (Ct - np.conj(Ct))
        # Dw = np.zeros(len(omegas), dtype=complex)
        # for T, Time in enumerate(time2s):
        #     weight = dT - 0.5*dT*(T==0 or T==len(time2s)-1)
        #     Dw += weight*(expi[:,T]*Dt[T])
 
        # return time2s, Ct, Ct_damp, omegas, Cw, Cw_damp, Dw
        return time2s, Ct, omegas, Cw