# -*- coding: utf-8 -*-

import scipy.io as sio
import numpy as np
import control.statesp as stsp
import matplotlib.pyplot as plt
import stsp_functions as f_stsp
import os


""" Parameter definitions """
# Parameters for frequency sweep
fmin = 0 ; fmax = 2
f = np.arange(fmin,fmax,.001)
wsys = 2 * np.pi * f

# PSS design parameters
Ks=20;
Tw=16;
Tn1=0.0847;
Td1=0.0000000001;
Tn2=3.474;
Td2=5.45;


""" Load initial System Data """
current_path = os.getcwd()
parent_path = os.path.dirname(current_path)
additional_path = f'/Assignment_data/system_q2.mat'
file_path = parent_path + additional_path

sys_data = sio.loadmat(file_path, squeeze_me=True) #squeeze_me=True gets rid of unnecessary nesting
A = sys_data['A']
B = sys_data['B']
C = sys_data['C']
D = sys_data['D']
strsps = sys_data['strsps']
statename = sys_data['StateName']

""" Form State space object of the original system """
sys = stsp.StateSpace(A,B,C,D)


"""###################### Your code starts here ############################"""


""" Create PSS state space object """
"""############################################################################

the function pss_stsp in stsp_functions.py can be used to construct a 
PSS state space object based on the provided settings

############################################################################"""


""" Find the location of relevant indices in the original system """
"""############################################################################

statename contains the location of the states:
statename == 'psif'
    returns a vector of size (no. of states) which is "True" for all indices related to the field flux linkage states
    
np.logical_or() can be used to return True for multiple states
np.logical_not() would return all states but the indicated ones (convenient for excluding states)

strsps contains the indices related to the in- and outputs of the system:
    strsps['pe'].item() : returns the indices associated with the electric power output
All indices occur in generator order, i.e the first index is related to machine 1, etc.

############################################################################"""

plt.figure()

# manipulate the system by removing the dynamic states
is_r = np.array([(name == 'angle') or (name == 'speed') for name in statename]) # state variables associated with rotor dynamics
is_e = ~is_r # state variables associated with electrical and control dynamics

a_ee = A[is_e, :][:, is_e] # extract submatrix

# frequency response of the system for all generators
for gen in range(0, 4): # for every generator
    vref_idx = strsps['vref'].item()[gen] - 1
    pe_idx   = strsps['pe'].item()[gen]   - 1
    b_e = B[is_e, vref_idx] # extract submatrix
    c_e  = C[pe_idx, is_e] # extract submatrix
    d_e  = D[pe_idx, vref_idx] # extract submatrix

    sys_red = stsp.StateSpace(a_ee, b_e, c_e, d_e) # reduced state space object
    
    mag_sys, phase_sys, omega_sys = sys_red.frequency_response(wsys, squeeze=None) # frequency response of the system
    
    phase_sys_deg = np.rad2deg(np.squeeze(phase_sys))
    phase_sys_deg = (phase_sys_deg + 180) % 360 - 180   
    
    ideal_lead_deg = -phase_sys_deg                     
    
    plt.plot(f, ideal_lead_deg, label=f'G{gen+1}')

# frequency response of the PSS
pss = f_stsp.pss_stsp(Ks, Tw, Tn1, Td1, Tn2, Td2)

mag_pss, phase_pss, omega_pss = pss.frequency_response(wsys, squeeze=None)
phase_pss_deg = np.rad2deg(np.squeeze(phase_pss))
phase_pss_deg = (phase_pss_deg + 180) % 360 - 180
plt.plot(f, phase_pss_deg, label = 'PSS')

plt.xlabel('Frequency [Hz]')
plt.ylabel('Phase [°]')
plt.grid(True)
plt.xlim(0, 2) # adjustment of the x-limits
plt.ylim(-10, 60) # adjustment of the y-limits
plt.legend()
plt.title('Ideal PSS Phase Lead for Generators G1-G4')
plt.tight_layout()
plt.show()
    
    






""" Determine the frequency response of the system and the PSS"""
"""############################################################################

mag, sys_phase, omega = sys_red.frequency_response(wsys, squeeze=None)
returns magnitude and phase of the state space system (single input single output) for the frequencies specified in wsys

############################################################################"""



"""Plot the frequency response curves"""
"""############################################################################

It might be a good idea to invert the system phase.
The ideal pss phase should fully compensate the phase of the system (pss_phase-sys_phase=0). 
Additionally you could plot the phase mismatch between the ideal compensation and the compensation using your parameters.

############################################################################"""




