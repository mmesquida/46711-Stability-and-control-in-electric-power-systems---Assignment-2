# -*- coding: utf-8 -*-
import scipy.io as sio
import numpy as np
import control.statesp as stsp
import stsp_functions as f_stsp
import os
from compute_results import get_eigenvalues
from print_results import print_eigenvalues

""" Load System Data """
current_path = os.getcwd()
parent_path = os.path.dirname(current_path)
additional_path = f'/Assignment_data/system_q2.mat'
file_path = parent_path + additional_path

sys_data = sio.loadmat(file_path, squeeze_me=True) # event. use squeeze_me=True to get rid of unnecessary nesting
A = sys_data['A']
B = sys_data['B']
C = sys_data['C']
D = sys_data['D']
strsps = sys_data['strsps']
statename = sys_data['StateName']

""" Form State space object of the original system """
sys = stsp.StateSpace(A,B,C,D)

""" Add a PSS to the system """
"""############################################################################

Form a PSS state space object "pss" based on the desired settings
Add the pss object to an existing state space system "sys"
The function "pss_sys=addPSS(sys,pss,strsps,gen)" in stsp_functions.py returns a new state space system "pss_sys" 
where the PSS is added to machine "gen".

If another PSS should be added the function can be called again using the system that already contains 
a pss as input. For example: "pss2_sys=addPSS(pss_sys,pss2,strsps,gen2)"

############################################################################"""
# right eigenvalues, damped frequency and damping ratio
# original system
lambda_A_sys, f_d_sys, zeta_sys = get_eigenvalues(sys.A)
print_eigenvalues(lambda_A_sys, f_d_sys, zeta_sys)

# added PSS with provided values
Ks=20;
Tw=10;
Tn1=0.05;
Td1=0.02;
Tn2=3.0;
Td2=5.4;
pss_0 = f_stsp.pss_stsp(Ks, Tw, Tn1, Td1, Tn2, Td2)

gen = 3
pss_sys_0 = f_stsp.addPSS(sys, pss_0, strsps,gen)

lambda_A_pss_sys_0, f_d_pss_sys_0, zeta_pss_sys_0 = get_eigenvalues(pss_sys_0.A)
print_eigenvalues(lambda_A_pss_sys_0, f_d_pss_sys_0, zeta_pss_sys_0)

# added PSS with developed values
Ks=20;
Tw=16;
Tn1=0.0847;
Td1=0.0000000001;
Tn2=3.474;
Td2=5.45;
pss_1 = f_stsp.pss_stsp(Ks, Tw, Tn1, Td1, Tn2, Td2)

gen = 3
pss_sys_1 = f_stsp.addPSS(sys, pss_1, strsps,gen)

lambda_A_pss_sys_1, f_d_pss_sys_1, zeta_pss_sys_1 = get_eigenvalues(pss_sys_1.A)
# print_eigenvalues(lambda_A_pss_sys_1, f_d_pss_sys_1, zeta_pss_sys_1)





""" Compare the systems """
"""############################################################################

Determine the damping and frequency of the different systems

For a start that could be:
    Original System
    System with PSS settings provided
    Improved PSS tuning

############################################################################"""





