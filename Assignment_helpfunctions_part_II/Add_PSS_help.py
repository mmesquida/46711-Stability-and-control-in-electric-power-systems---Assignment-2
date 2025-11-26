# -*- coding: utf-8 -*-
import scipy.io as sio
import numpy as np
import control.statesp as stsp
import stsp_functions as f_stsp
import os
from pathlib import Path
import matplotlib.pyplot as plt
from compute_results import get_eigenvalues, get_state_names, cat_of, gen_of, biorthonormalize
from scipy.linalg import eig as la_eig
from print_results import print_eigenvalues

""" Load System Data """
this_dir = Path(__file__).resolve().parent

# Project root: go one level up to .../46711-Stability-and-control-in-electric-power-systems---Assignment-2
project_root = this_dir.parent

# Now point to the MAT file in Assignment_data
file_path = project_root / "Assignment_data" / "system_q2.mat"


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
# -----------------------------------------------------------------------------
# original system
lambda_A_sys, f_d_sys, zeta_sys = get_eigenvalues(sys.A)
print('\nORIGINAL SYSTEM')
print_eigenvalues(lambda_A_sys, f_d_sys, zeta_sys)

# -----------------------------------------------------------------------------
# original system + PSS (provided initial PSS settings)
# added PSS with provided values
Ks=20;
Tw=10;
Tn1=0.05;
Td1=0.02;
Tn2=3.0;
Td2=5.4;
pss = f_stsp.pss_stsp(Ks, Tw, Tn1, Td1, Tn2, Td2)

gen = 1
pss_sys = f_stsp.addPSS(sys, pss, strsps, gen)

lambda_A_pss_sys, f_d_pss_sys, zeta_pss_sys = get_eigenvalues(pss_sys.A)
print(f'\nORIGINAL SYSTEM + PSS at G{gen} (provided initial PSS settings)')
print_eigenvalues(lambda_A_pss_sys, f_d_pss_sys, zeta_pss_sys)

# -----------------------------------------------------------------------------
# original system + PSS (developed PSS settings)
# added PSS with developed values

Ks=20;
Tw=16;
Tn1=0.0847;
Td1=0.0000000001;
Tn2=3.474;
Td2=5.45;

pss = f_stsp.pss_stsp(Ks, Tw, Tn1, Td1, Tn2, Td2)

gen = 1
pss_sys = f_stsp.addPSS(sys, pss, strsps, gen)

lambda_A_pss_sys, f_d_pss_sys, zeta_pss_sys = get_eigenvalues(pss_sys.A)
print(f'\nORIGINAL SYSTEM + PSS at G{gen} (developed PSS settings)')
print_eigenvalues(lambda_A_pss_sys, f_d_pss_sys, zeta_pss_sys)

# -----------------------------------------------------------------------------
# original system + 2 PSS (developed PSS settings)
pss_1 = f_stsp.pss_stsp(Ks, Tw, Tn1, Td1, Tn2, Td2)
gen_1 = 1
pss_sys = f_stsp.addPSS(sys, pss_1, strsps, gen_1)

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





