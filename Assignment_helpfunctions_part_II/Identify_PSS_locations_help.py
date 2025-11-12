# -*- coding: utf-8 -*-

import scipy.io as sio
import numpy as np
import os

from compute_results import get_eigenvalues, get_P_matrix, biorthonormalize
from print_results import print_eigenvalues, print_P_matrix
from scipy.linalg import eig as la_eig

""" Load System Data """
current_path = os.getcwd()
parent_path = os.path.dirname(current_path)
additional_path = f'/Assignment_data/system_q2.mat'
file_path = parent_path + additional_path

sys_data = sio.loadmat(file_path ,squeeze_me=True) #squeeze_me=True gets rid of unnecessary nesting
A = sys_data['A']
B = sys_data['B']
C = sys_data['C']
D = sys_data['D']
strsps = sys_data['strsps']
statename = sys_data['StateName']

""" Determine Eigenvalues and calculate Eigenvectors and participation matrix """
# get right eigenvalues, damped frequency and damping ratio
lambda_A, f_d, zeta = get_eigenvalues(A)

# get participation matrix, right and left eigenvector matrix
P, Phi, Psi = get_P_matrix(A)

# print participation matrix
row_headers = ['ΔδG1', 'ΔωG1', 'Δψ/f,G1', 'Δψ/kd,G1', 'Δψ/kq1,G1', 'Δψ/kq2,G1','Δv/m,exc,G1','ΔδG2', 'ΔωG2', 'Δψ/f,G2', 'Δψ/kd,G2', 'Δψ/kq1,G2', 'Δψ/kq2,G2','Δv/m,exc,G2','ΔδG3', 'ΔωG3', 'Δψ/f,G3', 'Δψ/kd,G3', 'Δψ/kq1,G3', 'Δψ/kq2,G3','Δv/m,exc,G3','ΔδG4', 'ΔωG4', 'Δψ/f,G4', 'Δψ/kd,G4', 'Δψ/kq1,G4', 'Δψ/kq2,G4','Δv/m,exc,G4']
#print_P_matrix(P, row_headers)

# getting the right and left eigenvalues
lam, VL, VR = la_eig(A, left=True, right=True)
VL, VR = biorthonormalize(VL, VR) 

""" Find the location of relevant indices in the original system """

""" reduce analysis to electromechanical modes """
"""############################################################################

To ease the analysis the participation matrix can be reduced to reflect only the states of interest in the modes of interest.
To provide a ranging the resulting participation matrix can further be normalized such that the
magitude of the generator with the highest participation is "1".   

############################################################################"""
modes_Interest = [16,17,18,19,20,21,22,23]
idx_states_interest = [0, 1, 7, 8, 14, 15, 21, 22]
P_reduced = np.abs(P[np.ix_(idx_states_interest, modes_Interest)])

""" Calculate transfer function residues  """
"""############################################################################

calculate the observability and controllability
reduce the matrices to observability & controllability of the relevant modes and the appropriate in- and outputs of the system 
 
res=em_w_obs*em_Vr_contr.T 

#should then provide a suitable format where the residues for all indicated modes appear column-wise 
and the rows reflect the different locations 

Normalization can be applied similar to the participation matrix.
############################################################################"""

avr_cont = [0,1,2,3]    
Cw = np.abs(C[idx_states_interest, :] @ VR[:, modes_Interest])
B_Vr = np.abs(VL[:, modes_Interest].T @ B[:, avr_cont]) 
residues = np.abs(Cw @ B_Vr)
    
print(residues)