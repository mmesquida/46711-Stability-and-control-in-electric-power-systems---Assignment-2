import numpy as np
import matplotlib.pyplot as plt

from load_data import load_data
from compute_results import get_eigenvalues
from compute_results import get_P_matrix
from print_results import print_eigenvalues
from print_results import print_P_matrix
import cmath
from scipy.io import loadmat
from scipy.linalg import eig as la_eig
import re
from Assignment_helpfunctions_part_I.plot_phasor_diagram_sss import plot_phasors

if __name__ == "__main__":

    output = 2
    
    match output:

        case 2:
    
            # - - - /// Assignment 2.1.1 /// - - - - - - - - - - - - - - - - - - - - - 
            print('- - - /// Assignment 2.1.1 /// - - - - - - - - - -')
            
            # load data
            data_q2 = load_data('q2')
            A_q2 = data_q2['A']
            
            # get right eigenvalues, damped frequency and damping ratio
            lambda_A_q2, f_d_q2, zeta_q2 = get_eigenvalues(A_q2)
            
            # print right eigenvalues, damped frequency and damping ratio
            print_eigenvalues(lambda_A_q2, f_d_q2, zeta_q2)
            
            # get participation matrix, right and left eigenvector matrix
            P_q2, Phi_q2, Psi_q2 = get_P_matrix(A_q2)
            
            # print participation matrix
            row_headers_q2 = ['ΔδG1', 'ΔωG1', 'Δψ/f,G1', 'Δψ/kd,G1', 'Δψ/kq1,G1', 
                              'Δψ/kq2,G1','Δv/m,exc,G1','ΔδG2', 'ΔωG2', 'Δψ/f,G2', 
                              'Δψ/kd,G2', 'Δψ/kq1,G2', 'Δψ/kq2,G2','Δv/m,exc,G2',
                              'ΔδG3', 'ΔωG3', 'Δψ/f,G3', 'Δψ/kd,G3', 'Δψ/kq1,G3', 
                              'Δψ/kq2,G3','Δv/m,exc,G3','ΔδG4', 'ΔωG4', 'Δψ/f,G4', 
                              'Δψ/kd,G4', 'Δψ/kq1,G4', 'Δψ/kq2,G4','Δv/m,exc,G4']
            print_P_matrix(P_q2, row_headers_q2)
            
      