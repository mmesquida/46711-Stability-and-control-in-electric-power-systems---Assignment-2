import numpy as np
import matplotlib.pyplot as plt

from load_data import load_data
from compute_results import get_eigenvalues
from compute_results import get_P_matrix, get_state_names, pair_modes, gen_of, top_generators, describe_modes, biorthonormalize, all_gen_ids, cat_of
from print_results import print_eigenvalues
from print_results import print_P_matrix
from scipy.io import loadmat
from scipy.linalg import eig as la_eig
from Assignment_helpfunctions_part_I.plot_phasor_diagram_sss import plot_phasors
from print_results import print_P_matrix, build_modes_table_latex, write_tex_file 
from Assignment_helpfunctions_part_I.P_matrix_write import latex_P_matrix


if __name__ == "__main__":

    output = 2_2
    
    match output:

        case 2_2:
    
            # - - - /// Assignment 2.2 /// - - - - - - - - - - - - - - - - - - - - - 
            
            print('- - - /// Assignment 2.2 /// - - - - - - - - - -')

            # load data
            data_q2 = load_data('q2')

            A_q2 = data_q2['A']

           # get right eigenvalues, damped frequency and damping ratio
            lambda_A_q2, f_d_q2, zeta_q2 = get_eigenvalues(A_q2)

            # get participation matrix, right and left eigenvector matrix
            P_q2, Phi_q2, Psi_q2 = get_P_matrix(A_q2)  