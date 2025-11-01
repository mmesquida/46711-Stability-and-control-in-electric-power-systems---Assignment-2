import numpy as np

from load_data import load_data
from print_results import print_eigenvalues
from print_results import print_P_matrix


if __name__ == "__main__":
    
    # choose the output to be shown
    output = 1
    
    match output:
        case 1:
            # - - - /// Assignment 1.1 /// - - - - - - - - - - - - - - - - - - 
            print('- - - /// Assignment 1.1 /// - - -')
            
            # load data
            data = load_data('q1a')
            A = data['A_q1a']
    
            # computation of right eigenvalues and right eigenvector matrix
            lambda_A, Phi = np.linalg.eig(A)
            
            # computation of damping and frequency
            sigma = np.real(lambda_A)
            omega_d = np.imag(lambda_A) # damped angular frequency
            f_d = np.abs(omega_d) / (2 * np.pi) # damped frequency
            zeta = -sigma / np.sqrt(sigma**2 + omega_d**2) # damping ratio
            
            # print eigenvalues in table format
            print_eigenvalues(lambda_A, f_d, zeta)
            
            # computation of left eigenvector matrix
            Psi = np.linalg.inv(Phi)
            
            # computation of participation matrix
            P = Phi * np.transpose(Psi)
            
            # print participation matrix in table format
            row_headers = ['Δδ', 'Δω', 'Δe_q\'', 'Δe_d\'', 'Δe_q\"', 'Δe_d\"']
            print_P_matrix(P, row_headers)
           
            
        case 2:
            # - - - /// Assignment 1.2 /// - - - - - - - - - - - - - - - - - -
            print('- - - /// Assignment 1.2 /// - - -')
            
            # load data
            data = load_data('q1b')
            A = data['A_q1b']
    
            # computation of right eigenvalues and right eigenvector matrix
            lambda_A, Phi = np.linalg.eig(A)
            
            # computation of damping and frequency
            sigma = np.real(lambda_A)
            omega_d = np.imag(lambda_A) # damped angular frequency
            f_d = np.abs(omega_d) / (2 * np.pi) # damped frequency
            zeta = -sigma / np.sqrt(sigma**2 + omega_d**2) # damping ratio
            
            # print eigenvalues in table format
            print_eigenvalues(lambda_A, f_d, zeta)
            
            # computation of left eigenvector matrix
            Psi = np.linalg.inv(Phi)
            
            # computation of participation matrix
            P = Phi * np.transpose(Psi)
            
            # print participation matrix in table format
            row_headers = ['Δδ', 'Δω', 'Δe_q\'', 'Δe_d\'', 'Δe_q\"', 'Δe_d\"', 'Δv_m_exc', 'Δv_r3_exc', 'Δv_f_exc']
            print_P_matrix(P, row_headers)
            
        
        case 3:
            # - - - /// Assignment 1.3 /// - - - - - - - - - - - - - - - - - -
            print('- - - /// Assignment 1.3 /// - - -')
            
            # load data
            data = load_data('q1c')
            A = data['A_q1c']
    
            # computation of right eigenvalues and right eigenvector matrix
            lambda_A, Phi = np.linalg.eig(A)
            
            # computation of damping and frequency
            sigma = np.real(lambda_A)
            omega_d = np.imag(lambda_A) # damped angular frequency
            f_d = np.abs(omega_d) / (2 * np.pi) # damped frequency
            zeta = -sigma / np.sqrt(sigma**2 + omega_d**2) # damping ratio
            
            # print eigenvalues in table format
            print_eigenvalues(lambda_A, f_d, zeta)
            
            # computation of left eigenvector matrix
            Psi = np.linalg.inv(Phi)
            
            # computation of participation matrix
            P = Phi * np.transpose(Psi)
            
            # print participation matrix in table format
            row_headers = ['Δδ', 'Δω', 'Δe_q\'', 'Δe_d\'', 'Δe_q\"', 'Δe_d\"',
                           'Δv_m_exc', 'Δv_r3_exc', 'Δv_f_exc', 
                           'Δv_1_pss', 'Δv_2_pss', 'Δv_3_pss', 'Δv_ss_pss',]
            print_P_matrix(P, row_headers)
            
            

