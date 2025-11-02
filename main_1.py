import numpy as np
import matplotlib.pyplot as plt

from load_data import load_data
from compute_results import get_eigenvalues
from compute_results import get_P_matrix
from print_results import print_eigenvalues
from print_results import print_P_matrix


if __name__ == "__main__":
    
    # - - - /// Assignment 1.1.1 /// - - - - - - - - - - - - - - - - - - - - - 
    print('- - - /// Assignment 1.1.1 /// - - - - - - - - - -')
    
    # load data
    data_q1a = load_data('q1a')
    A_q1a = data_q1a['A_q1a']
    
    # get right eigenvalues, damped frequency and damping ratio
    lambda_A_q1a, f_d_q1a, zeta_q1a = get_eigenvalues(A_q1a)
    
    # print right eigenvalues, damped frequency and damping ratio
    print_eigenvalues(lambda_A_q1a, f_d_q1a, zeta_q1a)
    
    # get participation matrix, right and left eigenvector matrix
    P_q1a, Phi_q1a, Psi_q1a = get_P_matrix(A_q1a)
    
    # print participation matrix
    row_headers_q1a = ['Δδ', 'Δω', 'Δe_q\'', 'Δe_d\'', 'Δe_q\"', 'Δe_d\"']
    print_P_matrix(P_q1a, row_headers_q1a)
    
    
    # - - - /// Assignment 1.2.1 /// - - - - - - - - - - - - - - - - - - - - - 
    print('- - - /// Assignment 1.2.1 /// - - - - - - - - - -')
    
    # load data
    data_q1b = load_data('q1b')
    A_q1b = data_q1b['A_q1b']
    
    # get right eigenvalues, damped frequency and damping ratio
    lambda_A_q1b, f_d_q1b, zeta_q1b = get_eigenvalues(A_q1b)
    
    # print right eigenvalues, damped frequency and damping ratio
    print_eigenvalues(lambda_A_q1b, f_d_q1b, zeta_q1b)
    
    # get participation matrix, right and left eigenvector matrix
    P_q1b, Phi_q1b, Psi_q1b = get_P_matrix(A_q1b)
    
    # print participation matrix
    row_headers_q1b = ['Δδ', 'Δω', 'Δe_q\'', 'Δe_d\'', 'Δe_q\"', 'Δe_d\"',
                       'Δv_m_exc', 'Δv_r3_exc', 'Δv_f_exc']
    print_P_matrix(P_q1b, row_headers_q1b)
    
    
    # - - - /// Assignment 1.3.1 /// - - - - - - - - - - - - - - - - - - - - - 
    print('- - - /// Assignment 1.3.1 /// - - - - - - - - - -')
    
    # load data
    data_q1c = load_data('q1c')
    A_q1c = data_q1c['A_q1c']
    
    # get right eigenvalues, damped frequency and damping ratio
    lambda_A_q1c, f_d_q1c, zeta_q1c = get_eigenvalues(A_q1c)
    
    # print right eigenvalues, damped frequency and damping ratio
    print_eigenvalues(lambda_A_q1c, f_d_q1c, zeta_q1c)
    
    # get participation matrix, right and left eigenvector matrix
    P_q1c, Phi_q1c, Psi_q1c = get_P_matrix(A_q1c)
    
    # print participation matrix
    row_headers_q1c = ['Δδ', 'Δω', 'Δe_q\'', 'Δe_d\'', 'Δe_q\"', 'Δe_d\"',
                       'Δv_m_exc', 'Δv_r3_exc', 'Δv_f_exc', 
                       'Δv_1_pss', 'Δv_2_pss', 'Δv_3_pss', 'Δv_ss_pss',]
    print_P_matrix(P_q1c, row_headers_q1c)
    
    
    # - - - /// Assignment 1.3.2 /// - - - - - - - - - - - - - - - - - - - - - 
    t = np.arange(0, 5, 0.001)
    
    delta_nul_deg = 5 
    x_nul_q1a = [np.deg2rad(delta_nul_deg), 0, 0, 0, 0, 0]
    x_nul_q1b = [np.deg2rad(delta_nul_deg), 0, 0, 0, 0, 0, 0, 0, 0]
    x_nul_q1c = [np.deg2rad(delta_nul_deg), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    
    x_t_q1a = np.zeros((len(x_nul_q1a), len(t)), dtype=np.complex64)
    x_t_q1b = np.zeros((len(x_nul_q1b), len(t)), dtype=np.complex64)
    x_t_q1c = np.zeros((len(x_nul_q1c), len(t)), dtype=np.complex64)
    
    for k in range(len(t)):
        x_t_q1a[:,k] = Phi_q1a.dot(np.exp(lambda_A_q1a*t[k])*Psi_q1a.dot(x_nul_q1a))
        
    for k in range(len(t)):
        x_t_q1b[:,k] = Phi_q1b.dot(np.exp(lambda_A_q1b*t[k])*Psi_q1b.dot(x_nul_q1b))
        
    for k in range(len(t)):
        x_t_q1c[:,k] = Phi_q1c.dot(np.exp(lambda_A_q1c*t[k])*Psi_q1c.dot(x_nul_q1c))
        
    delta_t_q1a_deg = np.rad2deg(np.real(x_t_q1a[0]))
    delta_t_q1b_deg = np.rad2deg(np.real(x_t_q1b[0]))
    delta_t_q1c_deg = np.rad2deg(np.real(x_t_q1c[0]))
        
    fig = plt.figure(1, figsize=(5, 4), dpi=300) # create the figure
    
    line_q1a, = plt.plot(t, delta_t_q1a_deg, color = 'blue', linewidth = 1)
    line_q1a.set_label('Manual excitation')
    
    line_q1b, = plt.plot(t, delta_t_q1b_deg, color = 'red', linewidth = 1)
    line_q1b.set_label('AVR only')
    
    line_q1c, = plt.plot(t, delta_t_q1c_deg, color = 'green', linewidth = 1)
    line_q1c.set_label('AVR and PSS')
    
    plt.xlim(0, 5) # adjustment of the x-limits
    plt.ylim(-12, 12) # adjustment of the y-limits
    plt.legend(loc = 'upper left', fontsize = 6)
    
    plt.title('Rotor angle response - Small signal stability', fontsize = 10, pad = 4)
    plt.xlabel('time [s]', fontsize = 8, labelpad = 0)
    plt.ylabel('δ - rotor angle [deg]', fontsize = 8, labelpad = 0)
    
    fig.tight_layout() # reduces white space around figures
    plt.margins(x=0, y=0) # results in that (0,0) is in lower left corner
    plt.tick_params(labelsize=8) # reduce the fontsize of the tickmarks
    plt.grid(linestyle='-.', linewidth=0.6) # grid line modifications.
    
    plt.show()
    
    
    
        
    
    
    
    
    
    
    
    
    
    
            
            

