import numpy as np

from load_data import load_data


if __name__ == "__main__":
    
    # choose the output to be shown
    output = 1
    
    match output:
        case 1:
            # - - - /// Assignment 1.1.1 /// - - - - - - - - - - - - - - - - -
            print('- - - /// Assignment 1.1.1 /// - - -')
            
            # load data
            data = load_data('q1a')
            A = data['A_q1a']
    
            # computation of right eigenvalues and right eigenvector matrix
            lambda_A, Phi = np.linalg.eig(A)
            
            # computation of damping and frequency of oscillatory modes
            sigma = np.real(lambda_A)
            omega_d = np.imag(lambda_A) # damped angular frequency
            f_d = np.abs(omega_d) / (2 * np.pi) # damped frequency
            
            zeta = -sigma / np.abs(lambda_A) # damping ratio
            
            omega_n = omega_d / np.sqrt(1 - zeta**2) # undamped natural angular frequency
            f_n = np.abs(omega_n) / (2 * np.pi) # undamped natural frequency
            
            # print right eigenvalues
            print('eigenvalues: ')
            for idx, l in enumerate(lambda_A):
                print(f'λ{idx+1} = {l:.4f}')
                
            # print damped frequency of oscillatory modes
            print('damped frequency of oscillatory modes: ')
            for idx, f in enumerate(f_d):
                if omega_d[idx] != 0:
                    print(f'f_d{idx+1} = {f:.4f}')
                    
            # print damping ratios of oscillatory modes
            print('damping ratio of oscillatory modes: ')
            for idx, z in enumerate(zeta):
                if omega_d[idx] != 0:
                    print(f'ζ{idx+1} = {z:.4f}')
            
            # print undamped natural frequency of oscillatory modes
            print('damped frequency of oscillatory modes: ')
            for idx, f in enumerate(f_n):
                if omega_d[idx] != 0:
                    print(f'f_n{idx+1} = {f:.4f}')
            

