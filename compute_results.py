import numpy as np

def get_eigenvalues(A):
    
    # computation of right eigenvalues and right eigenvector matrix
    lambda_A, Phi = np.linalg.eig(A)
    
    # computation of damping and frequency
    sigma = np.real(lambda_A)
    omega_d = np.imag(lambda_A) # damped angular frequency
    f_d = np.abs(omega_d) / (2 * np.pi) # damped frequency
    zeta = -sigma / np.sqrt(sigma**2 + omega_d**2) # damping ratio
    
    return lambda_A, f_d, zeta

def get_P_matrix(A):
    
    # computation of right eigenvalues and right eigenvector matrix
    lambda_A, Phi = np.linalg.eig(A)
    
    # computation of left eigenvector matrix
    Psi = np.linalg.inv(Phi)
    
    # computation of participation matrix
    P = Phi * np.transpose(Psi)
    
    return P, Phi, Psi

