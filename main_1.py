import numpy as np

from load_data import load_data


if __name__ == "__main__":
    
    # /// Assignment 1.1.1 ///
    data = load_data('q1a')
    A = data['A_q1a']
    
    #compute eigenvalues
    lambda_A, Phi = np.linalg.eig(A)
    
    #computing Frequency and damping ratio
    sigma = np.real(lambda_A)
    omega_d =  np.imag(lambda_A)
    
    zeta = -sigma/np.abs(lambda_A)
    Freq = np.abs(omega_d)/(2* np.pi)
    
    print("Eigenvalues of 1.1:" )
    for lam in lambda_A:
      print(f"{lam.real:+.4f}{'+' if lam.imag >= 0 else '-'}j{abs(lam.imag):.4f}")
    print("\nOscillating modes:")
    for i, lam in enumerate(lambda_A):
      if np.imag(lam) != 0:
        print(f"Mode {i+1}: damping = {zeta[i]:.3f},  f = {Freq[i]:.3f} Hz")

