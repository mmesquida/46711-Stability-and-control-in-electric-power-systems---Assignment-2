
"""
Loading the matlab files into python for part 1
"""

import numpy as np
import scipy.io as sio

data_folder = 'Assignment_data'

Case = 'q1a'


if Case == 'q1a': 
    file = f"{data_folder}/system_{Case}.mat"
    data = sio.loadmat(file)
    
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
