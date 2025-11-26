import scipy.io as sio
import numpy as np
import control.statesp as stsp
import stsp_functions as f_stsp
import os
from compute_results import get_eigenvalues, get_state_names, biorthonormalize
from print_results import print_eigenvalues
import matplotlib.pyplot as plt
from scipy.linalg import eig as la_eig
from pathlib import Path
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

pss_2 = f_stsp.pss_stsp(Ks, Tw, Tn1, Td1, Tn2, Td2)
gen_2 = 3
pss2_sys = f_stsp.addPSS(pss_sys, pss_2, strsps, gen_2)

lambda_A_pss2_sys, f_d_pss2_sys, zeta_pss2_sys = get_eigenvalues(pss2_sys.A)
print(f'\nORIGINAL SYSTEM + 2 PSS at G{gen_1} and G{gen_2} (developed PSS settings)')
print_eigenvalues(lambda_A_pss2_sys, f_d_pss2_sys, zeta_pss2_sys)

#Plot time serie for inter-area mode

n_states = pss2_sys.A.shape[0]
state_names = get_state_names(pss2_sys.A, n_states)
lam, VL, VR = la_eig(pss2_sys.A, left=True, right=True)
VL, VR = biorthonormalize(VL, VR)

for k in range(VR.shape[1]):
    den = VL[:, k].conj().T @ VR[:, k]
    if den != 0:
        VL[:, k] /= den.conj()


G = len(state_names)//7
delta_idx = [0 + 7*g for g in range(G)]
           
k_inter = 22

# partner + choose +imag member
def conj_partner(idx): return int(np.argmin(np.abs(lambda_A_pss2_sys - np.conj(lambda_A_pss2_sys[idx]))))
k_p = k_inter if np.real(k_inter) > 0 else conj_partner(k_inter)
k_m = conj_partner(k_p)

print(f"Inter-area: k+={k_p+1}, f={f_d_pss2_sys[k_p]:.3f} Hz, ζ={zeta_pss2_sys[k_p]:.3f}; partner k-={k_m+1}")
    
z0 = np.zeros_like(lambda_A_pss2_sys, dtype=complex)
z0[k_p] = 0.5; z0[k_m] = 0.5

t = np.arange(0.0, 10.0+0.01, 0.01)
exp_terms = np.exp(np.outer(lambda_A_pss2_sys, t))          
x_t = VR @ (exp_terms * z0[:, None])
x_t = np.real(x_t)

# --- select only Δδ states for generador
# --- number of generators and indices for Δδ states ---
G = len(state_names) // 7
delta_idx = [0 + 7*g for g in range(G)]

# --- reconstruct state trajectory for Δδ states ---
y = x_t[delta_idx, :]
m = np.max(np.abs(y));  
y = y/m if m > 0 else y

# --- labels and colors per generator ---
labels = [f"G{g+1}" for g in range(G)]
palette = plt.get_cmap('tab10').colors
COLORS = {lab: palette[i % len(palette)] for i, lab in enumerate(labels)}

plt.figure(figsize=(9, 4))
for g, lab in enumerate(labels):
    plt.plot(t, y[g, :], label=lab, color=COLORS[lab])

plt.xlabel("Time [s]")
plt.ylabel("Normalized $\\Delta\\delta$ [–]")
plt.title(f"Inter-area mode only: f={f_d_pss2_sys[k_p]:.3f} Hz, ζ={zeta_pss2_sys[k_p]:.3f}")
plt.grid(True, alpha=0.3)
plt.legend(ncols=4, loc="upper right")
plt.tight_layout()
plt.show()