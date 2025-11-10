import numpy as np
import matplotlib.pyplot as plt

from load_data import load_data
from compute_results import get_eigenvalues
from compute_results import get_P_matrix, get_state_names, pair_modes, gen_of, top_generators, describe_modes,normalize_P_matrix, biorthonormalize, all_gen_ids, cat_of
from print_results import print_eigenvalues
from print_results import print_P_matrix
import cmath
from scipy.io import loadmat
from scipy.linalg import eig as la_eig
import re
from Assignment_helpfunctions_part_I.plot_phasor_diagram_sss import plot_phasors
from print_results import print_P_matrix, build_modes_table_latex, write_tex_file 
from Assignment_helpfunctions_part_I.P_matrix_write import latex_P_matrix


if __name__ == "__main__":

    output = 2
    
    match output:

        case 2:
    
            # - - - /// Assignment 2.1.2 /// - - - - - - - - - - - - - - - - - - - - - 
            
            print('- - - /// Assignment 2.1.2 /// - - - - - - - - - -')

            # load data
            data_q2 = load_data('q2')

            A_q2 = data_q2['A']

           # get right eigenvalues, damped frequency and damping ratio
            lambda_A_q2, f_d_q2, zeta_q2 = get_eigenvalues(A_q2)

            # get participation matrix, right and left eigenvector matrix
            P_q2, Phi_q2, Psi_q2 = get_P_matrix(A_q2)           
             # print participation matrix
            row_headers_q2 = ['ΔδG1', 'ΔωG1', 'Δψ/f,G1', 'Δψ/kd,G1', 'Δψ/kq1,G1', 'Δψ/kq2,G1','Δv/m,exc,G1','ΔδG2', 'ΔωG2', 'Δψ/f,G2', 'Δψ/kd,G2', 'Δψ/kq1,G2', 'Δψ/kq2,G2','Δv/m,exc,G2','ΔδG3', 'ΔωG3', 'Δψ/f,G3', 'Δψ/kd,G3', 'Δψ/kq1,G3', 'Δψ/kq2,G3','Δv/m,exc,G3','ΔδG4', 'ΔωG4', 'Δψ/f,G4', 'Δψ/kd,G4', 'Δψ/kq1,G4', 'Δψ/kq2,G4','Δv/m,exc,G4']
            print_P_matrix(P_q2, row_headers_q2)
            latex_P_matrix(P_q2, np.array(row_headers_q2).reshape(-1,1), False, 'P_matrix_q2.tex', 4, 0.1)

            # - - - /// Assignment 2.1.2 /// - - - - - - - - - - - - - - - - - - - - - 
            print('- - - /// Assignment 2.1.2 /// - - - - - - - - - -')
            tol_imag = 1e-8
            top_k = 5
            n_states = A_q2.shape[0]
            state_names = get_state_names(data_q2, n_states)
            lam, VL, VR = la_eig(A_q2, left=True, right=True)
            VL, VR = biorthonormalize(VL, VR)
            groups = pair_modes(lam, tol_imag)
            names = get_state_names(data_q2, A_q2.shape[0])
            descriptions = describe_modes(lam, groups, P_q2, names, top_k=top_k)
            tex = build_modes_table_latex(lam, groups, descriptions,
                                        caption="Eigenvalues, oscillatory frequencies and damping, and dominant states (Q2.1.2).",
                                        label="\\label{tab:q2_modes}") 

            write_tex_file("modes_q2_table.tex", tex)

            # - - - /// Assignment 2.1.4 /// - - - - - - - - - - - - - - - - - - - - - 
            print('- - - /// Assignment 2.1.4 /// - - - - - - - - - -')
            # -- plot phasor diagram for mode k=21 (the first oscillatory mode)
            k_self = 21 - 1
            # If its eigenvalue has negative imaginary part, switch to its conjugate partner so we always plot the member with positive imag
            if np.iscomplex(lam[k_self]) and np.imag(lam[k_self]) < 0:
                cj = np.conj(lam[k_self])
                partner = int(np.argmin(np.abs(lam - cj)))
                if np.imag(lam[partner]) > 0:
                    k = partner

            # --- select only Δδ states for generador
            idx = []
            labels = []
            seen = set()
            for i, nm in enumerate(state_names):
                if cat_of(nm) == r'$\Delta \delta$':     # Used to categorizade states
                    g = gen_of(nm)                       # Get generator number (For each selected name, extract the generator id)     
                    if g and g not in seen:              # avoid duplicates
                        idx.append(i)
                        labels.append(f"G{g}")
                        seen.add(g)

            if not idx:
                raise RuntimeError("Could find Δδ states; review state_names.")

            # ---  build P matrix for selected states and normalize
            v = VR[:, k_self][idx]
            P = np.column_stack([np.real(v), np.imag(v)]) #Create a 2-column array with real and imaginary parts
            m = np.max(np.abs(P)) if P.size else 0.0
            if m > 0:
                P = P / m

            # --- colores simples y título con métricas
            C = ['tab:blue'] * len(labels)
            sig = float(np.real(lam[k_self])); omg = float(np.imag(lam[k_self]))
            f_hz = abs(omg)/(2*np.pi); zeta = -sig/np.hypot(sig, omg) if (sig or omg) else np.nan #Compute frequency and damping for the title

            plt.figure()
            plot_phasors(P, C, labels)
            plt.gca().set_aspect('equal', 'box'); plt.grid(True, alpha=0.3)
            plt.title(f"Mode k={k_self+1} — $\\Delta\\delta$ | f={f_hz:.3f} Hz, ζ={zeta:.3f}")
            plt.show()

            # - - - /// Assignment 2.1.5 /// - - - - - - - - - - - - - - - - - - - - - 

            "- - - - - /// Assignment 2.1.5 ///- - - - - -"

            # --- eigen-decomposition + bi-orthonormalization ---
            for k in range(VR.shape[1]):
                den = VL[:, k].conj().T @ VR[:, k]
                if den != 0:
                    VL[:, k] /= den.conj()

            # --- modal metrics ---
            sigma = np.real(lam); omega = np.imag(lam)
            fd = np.abs(omega)/(2*np.pi)
            zeta = -sigma/np.hypot(sigma, omega)

            # --- rotor-angle indices (Δδ) ---
            G = len(state_names)//7
            delta_idx = [0 + 7*g for g in range(G)]

            # --- pick inter-area: least damped among oscillatory modes ---
            osc = np.where(np.abs(omega) > 1e-8)[0] 
            k_inter = int(osc[np.argmin(zeta[osc])])

            # partner + choose +imag member
            def conj_partner(idx): return int(np.argmin(np.abs(lam - np.conj(lam[idx]))))
            k_p = k_inter if omega[k_inter] > 0 else conj_partner(k_inter)
            k_m = conj_partner(k_p)

            print(f"Inter-area: k+={k_p+1}, f={fd[k_p]:.3f} Hz, ζ={zeta[k_p]:.3f}; partner k-={k_m+1}")

            # --- excite only that pair (real response) ---
            z0 = np.zeros_like(lam, dtype=complex)
            z0[k_p] = 0.5; z0[k_m] = 0.5

            # --- time response ---
            t = np.arange(0.0, 12.0+0.01, 0.01)
            exp_terms = np.exp(np.outer(lam, t))          # (n_modes, N_t)
            x_t = VR @ (exp_terms * z0[:, None])
            x_t = np.real(x_t)

            # --- plot rotor angles Δδ (jointly normalized) ---
            y = x_t[delta_idx, :]
            m = np.max(np.abs(y));  y = y/m if m>0 else y
            labels = [f"G{g+1}" for g in range(G)]

            plt.figure(figsize=(9,4))
            for g in range(G):
                plt.plot(t, y[g, :], label=labels[g])
            plt.xlabel("Time [s]"); plt.ylabel("Normalized $\\Delta\\delta$ [–]")
            plt.title(f"Inter-area mode only: f={fd[k_p]:.3f} Hz, ζ={zeta[k_p]:.3f}")
            plt.grid(True, alpha=0.3); plt.legend(ncols=4, loc="upper right")
            plt.tight_layout(); plt.show()




            
            