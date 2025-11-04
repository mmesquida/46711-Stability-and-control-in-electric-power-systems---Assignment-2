import numpy as np
import matplotlib.pyplot as plt

from load_data import load_data
from compute_results import get_eigenvalues
from compute_results import get_P_matrix
from print_results import print_eigenvalues
from print_results import print_P_matrix
import numpy as np
import cmath
from scipy.io import loadmat
from scipy.linalg import eig as la_eig
import re
from Assignment_helpfunctions_part_I.plot_phasor_diagram_sss import plot_phasors

if __name__ == "__main__":

    output = 2
    
    match output:

        case 1:
    
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

        case 2:

             # - - - /// Assignment 2.1.1 /// - - - - - - - - - - - - - - - - - - - - - 
            print('- - - /// Assignment 2.1.1 /// - - - - - - - - - -')
            
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

            # - - - /// Assignment 2.1.2 /// - - - - - - - - - - - - - - - - - - - - - 
            print('- - - /// Assignment 2.1.2 /// - - - - - - - - - -')

            state_names = None
            for key in ['latex_names_q2a','StateName','names_q2a']:
                if key in data_q2:
                    arr = np.squeeze(data_q2[key])
                    try:
                        lst = arr.tolist()
                        state_names = [str(lst)] if isinstance(lst, (str,bytes)) else [str(x) for x in lst]
                    except Exception:
                        pass
                    break

            # biorthogonalize left and right eigenvectors
            lam, VL, VR = la_eig(A_q2, left=True, right=True)
            for k in range(VR.shape[1]):
                den = VL[:,k].conj().T @ VR[:,k]
                if den != 0:
                    VL[:,k] /= den

            # participation matrix normalized by column
            P = np.abs(VR * VL.conj())
            colsum = P.sum(axis=0, keepdims=True); colsum[colsum==0] = 1.0
            P = P / colsum

            # group modes (real and complex conjugate pairs)
            eps_imag = 1e-8
            n = len(lam)
            used = np.zeros(n, dtype=bool)
            groups = []
            for i in range(n):
                if used[i]: 
                    continue
                if abs(np.imag(lam[i])) <= eps_imag:
                    groups.append([i]); used[i] = True
                else:
                    cj = np.conj(lam[i]); j = None
                    for t in range(i+1, n):
                        if (not used[t]) and abs(lam[t] - cj) <= 1e-7:
                            j = t; break
                    if j is None:
                        diffs = [abs(lam[t]-cj) if not used[t] else 1e9 for t in range(n)]
                        j = int(np.argmin(diffs))
                    groups.append([i, j]); used[i] = True; used[j] = True

                # helper functions to classify states
            def _cat_of_name(nm_str):
                s = str(nm_str).replace(' ','').replace(',','')
                if ('psi_{f' in s) or ('ψ/f' in s):   return r'$\Delta \psi_{f}$'
                if ('psi_{kd' in s) or ('ψ/kd' in s): return r'$\Delta \psi_{kd}$'
                if ('psi_{kq1' in s) or ('ψ/kq1' in s): return r'$\Delta \psi_{kq1}$'
                if ('psi_{kq2' in s) or ('ψ/kq2' in s): return r'$\Delta \psi_{kq2}$'
                if (r'\Delta\delta' in s) or ('Δδ' in s) or ('delta' in s): return r'$\Delta \delta$'
                if (r'\Delta\omega' in s) or ('Δω' in s) or ('omega' in s): return r'$\Delta \omega$'
                if ('v_{m' in s) or ('v/m' in s):     return r'$\Delta v_{m,\mathrm{exc}}$'
                return r'$x$'

            def _gen_of_name(nm_str):
                s = str(nm_str).replace(' ','').replace(',','')
                m = re.search(r'G(\d+)', s)
                return m.group(1) if m else None

            # Construct descriptions of dominant states per mode
            descriptions = []
            for grp in groups:
                Pg = P[:, grp].mean(axis=1) if len(grp)==2 else P[:, grp[0]]
                order = np.argsort(-Pg)                 # descendente
                cum = np.cumsum(Pg[order])
                upto = max(np.searchsorted(cum, 0.80)+1, min(6, len(order)))
                take = order[:upto]

                gens_by_cat = {}
                for i in take:
                    nm = state_names[i] if state_names else f"x{i+1}"
                    cat = _cat_of_name(nm)
                    g   = _gen_of_name(nm)
                    gens_by_cat.setdefault(cat, set())
                    if g: gens_by_cat[cat].add(g)

                #Choose top 2 categories by total participation
                cat_scores = {}
                for cat in list(gens_by_cat.keys()):
                    tot = 0.0
                    for idx in range(len(Pg)):
                        nm2 = state_names[idx] if state_names else f"x{idx+1}"
                        if _cat_of_name(nm2) == cat:
                            tot += Pg[idx]
                    cat_scores[cat] = tot
                chosen = sorted(gens_by_cat.keys(), key=lambda c: -cat_scores[c])[:2]

                parts = []
                for cat in chosen:
                    gens = sorted(list(gens_by_cat[cat]), key=lambda x: int(x)) if gens_by_cat[cat] else []
                    parts.append(cat + ((" de G" + ", G".join(gens)) if gens else ""))
                desc = " y ".join(parts) if parts else "Dominant states not identified"
                descriptions.append(desc)

            # Export LaTeX table
            lines = []
            lines.append("\\begin{table}[H]")
            lines.append("\\centering")
            lines.append("\\scriptsize")
            lines.append("\\setlength{\\tabcolsep}{6pt}")
            lines.append("\\begin{tabular}{c|rr|r|r|l}")
            lines.append("\\hline")
            lines.append("\\textbf{No.} & \\multicolumn{2}{c|}{\\textbf{Eigenvalues}} & \\textbf{Frequency (Hz)} & \\textbf{Damping Ratio} & \\textbf{Dominant States}\\\\")
            lines.append(" & Real & Imaginary &  &  & \\\\ \\hline")

            no_counter = 1
            for gi, grp in enumerate(groups):
                if len(grp)==2 and abs(np.imag(lam[grp[0]]))>eps_imag:
                    k = grp[0]
                    sig = float(np.real(lam[k]))
                    omg = float(abs(np.imag(lam[k])))
                    fd  = omg/(2*np.pi)
                    zt  = -sig/np.sqrt(sig**2 + omg**2)
                    sig_s = f"{sig:.3f}"
                    omg_s = f"{omg:.3f}"
                    no_str = f"{no_counter},{no_counter+1}"
                    lines.append(f"{no_str} & {sig_s} & $\\pm{omg_s}$ & {fd:.3f} & {zt:.3f} & {descriptions[gi]}\\\\")
                    no_counter += 2
                else:
                    k = grp[0]
                    sig = float(np.real(lam[k])); sig_s = f"{sig:.3f}"
                    lines.append(f"{no_counter} & {sig_s} & -- & -- & -- & {descriptions[gi]}\\\\")
                    no_counter += 1

            lines.append("\\hline")
            lines.append("\\end{tabular}")
            lines.append("\\caption{Eigenvalues, oscillatory frequencies and damping, and dominant states (Q2.1.2).}")
            lines.append("\\label{tab:q2_modes}")
            lines.append("\\end{table}")

            with open("modes_q2_table.tex", "w", encoding="utf-8") as f:
                f.write("\n".join(lines))
            
             # - - - /// Assignment 2.1.4 /// - - - - - - - - - - - - - - - - - - - - - 
            print('- - - /// Assignment 2.1.4 /// - - - - - - - - - -')
            print("Total state names:", len(state_names))
            for i, nm in enumerate(state_names[:40]):  
                print(f"{i:02d}: {nm}")

            # Normalize names
            norm_names = []
            for nm in state_names:
                s = str(nm)
                s = s.replace("\\", "")
                s = s.replace("{", "").replace("}", "").replace(",", "")
                s = s.replace(" ", "")
                s = s.replace("Δ", "Delta")
                s = s.replace("δ", "delta")
                s = s.replace("ω", "omega")
                s = s.lower()
                norm_names.append(s)

            delta_idx = []
            omega_idx = []
            for i, s in enumerate(norm_names):
                if re.search(r'delta.*g\d+', s):
                    delta_idx.append(i)
                if re.search(r'omega.*g\d+', s):
                    omega_idx.append(i)

            print("Found delta_idx:", delta_idx)
            print("Found omega_idx:", omega_idx)

            if len(delta_idx) == 0 or len(omega_idx) == 0:
                print("No matches via name parsing -> using fallback pattern (blocks of 7 states per generator).")
                n = len(state_names)
                # try to detect number of generators for total states
                # asume 7 steates por generador
                G = n // 7
                delta_idx = [0 + 7*g for g in range(G)]       # Δδ
                omega_idx = [1 + 7*g for g in range(G)]       # Δω
                print("Fallback delta_idx:", delta_idx)
                print("Fallback omega_idx:", omega_idx)

            # Choose modes to plot: pick first 3 oscillatory modes with lowest damping ratio
            sigma = np.real(lam); omega = np.imag(lam)
            zeta = -sigma/np.sqrt(sigma**2 + omega**2)
            osc_idx = np.where(np.abs(omega) > 1e-6)[0]
            osc_idx = osc_idx[np.argsort(zeta[osc_idx])]  # from lower ζ to gratest ζ

            modes_to_plot = []
            used = set()
            for k in osc_idx:
                if k in used: 
                    continue
                modes_to_plot.append(k)
                cj = np.conj(lam[k])
                partner = None
                for t in range(len(lam)):
                    if t==k or t in used: 
                        continue
                    if abs(lam[t]-cj) < 1e-7:
                        partner = t; break
                if partner is not None:
                    used.add(partner)
                used.add(k)
            modes_to_plot = modes_to_plot[:3]
            print("Modes to plot (1-based):", [m+1 for m in modes_to_plot])

            # Aux function to get generator numbers and labels from state indices
            def _gens_from_idx(idxs):
                gens = []
                labels = []
                for i in idxs:
                    raw = state_names[i]
                    s = norm_names[i]
                    m = re.search(r'g(\d+)', s)
                    g = int(m.group(1)) if m else (len(labels)+1)
                    gens.append(g)
                    labels.append(f"G{g}")
                order = np.argsort(gens)
                return [idxs[o] for o in order], [labels[o] for o in order]

            # Plot phasor diagrams for selected modes
            for k in modes_to_plot:
                sig = float(np.real(lam[k])); omg = float(np.imag(lam[k]))
                fd_k = abs(omg)/(2*np.pi)
                zeta_k = -sig/np.sqrt(sig**2 + omg**2)

                if len(delta_idx)==0 and len(omega_idx)==0:
                    print(f"[WARN] No Δδ/Δω indices found for mode k={k+1}. Skipping plots.")
                    continue

                # Δδ
                if len(delta_idx) > 0:
                    delta_idx_sorted, L_delta = _gens_from_idx(delta_idx)
                    v = VR[:, k]
                    v_delta = v[delta_idx_sorted]
                    if v_delta.size > 0 and np.max(np.abs(v_delta))>0:
                        P_delta = np.column_stack([np.real(v_delta), np.imag(v_delta)])
                        P_delta = P_delta / np.max(np.abs(P_delta)) 
                        ang = np.arctan2(P_delta[:,1], P_delta[:,0])
                        if len(ang) >= 2:
                            ref = np.angle(np.mean(np.exp(1j*ang)))
                            grp = np.sign(np.cos(ang - ref))
                            C_delta = ['tab:blue' if g>=0 else 'tab:red' for g in grp]
                        else:
                            C_delta = ['k']*len(L_delta)
                        plt.figure()
                        plot_phasors(P_delta, C_delta, L_delta)
                        plt.gca().set_aspect('equal', 'box'); plt.grid(True, alpha=0.3)
                        mode_type = "inter-area" if (len(set(np.sign(np.cos(ang - np.mean(ang)))))>1) else "local"
                        plt.title(f"Mode k={k+1} — $\\Delta\\delta$ | f={fd_k:.3f} Hz, ζ={zeta_k:.3f} ({mode_type})")
                    else:
                        print(f"[INFO] Δδ empty or nearly zero for mode k={k+1}")
                else:
                    print(f"[INFO] No Δδ indices. Skipping Δδ plot for mode k={k+1}")

                # Δω
                if len(omega_idx) > 0:
                    omega_idx_sorted, L_omega = _gens_from_idx(omega_idx)
                    v = VR[:, k]
                    v_omega = v[omega_idx_sorted]
                    if v_omega.size > 0 and np.max(np.abs(v_omega))>0:
                        P_omega = np.column_stack([np.real(v_omega), np.imag(v_omega)])
                        P_omega = P_omega / np.max(np.abs(P_omega))
                        C_omega = ['k']*len(L_omega)
                        plt.figure()
                        plot_phasors(P_omega, C_omega, L_omega)
                        plt.gca().set_aspect('equal', 'box'); plt.grid(True, alpha=0.3)
                        plt.title(f"Mode k={k+1} — $\\Delta\\omega$ | f={fd_k:.3f} Hz, ζ={zeta_k:.3f}")
                    else:
                        print(f"[INFO] Δω empty or nearly zero for mode k={k+1}")
                else:
                    print(f"[INFO] No Δω indices. Skipping Δω plot for mode k={k+1}")

            plt.show()

        # - - - /// Assignment 2.1.5 /// - - - - - - - - - - - - - - - - - - - - - 
            "- - - - - /// Assignment 2.1.5 ///- - - - - -"
            lam, VL, VR = la_eig(A_q2, left=True, right=True)
            for k in range(VR.shape[1]):
                den = VL[:,k].conj().T @ VR[:,k]
                if den != 0:
                    VL[:,k] /= den

            sigma = np.real(lam)
            omega = np.imag(lam)
            fd    = np.abs(omega)/(2*np.pi)
            zeta  = -sigma/np.sqrt(sigma**2 + omega**2)

            G = len(state_names)//7
            delta_idx = [0 + 7*g for g in range(G)]   # Δδ: 0,7,14,21
            omega_idx = [1 + 7*g for g in range(G)]   # Δω: 1,8,15,22
            osc = np.where(np.abs(omega) > 1e-6)[0]
            k_inter = osc[np.argmin(zeta[osc])]
         
            cj = np.conj(lam[k_inter])
            k_partner = None
            for t_ in range(len(lam)):
                if t_ != k_inter and abs(lam[t_] - cj) < 1e-7:
                    k_partner = t_; break

            print(f"[2.1.5] Inter-area mode -> k={k_inter+1}, f={fd[k_inter]:.3f} Hz, ζ={zeta[k_inter]:.3f}, partner={(k_partner+1) if k_partner is not None else 'N/A'}")

            z0 = np.zeros_like(lam, dtype=complex)
            if k_partner is not None:
                z0[k_inter]   = 0.5
                z0[k_partner] = 0.5
            else:
                z0[k_inter]   = 1.0

            t_end = 12.0     # s (≈ 7 ciclos a 0.6 Hz)
            dt    = 0.01
            t     = np.arange(0.0, t_end+dt, dt)

            exp_terms = np.exp(np.outer(lam, t))         # (n_modes x N_t)
            x_t = VR @ (exp_terms * z0[:,None])          # (n_states x N_t)
            x_t = np.real(x_t)                           # estados reales

            labels = [f"G{g+1}" for g in range(G)]
            plt.figure(figsize=(9,4))
            for g, i in enumerate(delta_idx):
                y = x_t[i, :]
                ymax = np.max(np.abs(y))
                if ymax > 0: y = y / ymax
                plt.plot(t, y, label=labels[g])

            plt.xlabel("Time [s]")
            plt.ylabel("Normalized $\\Delta\\delta$ [\\,]")
            plt.title(f"Linearized response — inter-area mode only: f={fd[k_inter]:.3f} Hz, ζ={zeta[k_inter]:.3f}")
            plt.grid(True, alpha=0.3)
            h, l = plt.gca().get_legend_handles_labels()
            by_lab = {}
            for hh, ll in zip(h, l):
                if ll not in by_lab: by_lab[ll] = hh
            plt.legend(by_lab.values(), by_lab.keys(), ncols=4, loc="upper right")
            plt.tight_layout()
            plt.show()

            plt.figure(figsize=(9,4))
            for g, i in enumerate(omega_idx):
                y = x_t[i, :]
                ymax = np.max(np.abs(y))
                if ymax > 0: y = y / ymax
                plt.plot(t, y, label=labels[g])
            plt.xlabel("Time [s]")
            plt.ylabel("Normalized $\\Delta\\omega$ [\\,]")
            plt.title(f"Linearized response — inter-area mode only (speeds): f={fd[k_inter]:.3f} Hz, ζ={zeta[k_inter]:.3f}")
            plt.grid(True, alpha=0.3)
            h, l = plt.gca().get_legend_handles_labels()
            by_lab = {}
            for hh, ll in zip(h, l):
                if ll not in by_lab: by_lab[ll] = hh
            plt.legend(by_lab.values(), by_lab.keys(), ncols=4, loc="upper right")
            plt.tight_layout()
            plt.show()