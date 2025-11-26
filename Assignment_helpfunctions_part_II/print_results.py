from tabulate import tabulate
import numpy as np
from pathlib import Path

def print_eigenvalues(lambda_A, f_d, zeta):
    
    column_headers = ['No.', 'Re{λ}', 'Im{λ}', 'f_d [Hz]', 'ζ']
    
    data = []
    for i in range(len(lambda_A)):
        row = [i+1, 
               f"{np.real(lambda_A[i]):.4f}",
               f"{np.imag(lambda_A[i]):.4f}" if np.imag(lambda_A[i]) != 0 else '−',
               f"{f_d[i]:.4f}" if np.imag(lambda_A[i]) != 0 else '−',
               f"{zeta[i]:.4f}" if np.imag(lambda_A[i]) != 0 else '−']
        data.append(row)
    
    table = tabulate(data, headers = column_headers, tablefmt = "fancy_grid")
    
    print('Eigenvalues, Frequency and Damping Ratio:')
    print(table)
 
    
def print_P_matrix(P, row_headers,decimals=4):

    column_headers = []
    for i in range(P.shape[1]):
        column_headers.append(f'λ{i+1}')
        
    data = np.empty(P.shape)
    for i in range(P.shape[0]):
        for j in range(P.shape[1]):
            data[i,j] = np.abs(P[i,j])
    
    table = tabulate(data, headers = column_headers, showindex = row_headers, tablefmt="fancy_grid", floatfmt = '.2f')
    
    print('Participation Matrix (Absolute Values):')
    print(table)

def build_modes_table_latex(lam, groups, descriptions, *,
                            caption="Eigenvalues, oscillatory frequencies and damping, and dominant states (Q2.1.2).",
                            label="tab:q2_modes",
                            tabcolsep_pt=6,
                            eps_imag=1e-8):
    
    lines = []
    lines.append("\\begin{table}[H]")
    lines.append("\\centering")
    lines.append("\\scriptsize")
    lines.append(f"\\setlength{{\\tabcolsep}}{{{tabcolsep_pt}pt}}")
    lines.append("\\begin{tabular}{c|rr|r|r|l}")
    lines.append("\\hline")
    lines.append("\\textbf{No.} & \\multicolumn{2}{c|}{\\textbf{Eigenvalues}} & \\textbf{Frequency (Hz)} & \\textbf{Damping Ratio} & \\textbf{Dominant States}\\\\")
    lines.append(" & Real & Imaginary &  &  & \\\\ \\hline")

    no_counter = 1
    for gi, grp in enumerate(groups):
        desc = descriptions[gi] if gi < len(descriptions) else ""
        if len(grp) == 2 and abs(np.imag(lam[grp[0]])) > eps_imag:
            k = grp[0]
            sig = float(np.real(lam[k]))
            omg = float(abs(np.imag(lam[k])))
            fd  = omg/(2*np.pi)
            zt  = -sig/np.sqrt(sig**2 + omg**2)
            sig_s = f"{sig:.3f}"
            omg_s = f"{omg:.3f}"
            no_str = f"{no_counter},{no_counter+1}"
            lines.append(f"{no_str} & {sig_s} & $\\pm{omg_s}$ & {fd:.3f} & {zt:.3f} & {desc}\\\\")
            no_counter += 2
        else:
            k = grp[0]
            sig = float(np.real(lam[k]))
            sig_s = f"{sig:.3f}"
            lines.append(f"{no_counter} & {sig_s} & -- & -- & -- & {desc}\\\\")
            no_counter += 1

    lines.append("\\hline")
    lines.append("\\end{tabular}")
    lines.append(f"\\caption{{{caption}}}")
    lines.append(f"\\label{{{label}}}")
    lines.append("\\end{table}")
    return "\n".join(lines)


def write_tex_file(path, text):
    Path(path).write_text(text, encoding="utf-8", newline="\n")
    