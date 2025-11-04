import numpy as np
from Assignment_helpfunctions_part_I.P_matrix_write import latex_P_matrix, excel_P_matrix
from load_data import load_data

# load data
data_q2 = load_data('q2')
# --- 1) Nombres de estados (usa latex_names_q2a si está) ---
def flatten_names(arr):
    a = np.atleast_1d(arr).squeeze()
    out = []
    for x in a:
        try:
            out.append(str(np.array(x).item()))
        except Exception:
            out.append(str(x))
    return out

names_raw = (data_q2.get('latex_names_q2a') or
             data_q2.get('names_q2a') or
             data_q2.get('StateName'))
state_names = flatten_names(names_raw)

# --- 2) Eigen-análisis y métricas de modo ---
def get_modes(A):
    lam, Phi = np.linalg.eig(A)   # vectores propios derechos
    Psi = np.linalg.inv(Phi)      # asegura Psi.T @ Phi = I (normalización correcta)
    f_d  = np.abs(np.imag(lam)) / (2*np.pi)           # Hz
    zeta = -np.real(lam) / np.abs(lam)                # razón de amortiguamiento
    order = np.argsort(zeta)                          # peor amortiguados primero
    return lam[order], Phi[:, order], Psi[order, :], f_d[order], zeta[order]

lam, Phi, Psi, f_d, zeta = get_modes(A_q2)

# --- 3) Matriz de participación ---
P = Phi * Psi.T   # Hadamard: columna i de Phi con fila i de Psi^T

# --- 4) Resumen por consola: modos y top estados por modo ---
def print_modes(lam, f_d, zeta, k=10, title="Least-damped modes"):
    print("\n" + title)
    print("-"*64)
    m = min(k, len(lam))
    for i in range(m):
        a, b = np.real(lam[i]), np.imag(lam[i])
        print(f"λ{i+1:02d}: {a:+.5f} {b:+.5f}j   f_d={f_d[i]:.3f} Hz   ζ={100*zeta[i]:.2f}%")
    print("-"*64)

def print_P_top(P, names, mode_idx=0, top=6):
    col = P[:, mode_idx]
    mag = np.abs(col)
    idx = np.argsort(mag)[::-1][:top]
    print(f"\nTop {top} participating states in mode {mode_idx+1}:")
    for r in idx:
        ang = np.angle(col[r], deg=True)
        print(f"  {names[r]:24s}  |P|={mag[r]:.3f}   ∠{ang:+.1f}°")

print_modes(lam, f_d, zeta, k=10)
for m in range(3):   # inspecciona 3 modos peor amortiguados
    print_P_top(P, state_names, mode_idx=m, top=6)

# --- 5) Exporta tabla completa a LaTeX/Excel (opcional) ---
try:
    # Como tienes latex_names_q2a (cell de MATLAB), pasa Matlab_array=True
    latex_P_matrix(P, names_raw, Matlab_array=True,
                   filename="P_q2a.tex", maximum_col=4, bold_limit=0.10)
    excel_P_matrix(P, names_raw, Matlab_array=True,
                   filename="P_q2a.xls", bold_limit=0.10)
    print("Exported P to P_q2a.tex and P_q2a.xls")
except ImportError:
    print("P_matrix_write.py not found; skipping export.")
