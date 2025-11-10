import numpy as np
import re

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
    
    # normalize participation matrix
    P_abs = np.abs(P)
    max_abs_per_col = P_abs.max(axis=0)
    P = P / max_abs_per_col
    
    return P, Phi, Psi

#- - - - - - - - - - - functions for PART 2.1 - - - - - - - - - - - 
def get_state_names(data, n_states):
    for key in ('latex_names_q2a', 'StateName', 'names_q2a', 'state_names', 'names'):
        if key in data:
            arr = np.squeeze(data[key])
            try:
                lst = arr.tolist() if hasattr(arr, "tolist") else arr
                if isinstance(lst, (list, tuple, np.ndarray)):
                    return [str(x) for x in lst]
                return [str(lst)]
            except Exception:
                pass
    return [f"x{i+1}" for i in range(n_states)]


"""Scale columns so that VL[:,k]^H VR[:,k] = 1"""
def biorthonormalize(VL, VR):
    for k in range(VR.shape[1]):
        den = VL[:, k].conj().T @ VR[:, k]
        if den != 0:
            VL[:, k] /= den.conj()
    return VL, VR

"""Group real modes [i] and complex-conjugate pairs [i,j]."""
def pair_modes(lam, tol_imag=1e-8):
    lam = np.asarray(lam)
    n = len(lam)
    used = np.zeros(n, dtype=bool)
    groups = []
    for i in range(n):
        if used[i]:
            continue
        if abs(lam[i].imag) <= tol_imag:
            groups.append([i])
            used[i] = True
        else:
            cj = lam[i].conjugate()
            cand = np.where(~used)[0]
            # choose closest to conjugate (among unused)
            j = cand[np.argmin(np.abs(lam[cand] - cj))]
            groups.append([i, int(j)])
            used[i] = True
            used[int(j)] = True
    return groups

def gen_of(name):
    _re_G = re.compile(r'G[_\-]?(\d+)')  # matches G1, G_1, G-1
    m = _re_G.search(str(name))
    return m.group(1) if m else None

def all_gen_ids(state_names):
    ids = set()
    for nm in state_names:
        g = gen_of(nm)
        if g:
            ids.add(int(g))
    return sorted(ids)

def cat_of(name):
    s = str(name).replace(' ', '')
    # δ / ω (accept LaTeX \Delta, Greek Δ, plain words)
    if r'\Delta\delta' in s or 'Δδ' in s or 'delta' in s:
        return r'$\Delta \delta$'
    if r'\Delta\omega' in s or 'Δω' in s or 'omega' in s:
        return r'$\Delta \omega$'
    # field flux linkages (generic)
    if 'psi_' in s or 'ψ' in s:
        return r'$\Delta \psi$'
    # exciter voltage
    if 'v_{m' in s or 'v_m' in s or 'vm' in s:
        return r'$\Delta v_{m,\mathrm{exc}}$'
    return r'$x$'

"""Return top-m generators by total participation in Pg (global per mode)."""
def top_generators(Pg, state_names, m=2):
    sums = {}
    for i, nm in enumerate(state_names):
        g = gen_of(nm)
        if g:
            sums[g] = sums.get(g, 0.0) + float(Pg[i])
    if not sums:
        return []
    order = sorted(sums.items(), key=lambda kv: (-kv[1], int(kv[0])))
    return [g for g, _ in order[:m]]


def describe_modes(lam, groups, P, state_names, top_k=5):
    G_IDS = all_gen_ids(state_names)
    G_total = len(G_IDS)

    descs = []
    for grp in groups:
        Pg = P[:, grp].mean(axis=1) if len(grp) == 2 else P[:, grp[0]]
        idx = np.argsort(Pg)[-top_k:][::-1]

        cat_sets = {}
        cat_scores = {}
        for i in idx:
            c = cat_of(state_names[i])
            g = gen_of(state_names[i])
            cat_sets.setdefault(c, set())
            if g:
                cat_sets[c].add(g)
            cat_scores[c] = cat_scores.get(c, 0.0) + float(Pg[i])

        cats_top = sorted(cat_scores, key=cat_scores.get, reverse=True)[:2]

        parts = []
        for c in cats_top:
            gs = sorted(cat_sets.get(c, set()), key=int)
            if G_total > 0 and len(gs) == G_total:
                who = "all"                       
            elif gs:
                who = "G" + ", G".join(gs)
            else:
                who = ""
            parts.append(c + (f" from {who}" if who else ""))

        desc = " and ".join(parts) if parts else "Dominant states not identified"
        descs.append(desc)

    return descs


import re
import numpy as np

def _normalize_name(s: str) -> str:
    # strip LaTeX/greek and unify
    s = str(s)
    s = s.replace("\\", "").replace("{","").replace("}","").replace(",","").replace(" ", "")
    s = s.replace("Δ", "delta").replace("δ","delta").replace("ω","omega").lower()
    s = s.replace("delta", "delta").replace("omega","omega")
    return s

def build_delta_omega_indices(state_names, states_per_gen=7):
    """
    Try to find Δδ and Δω by parsing names. If not found, fall back to
    [0, 7, 14, ...] for Δδ and [1, 8, 15, ...] for Δω.
    """
    norm = [_normalize_name(nm) for nm in state_names]

    # regex: look for 'delta'/'omega' and a generator tag like g1, g2, ...
    has_g = [re.search(r"g(\d+)", s) for s in norm]
    delta_idx = [i for i, s in enumerate(norm) if ("delta" in s and re.search(r"g\d+", s))]
    omega_idx = [i for i, s in enumerate(norm) if ("omega" in s and re.search(r"g\d+", s))]

    # fallback if nothing matched
    if not delta_idx or not omega_idx:
        G = len(state_names) // states_per_gen
        delta_idx = [0 + states_per_gen*g for g in range(G)]
        omega_idx = [1 + states_per_gen*g for g in range(G)]

    return delta_idx, omega_idx

def _gens_from_idx(idxs, state_names):
    """
    Given indices of states for each generator (e.g., δ or ω),
    return indices sorted by generator number and the labels ['G1','G2',...].
    """
    gens = []
    labels = []
    for i in idxs:
        s = _normalize_name(state_names[i])
        m = re.search(r"g(\d+)", s)
        if m:
            g = int(m.group(1))
        else:
            g = len(labels) + 1  # fallback numbering
        gens.append(g)
        labels.append(f"G{g}")
    order = np.argsort(gens)
    return [idxs[o] for o in order], [labels[o] for o in order]
