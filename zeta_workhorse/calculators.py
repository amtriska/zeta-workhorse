import numpy as np
import sympy as sp
import networkx as nx

def ihara_zeta_poly(A):
    unweighted_A = np.where(A != 0, 1, 0)
    x = np.sum(unweighted_A, axis=1)
    D = np.diag(x)
    V = A.shape[0]
    E = np.sum(unweighted_A) // 2
    r = sp.Integer(E - V + 1)
    I = np.identity(V, dtype=int)
    u = sp.symbols('u')
    M = I - (u * unweighted_A) + ((u ** 2) * (D - I))
    det_M = sp.Matrix(M).det()
    inverse_zeta = (1 - u**2)**(r - 1) * det_M
    return inverse_zeta

def ihara_zeta_roots(A):
    unweighted_A = np.where(A != 0, 1, 0)
    x = np.sum(unweighted_A, axis=1)
    D = np.diag(x)
    V = unweighted_A.shape[0]
    I = np.identity(V, dtype=int)
    zero_matrix = np.zeros_like(I)
    B = np.block([
        [unweighted_A, I - D],
        [I, zero_matrix]
    ])
    eigenvalues = np.linalg.eigvals(B)
    roots = 1.0 / eigenvalues[eigenvalues != 0]
    return roots
    
def bartholdi_zeta_poly(A, t_val=None):
    unweighted_A = np.where(A != 0, 1, 0)
    x = np.sum(unweighted_A, axis=1)
    D = np.diag(x)
    V = unweighted_A.shape[0]
    E = np.sum(unweighted_A) // 2
    r = sp.Integer(E - V + 1) 
    I = np.identity(V, dtype=int)
    u, t = sp.symbols('u t')
    M = I - (u * unweighted_A) + ((u**2) * (D - ((1 - t) * I)))
    det_M = sp.Matrix(M).det()
    inverse_zeta_b = ((1 - ((1 - t)**2) * (u**2))**(r - 1)) * det_M
    if t_val is not None:
        inverse_zeta_b = inverse_zeta_b.subs(t, t_val)
    return sp.expand(inverse_zeta_b)
    
def bartholdi_zeta_roots(A, t_val):
    unweighted_A = np.where(A != 0, 1, 0)
    x = np.sum(unweighted_A, axis=1)
    D = np.diag(x)
    V = unweighted_A.shape[0]
    I = np.identity(V, dtype=int)
    zero_matrix = np.zeros_like(I)
    B = np.block([
        [unweighted_A, (1 - t_val) * (I - D)],
        [I, zero_matrix]
    ])
    eigenvalues = np.linalg.eigvals(B)
    bart_roots = 1.0 / eigenvalues[eigenvalues != 0]
    return bart_roots

def spectral_zeta_poly(A):
    x = np.sum(A, axis=1)
    D = np.diag(x)
    L = sp.Matrix((D - A).astype(int))
    eigenvals = L.eigenvals()
    s = sp.symbols('s')
    zeta_spec_poly = sp.Integer(0)
    for val, mult in eigenvals.items():
        if val != 0:
            zeta_spec_poly += mult * val**(-s)
    return zeta_spec_poly
    
def spectral_zeta_val(A, s):
    x = np.sum(A, axis=1)
    D = np.diag(x)
    L = D - A
    eigenvals = np.linalg.eigvals(L)  
    zeta_spec_val = 0
    for val in eigenvals:
        if abs(val) > 1e-10:
            zeta_spec_val += val**(-s)
    return zeta_spec_val

def bowen_lanford_zeta_roots(A):
    eigenvals = np.linalg.eigvals(A)
    inverse_zeta_bl_roots = 1.0 / eigenvals[eigenvals != 0]    
    return inverse_zeta_bl_roots
    
def bowen_lanford_zeta_poly(A):
    V = A.shape[0]
    I = np.identity(V, dtype = int)
    u = sp.symbols('u')    
    M = I - (u * A)
    inverse_zeta_bl_poly = sp.Matrix(M).det()
    return inverse_zeta_bl_poly

def ihara_gap(roots):
    unique_mags = np.unique(np.round(np.abs(roots),decimals=5))
    if len(unique_mags) < 2:
        raise ValueError("Not enough distinct roots to calculate a spectral gap.")
    gap = unique_mags[1] - unique_mags[0]
    return float(gap)

def fiedler_value(A):
    x = np.sum(A, axis = 1)
    D = np.diag(x)
    L = D - A
    eigenvals = np.linalg.eigvals(L)
    zero_eigenvals = eigenvals[abs(eigenvals) < 1e-10]
    if len(zero_eigenvals) > 1:
        raise ValueError("Error. Graph disconnected.")
    filtered_eigenvals = eigenvals[eigenvals > 1e-10]
    sorted_values = np.sort(filtered_eigenvals)
    if len(sorted_values) < 1:
        raise ValueError("Error. Not enough non-zero eigenvalues to calculate a spectral gap.")
    fiedler = sorted_values[0]
    return float(fiedler)

def laplacian_gap(A, k):
    x = np.sum(A, axis = 1)
    D = np.diag(x)
    L = D - A
    eigenvals = np.linalg.eigvals(L)
    zero_eigenvals = eigenvals[abs(eigenvals) < 1e-10]
    if len(zero_eigenvals) > 1:
        raise ValueError("Error. Graph disconnected.")
    filtered_eigenvals = eigenvals[eigenvals > 1e-10]
    sorted_values = np.sort(filtered_eigenvals)
    if len(sorted_values) < 2:
        raise ValueError("Error. Not enough non-zero eigenvalues to calculate a spectral gap.")
    if k < 2:
        raise ValueError("Error. k must be at least 2.")
    if k-1 >= len(sorted_values):
        raise ValueError("Error. k exceeds the number of non-zero eigenvalues.")
    lap_gap = sorted_values[k-1] - sorted_values[k-2]
    return float(lap_gap)

def tropical_multiply_min(A, B):
    n = A.shape[0]
    C = np.zeros((n, n))
    for i in range(n):
        C[i, :] = np.min(A[i, :, np.newaxis] + B, axis=0)
    return C
    
def tropical_multiply_max(A, B):
    n = A.shape[0]
    C = np.zeros((n, n))
    for i in range(n):
        C[i, :] = np.max(A[i, :, np.newaxis] + B, axis=0)
    return C

def tropical_trace(A, max_k, mode = "min"):
    sequence = []
    if mode == "min":
        tropical_A = np.where(A != 0, A, np.inf)
        current_power_matrix = tropical_A
        for i in range(1, max_k + 1):
            diag = np.diag(current_power_matrix)
            sequence.append(float(np.min(diag)))
            if i < max_k:
                current_power_matrix = tropical_multiply_min(current_power_matrix, tropical_A)
    elif mode == "max":
        tropical_A = np.where(A != 0, A, -np.inf)
        current_power_matrix = tropical_A
        for i in range(1, max_k + 1):
            diag = np.diag(current_power_matrix)
            sequence.append(float(np.max(diag)))
            if i < max_k:
                current_power_matrix = tropical_multiply_max(current_power_matrix, tropical_A)
    else:
        raise ValueError('Invalid mode. Must use "min" or "max."')
    return sequence