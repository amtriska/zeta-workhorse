import numpy as np
import sympy as sp
import networkx as nx

def ihara_zeta(A):
    x = np.sum(A, axis=1)
    D = np.diag(x)
    V = A.shape[0]
    E = np.sum(x) // 2
    r = sp.Integer(E - V + 1)
    I = np.identity(V, dtype = int)
    u = sp.symbols('u')
    M = I - (u * A) + ((u ** 2) * (D - I))
    det_M = sp.Matrix(M).det()
    inverse_zeta = (1 - u**2)**(r - 1) * det_M
    return inverse_zeta
    
def bartholdi_zeta(A, t_val=None):
    x = np.sum(A, axis=1)
    D = np.diag(x)
    V = A.shape[0]
    E = np.sum(x) // 2
    r = sp.Integer(E - V + 1) 
    I = np.identity(V, dtype=int)
    u, t = sp.symbols('u t')
    M = I - (u * A) + ((u**2) * (D - ((1 - t) * I)))
    det_M = sp.Matrix(M).det()
    inverse_zeta_b = ((1 - ((1 - t)**2) * (u**2))**(r - 1)) * det_M
    if t_val is not None:
        inverse_zeta_b = inverse_zeta_b.subs(t, t_val)
    return sp.expand(inverse_zeta_b)
    
def spectral_zeta(A):
    x = np.sum(A, axis=1)
    D = np.diag(x)
    L = sp.Matrix(D - A)
    eigenvals = L.eigenvals()    
    s = sp.symbols('s')
    zeta_spec = 0
    for val, mult in eigenvals.items():
        if val != 0:
            zeta_spec += mult * (val**(-s))
    return zeta_spec

def bowen_lanford_zeta(A):
    V = A.shape[0]
    I = np.identity(V, dtype = int)
    u = sp.symbols('u')    
    M = I - (u * A)
    inverse_zeta_bl = sp.Matrix(M).det()
    return inverse_zeta_bl
    
def ihara_gap(zeta_expr):
    u = sp.symbols('u')
    raw_coeffs = sp.Poly(zeta_expr, u).all_coeffs()
    coeffs = np.array(raw_coeffs, dtype = float)
    roots_array = np.roots(coeffs)
    unique_mags = np.unique(np.round(np.abs(roots_array),decimals=5))
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

def tropical_multiply(A, B):
    min_combined = np.min(A[:, :, np.newaxis] + B[np.newaxis, :, :], axis=1)
    return min_combined

def tropical_power(A, k):
    result = A
    for _ in range(k-1):
        result = tropical_multiply(result, A)
    return result

def tropical_zeta(A, max_k):
    zeta_sequence = []
    tropical_A = np.where(A != 0, A, np.inf)
    for i in range(1, max_k + 1):
        current_power_matrix = tropical_power(tropical_A, i)
        diag = np.diag(current_power_matrix)
        zeta_sequence.append(float(np.min(diag)))
    return zeta_sequence