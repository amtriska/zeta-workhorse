import numpy as np
import sympy as sp
import networkx as nx

def load_csv(file_path, data_type):
    data = np.loadtxt(fname = file_path, delimiter=',', dtype = data_type)
    if data.ndim != 2:
        raise ValueError("Error. Adjacency matrices must have two axes.")
    if data.shape[0] != data.shape[1]:
        raise ValueError("Error. Adjacency matrices must be square.")
    return data

def load_gml(gml_file_path, data_type, weight_field = None):
    target_attribute = weight_field
    G = nx.read_gml(gml_file_path, label="id")
    data = nx.to_numpy_array(G, weight = target_attribute, dtype = data_type)
    if data.ndim != 2:
        raise ValueError("Error. Adjacency matrices must have two axes.")
    if data.shape[0] != data.shape[1]:
        raise ValueError("Error. Adjacency matrices must be square.")
    return data