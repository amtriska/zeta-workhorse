# zeta_workhorse

**Overview**
This library calculates graph zeta functions and spectral gaps using the graph's adjacency matrix. It can directly parse a comma-separated adjacency matrix (.CSV) or accept an undirected graph file in GraphML (.GML) format. All zeta functions treat input graphs as unweighted. Weighted graphs are supported by `tropical_trace`, `spectral_zeta_val`, `fiedler_value`, and `laplacian_gap`, which operate directly on the Laplacian or adjacency matrix.

**Relevant mathematics**
The functions in this library compute the reciprocal (inverse) zeta polynomials. The specific formulas implemented are:

* **Ihara zeta (inverse)**: $(1 - u^2)^{r - 1} \det(I - uA + u^2(D - I))$
  Useful for counting non-backtracking cycles and analyzing the topological structure of a graph.
* **Bartholdi zeta (inverse)**: $(1 - (1 - t)^2 u^2)^{r - 1} \det(I - uA + u^2(D - (1 - t)I))$
  Introduces the variable $t$ to control the weight of backtracking steps. Setting $t=0$ yields the non-backtracking Ihara zeta function, and setting $t=1$ allows unlimited backtracking.
* **Bowen-Lanford zeta (inverse)**: $\det(I - uA)$
  Useful for analyzing standard closed walks and path counting.
* **Spectral zeta**: $\zeta(s) = \sum \lambda_i^{-s}$
  Calculated from the non-zero eigenvalues ($\lambda_i$) of the graph Laplacian. It serves as a graph-theoretic analog to the Riemann zeta function and is used to compute structural invariants.
* **Tropical trace**: $Z_{trop}(k) = \min_{v} (A^{\otimes k})_{v,v}$
  Evaluated using min-plus or max-plus algebra (depending on whether min or max mode is selected). By calculating tropical matrix powers ($A^{\otimes k}$), this function determines the minimum/maximum cost of cycles up to a specified length for network optimization.
* **Ihara gap**: $|u_1| - |u_0|$
  Isolates the roots of the reciprocal Ihara zeta polynomial, determines their absolute values (magnitudes), sorts them in ascending order, and finds the difference between the two smallest distinct values. The gap is useful for measuring expansion properties based specifically on non-backtracking walks. Networks with an optimally large Ihara gap are known as Ramanujan graphs, which represent the theoretical ideal of a network being completely sparse yet highly connected.
* **Fiedler value**
  The smallest non-zero eigenvalue ($\lambda_2$) of the graph Laplacian. It measures algebraic connectivity. A value approaching zero indicates a structural bottleneck (i.e. the network can easily be cut into two disconnected subgraphs).
* **Laplacian gap (eigengap)**
  The difference between consecutive non-zero eigenvalues ($\lambda_{k+1} - \lambda_k$). This gap is used as a heuristic in spectral clustering to identify the optimal number of communities within a network. A disproportionately large gap at a specific $k$ suggests the data naturally separates into exactly $k$ dense clusters.

**A note on `_poly` vs. `_roots` functions**
Each zeta function is available in two variants. The `_poly` variants use SymPy's exact symbolic eigensolver to return an exact polynomial or rational expression. The `_roots` variants use NumPy's numerical eigenvalue routines to directly return the roots as a NumPy array of complex numbers, which is much faster and scales to larger graphs. For most practical applications, the `_roots` variants are recommended. The `_poly` variants are best suited for small graphs where exact symbolic output is needed (e.g. for further symbolic manipulation or verification).

**Performance and graph size**
**`_poly` functions** (`ihara_zeta_poly`, `bartholdi_zeta_poly`, `bowen_lanford_zeta_poly`, `spectral_zeta_poly`) compute a symbolic determinant of a $V \times V$ matrix, which runs in $O(V^4)$ time with symbolic arithmetic. These functions become impractical on graphs with more than approximately 30–40 nodes and will hang indefinitely on larger inputs.

**`_roots` functions and eigenvalue-based functions** (`ihara_zeta_roots`, `bartholdi_zeta_roots`, `bowen_lanford_zeta_roots`, `spectral_zeta_val`, `fiedler_value`, `laplacian_gap`) use dense matrix eigenvalue routines running in $O(V^3)$ time and require $O(V^2)$ memory. The Ihara and Bartholdi `_roots` functions operate on a $2V \times 2V$ companion matrix, making them approximately eight times more memory-intensive than the others. These functions become impractical on graphs with more than a few thousand nodes on typical hardware.

**`tropical_trace`** runs in O(max_k · V³) time. For large graphs, keep max_k as small as possible.

**Installation**
```bash
pip install zeta-workhorse
```

**Quickstart**
```python
from zeta_workhorse import load_gml, ihara_zeta_roots, ihara_gap

# Load an unweighted network
data = load_gml("karate.gml", data_type=float)

# Calculate the roots of the Ihara zeta polynomial numerically
roots = ihara_zeta_roots(data)

# Calculate the Ihara spectral gap from the roots
gap = ihara_gap(roots)
print(gap)
```

**API Reference**

**Parsers**
* **`load_csv(file_path, data_type)`**: Requires `file_path` (string) and `data_type` (type). Returns a NumPy array.
* **`load_gml(gml_file_path, data_type, weight_field=None)`**: Requires `gml_file_path` (string) and `data_type` (type). If `weight_field` is omitted or `None`, the graph is treated as unweighted (all edge weights set to 1). If `weight_field` is provided as a string (e.g. `weight_field="value"`), that named edge attribute is used as the weight. Returns a NumPy array.

**Ihara zeta**
* **`ihara_zeta_poly(A)`**: Requires `A` (NumPy array). Treats the graph as unweighted. Returns an exact SymPy polynomial in terms of `u`. Not recommended for graphs with more than ~30–40 nodes.
* **`ihara_zeta_roots(A)`**: Requires `A` (NumPy array). Treats the graph as unweighted. Returns the roots of the Ihara zeta polynomial as a NumPy array of complex numbers, computed numerically via a $2V \times 2V$ companion matrix.

**Bartholdi zeta**
* **`bartholdi_zeta_poly(A, t_val=None)`**: Requires `A` (NumPy array). Treats the graph as unweighted. Returns an exact SymPy polynomial in terms of `u` and `t`. If `t_val` (float or integer) is provided, substitutes that value for `t` and returns a polynomial in `u` only. Setting `t_val=0` recovers the Ihara zeta. Not recommended for graphs with more than ~30–40 nodes.
* **`bartholdi_zeta_roots(A, t_val)`**: Requires `A` (NumPy array) and `t_val` (float or integer). Treats the graph as unweighted. Returns the roots as a NumPy array of complex numbers, computed numerically via a $2V \times 2V$ companion matrix.

**Bowen-Lanford zeta**
* **`bowen_lanford_zeta_poly(A)`**: Requires `A` (NumPy array). Returns an exact SymPy polynomial in terms of `u`. Not recommended for graphs with more than ~30–40 nodes.
* **`bowen_lanford_zeta_roots(A)`**: Requires `A` (NumPy array). Returns the roots as a NumPy array of complex numbers, computed from the eigenvalues of `A`.

**Spectral zeta**
* **`spectral_zeta_poly(A)`**: Requires `A` (NumPy array). Returns an exact SymPy expression in terms of `s`, with integer coefficients representing eigenvalue multiplicities. Uses SymPy's symbolic eigensolver; unlike the other `_poly` functions, the output is not a polynomial but a sum of power terms $\lambda_i^{-s}$. Not recommended for graphs with more than ~30–40 nodes.
* **`spectral_zeta_val(A, s)`**: Requires `A` (NumPy array) and `s` (float or integer). Evaluates the spectral zeta numerically at the given value of `s`. Returns a float.

**Gap and connectivity measures**
* **`ihara_gap(roots)`**: Requires `roots` (NumPy array of complex numbers, as returned by `ihara_zeta_roots`). Returns a float representing the difference between the two smallest distinct root magnitudes.
* **`fiedler_value(A)`**: Requires `A` (NumPy array). Returns a float representing the Fiedler value — the smallest non-zero eigenvalue of the graph Laplacian ($\lambda_2$), which measures algebraic connectivity. Raises a `ValueError` if the graph is disconnected. Filters eigenvalues smaller than 1e-10.
* **`laplacian_gap(A, k)`**: Requires `A` (NumPy array) and `k` (integer, must be ≥ 2). Returns a float representing the eigengap $\lambda_{k+1} - \lambda_k$ between the $k$-th and $(k+1)$-th smallest non-zero Laplacian eigenvalues. Raises a `ValueError` if the graph is disconnected or if `k` exceeds the number of available eigenvalues. Filters eigenvalues smaller than 1e-10.

**Tropical trace**
* **`tropical_trace(A, max_k, mode="min")`**: Requires `A` (NumPy array) and `max_k` (integer). The optional `mode` parameter accepts `"min"` (default, min-plus algebra, finds minimum-cost cycles) or `"max"` (max-plus algebra, finds maximum-cost cycles). Evaluates the tropical trace sequence from k=1 to k=`max_k`. Structural zeros are automatically converted to the appropriate identity element (+∞ for min, −∞ for max). Returns a list of floats. Runtime scales as O(`max_k` · V³); keep `max_k` small for large graphs.