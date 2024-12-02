import numpy as np

# Example matrix (A)
A = np.array([[3, 2], [2, 3]])

# Perform Singular Value Decomposition
U, S, Vt = np.linalg.svd(A)

# U, S, and Vt are the decomposition results:
# U: Left singular vectors
# S: Singular values (as a 1D array)
# Vt: Right singular vectors (transpose of V)

print("Matrix A:")
print(A)

print("\nLeft singular vectors (U):")
print(U)

print("\nSingular values (S):")
print(S)

print("\nRight singular vectors (Vt):")
print(Vt)
