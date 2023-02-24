# Imports

using LinearAlgebra
using PyCall
using PyPlot

# Settings
# --------------------------------------

N = 1000  # Number of discrete points
m = 1     # Mass
omega = 1 # Angular velocity
Lbound = -5 # Left x - bound
Rbound = 5  # Right x - bound
deltax = (Rbound - Lbound)/N       # x - step
xs = LinRange(Lbound, Rbound, N+1) # Create an array of x values

# Create kinetic energy matrix
# --------------------------------------

diagonal = vec(2 * ones(1, N+1))   # vector containing the diagonal entries in the matrix
subdiagonal = vec(-1 * ones(1, N)) # vector containing the subdiagonal entries in the matrix

T = 1/(2 * m * deltax ^ 2) * SymTridiagonal(diagonal, subdiagonal) # Kinetic energy matrix

# Create potential energy matrix
# --------------------------------------

function v(x)
    form = x 
    return form
end

diagonal = ones(1, N+1)

for i = 1:N+1 
    diagonal[i] = v(xs[i])
end

U = Diagonal(vec(diagonal))

# Create Hamiltonian matrix and calculation eigenvalues and eigenfunctions
# --------------------------------------

H = T + U

eigenvalues = eigvals(H)
eigenfunctions = eigvecs(H)

