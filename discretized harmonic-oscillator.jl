using LinearAlgebra
Pkg,add("Plots")
using Plots

# Settings
# --------------------------------------

N = 10000    # Number of discrete points
m = 1       # Mass
omega = 1 
Lbound = -3 # Left x - bound
Rbound = 3  # Right x - bound
deltax = (Rbound - Lbound)/N # x - step

# Create kinetic energy matrix
# --------------------------------------

diagonal = vec(2 * ones(1, N+1))
subdiagonal = vec(-1 * ones(1, N))

T = 1/(2 * m * deltax ^ 2) * SymTridiagonal(diagonal, subdiagonal) # Kinetic energy matrix

#println(T)

# Create potential enery matrix
# --------------------------------------

diagonal = ones(1, N+1)
x = Lbound
for i = 1:N+1 
    diagonal[i] = x ^ 2
    global x = x + deltax
end

U = ((m * omega ^ 2)/2) * Diagonal(vec(diagonal))

# Create Hamiltonian matrix
# --------------------------------------

H = T + U # Hamiltonian matrix

# Calculate eigenvalues and eigenfunctions

eigenvalues = eigvals(H)
eigenvectors = eigvecs(H)

println(eigenvalues[1:10])
println(' ')
#println(eigenvectors)