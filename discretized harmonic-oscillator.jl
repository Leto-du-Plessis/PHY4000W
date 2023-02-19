import LinearAlgebra.Tridiagonal as Tridiagonal
import LinearAlgebra.Diagonal as Diagonal

# Settings
# --------------------------------------

N = 10    # Number of discrete points
m = 1       # Mass
omega = 1 
Lbound = -3 # Left x - bound
Rbound = 3  # Right x - bound
deltax = (Rbound - Lbound)/N # x - step

# Create kinetic energy matrix
# --------------------------------------

diagonal = vec(2 * ones(1, N+1))
subdiagonal = vec(-1 * ones(1, N))

T = 1/(2 * m * deltax ^ 2) * Tridiagonal(subdiagonal, diagonal, subdiagonal) # Kinetic energy matrix

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

println(eigenvalues)