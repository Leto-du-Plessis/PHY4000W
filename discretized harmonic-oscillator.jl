import LinearAlgebra.Tridiagonal as Tridiagonal
import LinearAlgebra.Diagonal as Diagonal

# Settings
# --------------------------------------

N = 1000    # Number of discrete points
m = 1       # Mass
omega = 1 
Lbound = -3 # Left x - bound
Rbound = 3  # Right x - bound
deltax = (Rbound - Lbound)/N # x - step

# Create kinetic energy matrix
# --------------------------------------

diagonal = vec(2 * ones(1, N))
subdiagonal = vec(-1 * ones(1, N-1))

T = 1/(2 * m * deltax ^ 2) * Tridiagonal(subdiagonal, diagonal, subdiagonal) # Kinetic energy matrix

# Create potential enery matrix
# --------------------------------------

diagonal = ones(1, N)
x = Lbound
for i = 1:N 
    diagonal[i] = x ^ 2
    global x = x + deltax
end

U = ((m * omega ^ 2)/2) * Diagonal(vec(diagonal))

println(U)