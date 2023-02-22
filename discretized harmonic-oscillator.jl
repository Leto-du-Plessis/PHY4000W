using LinearAlgebra
using PyCall
pygui(:tk) # Set pyploy backend to tk
using PyPlot
pygui(true)

# Settings
# --------------------------------------

N = 1000   # Number of discrete points
m = 1     # Mass
omega = 1 # Angular momentum
Lbound = -5 # Left x - bound
Rbound = 5  # Right x - bound
deltax = (Rbound - Lbound)/N       # x - step
xs = LinRange(Lbound, Rbound, N+1) # Create an array of x values

# Create kinetic energy matrix
# --------------------------------------

diagonal = vec(2 * ones(1, N+1))   # vector containing the diagonal entries in the matrix
subdiagonal = vec(-1 * ones(1, N)) # vector containing the subdiagonal entries in the matrix

T = 1/(2 * m * deltax ^ 2) * SymTridiagonal(diagonal, subdiagonal) # Kinetic energy matrix

#println(T)

# Create potential energy matrix
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

# Outputs
# --------------------------------------

# Text 

println("Calculated eigenvalues:")
println(eigenvalues[1:10])
println(' ')
#println(eigenvectors)

# Plots

for i = 1:5
    plot(xs, eigenvectors[:,i], label = "psi_" * repr(i))
end
title("First five eigenfunctions of the quantum harmonic oscillator", size = 30)
ylabel("psi(x)", size = 25)
xlabel("x", size = 25)
legend()


show()
#plot(range(Lbound,Rbound,N), toBePlot[1])