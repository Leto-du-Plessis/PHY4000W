# Imports

using LinearAlgebra
using PyCall
using PyPlot

pygui(:tk) # Set pyploy backend to tk
pygui(true) # Enable python backend use

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

function V(x)
    """
    # Arguments
    - 'x::Float': positional argument to the potential
    # Outputs
    - 'potential:: Float': Potential V(X) evaluated at position x
    # Usage
    - Functional form should be modified to represent the potential of the system under consideration.
    """
    potential = cos(x)
    return potential
end

diagonal = ones(1, N+1)

for i = 1:N+1 
    diagonal[i] = V(xs[i])
end

U = Diagonal(vec(diagonal))

# Create Hamiltonian matrix and calculation eigenvalues and eigenfunctions
# --------------------------------------

H = T + U

eigenvalues = eigvals(H)
eigenvectors = eigvecs(H)

# Outputs

println("First 10 eigenvalues:")
println(eigenvalues[1:10])

# Plots
# Plots Settings
plotset = true

function outputplots(rangemin, rangemax, rangestep, xs, eigenvectors)

    for i = rangemin:rangestep:rangemax
        plot(xs, eigenvectors[:,i], label = "psi_" * repr(i))
    end
    vlines(0, minimum(eigenvectors), maximum(eigenvectors), color = "k")
    hlines(0, xs[1], xs[size(xs)[1]-1], color = "k")
    title("First couple of eigenfunctions for the specified potentials", size = 30)
    ylabel("psi(x)", size = 25)
    xlabel("x", size = 25)
    xticks(size = 20)
    yticks(size = 20)
    legend(fontsize = 20)

end

if plotset == true
    outputplots(1, 5, 1, xs, eigenvectors)
end
