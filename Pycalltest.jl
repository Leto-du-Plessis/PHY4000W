using PyCall

function __init__()
    py"""
    import numpy as np 

    def one(x):
        return np.sin(x)

    """
end

x = 5
println(py"one"(x))

