using ITensors
using Plots

# Constants
N = 3
dt = 0.1
dx = 1/(2^N-1)
A = dt/(2*pi^2*dx)
B = 1-dt/(pi^2*dx)
n_steps = 3

function get_mpo(sites)

    # Define links indices and initialize MPO of time evolution
    links = [Index(3, "Link,l=$n") for n in 1:N-1]
    mpo = MPO(sites)

    # First tensor of MPO
    s1, s2, l1 = prime(sites[1]), dag(sites[1]), links[1]
    mpo[1] = ITensor(l1, s1, s2)
    # Id
    mpo[1][l1 => 1, s1 => 1, s2 => 1] = 1 
    mpo[1][l1 => 1, s1 => 2, s2 => 2] = 1
    # sigma_01 = [[0 1], [0 0]] = Up
    mpo[1][l1 => 2, s1 => 1, s2 => 2] = 1
    # sigma_10 = [[0 0], [1 0]] = Down
    mpo[1][l1 => 3, s1 => 2, s2 => 1] = 1

    # Set the bulk MPO 
    for i in 2:N-1

        s1, s2, l1, l2 = prime(sites[i]), dag(sites[i]), links[i-1], links[i]
        mpo[i] = ITensor(l1, l2, s1, s2)
        # Id
        mpo[i][l1 => 1, l2 => 1, s1 => 1, s2 => 1] = 1
        mpo[i][l1 => 1, l2 => 1, s1 => 2, s2 => 2] = 1
        # Up
        mpo[i][l1 => 1, l2 => 2, s1 => 1, s2 => 2] = 1
        # Down 
        mpo[i][l1 => 1, l2 => 3, s1 => 2, s2 => 1] = 1
        # Up
        mpo[i][l1 => 1, l2 => 3, s1 => 1, s2 => 2] = 1
        # Up
        mpo[i][l1 => 3, l2 => 3, s1 => 1, s2 => 2] = 1
        # Down
        mpo[i][l1 => 2, l2 => 2, s1 => 2, s2 => 1] = 1

    end

    # Last tensor of MPO
    s1, s2, l1 = prime(sites[N]), dag(sites[N]), links[N-1]
    mpo[N] = ITensor(l1, s1, s2)
    # Id
    mpo[N][l1 => 1, s1 => 1, s2 => 1] = B
    mpo[N][l1 => 1, s1 => 2, s2 => 2] = B
    # Up
    mpo[N][l1 => 1, s1 => 1, s2 => 2] = A
    # Down
    mpo[N][l1 => 1, s1 => 2, s2 => 1] = A
    # Down
    mpo[N][l1 => 2, s1 => 2, s2 => 1] = A
    # Up
    mpo[N][l1 => 3, s1 => 1, s2 => 2] = A

    return mpo

end

function decimal_to_padded_binary_list(decimal, bit_length)
    binary_list = Int[]

    while decimal > 0 || length(binary_list) < bit_length
        pushfirst!(binary_list, (decimal % 2) + 1)
        decimal = div(decimal, 2)
    end

    # Pad with leading zeros if needed
    while length(binary_list) < bit_length
        pushfirst!(binary_list, 0)
    end

    return binary_list 
end

function get_initial_MPS(sites)

    # Prepare the N-dimensional tensor to be SVDd into an MPS
    tmp = ITensor(sites)
    tmp_x = 0
    for i in 1:2^N
        # Convert a decimal into a list of integers that represent the number in binary 
        # (instead of 0 and 1 we have 1 and 2 in this list to fit with Julia)
        binary_list = decimal_to_padded_binary_list(i-1, N) 
        tmp[binary_list...] = sin(pi*(tmp_x))/pi^2
        tmp_x += dx
    end

    # Perform SVD N-1 times to convert the N-dimensional tensor into an MPS of N sites
    mps = MPS(sites)
    for i in 1:N-1
        tmp_inds = inds(tmp)
        if i == 1
            U, S, V = svd(tmp, tmp_inds[1], lefttags = "Link,l=$(i)")
        else
            U, S, V = svd(tmp, tmp_inds[1:2]..., lefttags = "Link,l=$(i)", righttags = "Link,l=$(i)")
        end
        V = S*V
        tmp = V
        mps[i] = U
    end
    mps[N] = tmp

    return mps
    
end

function time_evolve(n_steps, mpo, mps)

    res = [mps]
    for _ in 1:n_steps
        mps = apply(mpo, mps)
        push!(res, mps)
    end
    return res
end

function mps_to_list(mps)
    res = []
    for i in 1:2^N
        basis_state = decimal_to_padded_binary_list(i-1, N)
        coefficient = 1.0
        for j in 1:N
            coefficient *= mps[j][basis_state[j]]
        end
        push!(res, coefficient)
    end
    return res
end

function plot_results(mps_list)

    p = plot()
    for (t, mps) in enumerate(mps_list)
        plot!(mps_to_list(mps), label = "t = $t")
    end
    savefig("evolution.png")

end

# Initialize the site indices
sites = siteinds("S=1/2", N)

# Get the time evolution MPO
mpo = get_mpo(sites)

# Get the initial MPS according to the initial condition U(t = 0, x) = sin(pi * x)/pi^2
mps = get_initial_MPS(sites)

# Perform time evolution according to the heat equation
mps_list = time_evolve(n_steps, mpo, mps)

# Plot the results
plot_results(mps_list)

