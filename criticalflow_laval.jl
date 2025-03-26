using LinearAlgebra
using Plots

using LinearAlgebra
using Plots

function solve_laval_nozzle(u_exit)
    # Setup parameters
    ix = 211
    Δx = 1/(ix-1)
    x = range(0, 1, length=ix) # Create x grid points
    
    # Initialize solution
    u = zeros(ix)
    u[1] = u_exit # Set inlet boundary condition
    u[ix] = u_exit # Set outlet boundary condition
    
    # set all interior points to 0 
    for i in 2:ix-1
        u[i] = 0.0
    end
    
    # Functions for nozzle geometry
    function g(x)
        return 1 + 4*(x-0.5)^2
    end
    
    function dg(x)
        return 8*(x-0.5)
    end
    
    # Functions to calculate half-point values
    function u_halffwd(u, i)
        return (u[i+1] + u[i])/2
    end
    
    function u_halfbwd(u, i)
        return (u[i] + u[i-1])/2
    end
    
    # Main time-stepping loop
    uold = copy(u)
    converged = false
    
    while !converged
        # Calculate time step for stability
        Δt = Δx/(-u_exit)
        
        # Save current solution 
        uold .= u
        
        # Update solution for interior points
        for i in 2:ix-1
            # Calculate RHS (right-hand side) term
            rhs = (1 - (u[i])^2/2) * dg(x[i])/g(x[i])
            
            # Determine flow regime and apply appropriate scheme
            if (u_halfbwd(u, i) > 0) && (u_halffwd(u, i) > 0)
                # Supersonic point
                u[i] = u[i] + Δt * (rhs - u_halfbwd(u, i) * (u[i] - u[i-1])/Δx)
            elseif (u_halfbwd(u, i) > 0) && (u_halffwd(u, i) < 0)
                # Shock point
                u[i] = u[i] + Δt * (rhs - u_halfbwd(u, i) * (u[i] - u[i-1])/Δx - 
                                   u_halffwd(u, i) * (u[i+1] - u[i])/Δx)
            elseif (u_halfbwd(u, i) < 0) && (u_halffwd(u, i) > 0)
                # Sonic point
                u[i] = (u[i] + Δt * rhs) / (1 + Δt/(2*Δx) * (u[i+1] - u[i-1]))
            elseif (u_halfbwd(u, i) < 0) && (u_halffwd(u, i) < 0)
                # Subsonic point
                u[i] = u[i] + Δt * (rhs - u_halffwd(u, i) * (u[i+1] - u[i])/Δx)
            end
        end
        
        if maximum(abs.(u - uold)) <= 1.0e-6
            converged = true
        end
    end
    
    return u, x
end

# Solve for the first case: u_ix = -sqrt(1.5)
u1, x = solve_laval_nozzle(-sqrt(1.5))

# Solve for the second case: u_ix = -1
u2, x = solve_laval_nozzle(-1)

# Function to calculate exact solution
function u_exact(x, u_boundary)
    c = 4 - 2*u_boundary^2
    return -sqrt(2 - c/(1 + 4*(x-0.5)^2))
end

# Calculate exact solutions for comparison
u1_exact = [u_exact(xi, -sqrt(1.5)) for xi in x]
u2_exact = [u_exact(xi, -1) for xi in x]

# Compare numerical and exact solutions
error1 = maximum(abs.(u1 - u1_exact))
error2 = maximum(abs.(u2 - u2_exact))
println("Maximum error for u_ix = -sqrt(1.5): $error1")
println("Maximum error for u_ix = -1: $error2")

# Create subplot layout with 2 rows, 1 column
p = plot(
    layout = (2, 1),
    size = (800, 600),
    legend = :topright
)

# Plot for case 1: u_ix = -sqrt(1.5)
plot!(
    p[1],
    x,
    [u1, u1_exact],
    label = ["Numerical" "Exact"],
    title = "u_outlet = -sqrt(1.5)",
    xlabel = "x",
    ylabel = "u"
)

# Plot for case 2: u_ix = -1
plot!(
    p[2],
    x,
    [u2, u2_exact],
    label = ["Numerical" "Exact"],
    title = "u_outlet = -1.0",
    xlabel = "x",
    ylabel = "u"
)
display(p)