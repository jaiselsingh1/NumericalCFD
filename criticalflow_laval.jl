using LinearAlgebra
using Plots

function laval_nozzle_numerical_solver(u_outlet, ix=210)
    # Parameters
    # iter = 0
    # max_iter = 500
    ε = 1e-3
    Δx = 1/ix
    Δt = Δx/(-u_outlet)
    x = range(0, 1, length=ix+1)
    
    # Laval Nozzle Geometry
    g = zeros(ix+1)
    g′ = zeros(ix+1)
    
    for i in 1:ix+1
        ξ = x[i]
        g[i] = 1 + 4*(ξ - 0.5)^2
        g′[i] = 8*(ξ - 0.5)
    end
    
    # Initial Conditions
    u = zeros(ix+1)
    u_prev = copy(u)
    
    # Boundary Conditions
    u[1] = u_outlet
    u[ix+1] = u_outlet
    

    converged = false
    while !converged
        u_prev .= u
        
        # Create new solution array
        u_next = zeros(ix+1)
        u_next[1] = u_outlet
        u_next[ix+1] = u_outlet
        
        for i in 2:ix
            u_i₋½ = (u[i-1] + u[i])/2
            u_i₊½ = (u[i+1] + u[i])/2
            
            RHS = (1-(u[i]^2)/2)*g′[i]/g[i]
            
            if (u_i₋½ > 0) && (u_i₊½ > 0)
                # Supersonic point
                u_next[i] = u[i] + Δt * (RHS - u_i₋½*(u[i] - u[i-1])/Δx)
            elseif (u_i₋½ > 0) && (u_i₊½ < 0)
                # Shock point
                u_next[i] = u[i] + Δt * (RHS - u_i₋½*(u[i] - u[i-1])/Δx - 
                                         u_i₊½*(u[i+1] - u[i])/Δx)
            elseif (u_i₋½ ≤ 0) && (u_i₊½ ≥ 0)
                # Sonic point
                u_next[i] = (u[i] + Δt*RHS)/(1 + Δt/(2*Δx)*(u[i+1] - u[i-1]))
            elseif (u_i₋½ < 0) && (u_i₊½ < 0)
                # Subsonic point
                u_next[i] = u[i] + Δt * (RHS - u_i₊½*(u[i+1] - u[i])/Δx)
            end
        end
        
        # Update solution
        u .= u_next
        # iter += 1 
        
        # Check for convergence
        if norm(u - u_prev) < ε
            converged = true
        end

        #=
        if iter == max_iter
            converged = true 
        end
        =# 
    end
    
    return u, x
end

function laval_nozzle_exact_solver(u_outlet, ix=210)
    Δx = 1/ix
    x = range(0, 1, length=ix+1)
    u = zeros(ix+1)
    
    for i in 1:ix+1
        ξ = x[i]
        u[i] = -sqrt(2 - (4 - 2*u_outlet^2)/(1 + 4*(ξ - 0.5)^2))
    end
    
    return u, x
end

# Run both cases
ix = 210
u_num_a, x = laval_nozzle_numerical_solver(-sqrt(1.5), ix)
u_exact_a, _ = laval_nozzle_exact_solver(-sqrt(1.5), ix)
u_num_b, _ = laval_nozzle_numerical_solver(-1.0, ix)
u_exact_b, _ = laval_nozzle_exact_solver(-1.0, ix)

# Compare numerical and exact solutions
error1 = maximum(abs.(u_num_a - u_exact_a))
error2 = maximum(abs.(u_num_b - u_exact_b))
println("Maximum error for u_outlet = -sqrt(1.5): $error1")
println("Maximum error for u_outlet = -1.0: $error2")

# Create plots
p = plot(
    layout = (2, 1),
    size = (800, 600),
    legend = :topright
)

# Plot for case 1: u_outlet = -sqrt(1.5)
plot!(
    p[1],
    x,
    [u_num_a, u_exact_a],
    label = ["Numerical" "Exact"],
    title = "u_outlet = -sqrt(1.5)",
    xlabel = "x",
    ylabel = "u"
)

# Plot for case 2: u_outlet = -1.0
plot!(
    p[2],
    x,
    [u_num_b, u_exact_b],
    label = ["Numerical" "Exact"],
    title = "u_outlet = -1.0",
    xlabel = "x",
    ylabel = "u"
)

display(p)