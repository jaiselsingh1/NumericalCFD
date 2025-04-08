using LinearAlgebra 
using Plots 

# Domain specific 
W = 1.0 
L = 1.0 

# Grid sizes
ix1 = jx1 = 50 
ix2 = jx2 = 100 

# Number of terms for the exact solution
n_terms = 50 

# Grid spacings
Δx1 = L / ix1
Δy1 = W / jx1
Δx2 = L / ix2
Δy2 = W / jx2

# Convergence criteria
ϵ = 1e-5 

function exact_solution(x, y, L, W, n_terms)
    result = 0.0 
    for n in 1:n_terms 
        # Calculate the coefficient
        coef = ((-1)^(n+1) + 1)
        
        # Only proceed if coefficient is non-zero
        if coef != 0
            term = (2.0 / π) * (coef / n)
            term *= sin(n * π * x / L)
            term *= sinh(n * π * y / L) / sinh(n * π * W / L)
            result += term
        end
    end 
    return result 
end 

function jacobi_method(ix, jx, Δx, Δy, ϵ)
    ϕ = zeros(ix+1, jx+1)
    
    # Apply boundary conditions
    ϕ[:, 1] .= 0.0         # ϕ(x,0) = 0
    ϕ[:, jx+1] .= 1.0      # ϕ(x,W) = 1
    ϕ[1, :] .= 0.0         # ϕ(0,y) = 0
    ϕ[ix+1, :] .= 0.0      # ϕ(L,y) = 0
    
    ϕ_next = copy(ϕ)
    
    # Track residuals for plotting
    residuals = Float64[]
    
    # Iteration loop
    converged = false
    
    while !converged
        # Interior points update using the Jacobi formula from Image 3
        for i in 2:ix
            for j in 2:jx
                ϕ_next[i, j] = (
                    (ϕ[i+1, j] + ϕ[i-1, j]) / (Δx^2) + 
                    (ϕ[i, j+1] + ϕ[i, j-1]) / (Δy^2)
                ) / (2 * (1/(Δx^2) + 1/(Δy^2)))
            end
        end
        
        # Calculate max residual
        max_res = maximum(abs.(ϕ_next - ϕ))
        push!(residuals, max_res)
        
        # Check convergence
        if max_res < ϵ
            converged = true
        end
        
        # Update solution for next iteration
        ϕ .= ϕ_next
    end
    
    return ϕ, residuals
end

# Solve using two different grid sizes
ϕ1, residuals1 = jacobi_method(ix1, jx1, Δx1, Δy1, ϵ)
ϕ2, residuals2 = jacobi_method(ix2, jx2, Δx2, Δy2, ϵ)

# Create x and y coordinates for both grids
x_values1 = range(0, L, length=ix1+1)
y_values1 = range(0, W, length=jx1+1)
x_values2 = range(0, L, length=ix2+1)
y_values2 = range(0, W, length=jx2+1)

# Calculate exact solution on both grids
ϕ_exact1 = [exact_solution(x, y, L, W, n_terms) for y in y_values1, x in x_values1]
ϕ_exact2 = [exact_solution(x, y, L, W, n_terms) for y in y_values2, x in x_values2]

# Plot the max residual vs iteration number
plt1 = plot(1:length(residuals1), residuals1, 
           xlabel="Iteration", ylabel="Max Residual", 
           title="Convergence: 50×50 vs 100×100 Grid",
           label="50×50 Grid", yscale=:log10, legend=:topright)
plot!(plt1, 1:length(residuals2), residuals2, 
     label="100×100 Grid")

# Plot ϕ(x,W/2) for both grids along with exact solution
mid_index1 = div(jx1, 2) + 1
mid_index2 = div(jx2, 2) + 1

plt2 = plot(x_values1, ϕ_exact1[mid_index1, :], 
           xlabel="x", ylabel="ϕ(x,W/2)", 
           title="Solution at y=W/2",
           label="Exact Solution", legend=:topright)
plot!(plt2, x_values1, ϕ1[mid_index1, :], 
     linestyle=:dash, label="Numerical (50×50)")
plot!(plt2, x_values2, ϕ2[mid_index2, :], 
     linestyle=:dot, label="Numerical (100×100)")

# Contour plot of the exact solution
plt3 = contour(x_values2, y_values2, ϕ_exact2, 
              title="Contour Plot of Exact Solution",
              xlabel="x", ylabel="y",
              color=:turbo, levels=20)

# Display the plots
plot(plt1, plt2, plt3, layout=(3,1), size=(800, 900))








