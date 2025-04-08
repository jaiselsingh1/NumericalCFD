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
        term =  ((-1)^(n+1) + 1) / n 
        term *= sin(n * π * x / L)
        term *= sinh(n * π * y / L) / sinh(n * π * W / L)
        result += term
    end 
    result *= (2.0 / π)  
    return result 
end 

# Residual calculation for Jacobi method
function calculate_jacobi_residual(ϕ, ix, jx, Δx, Δy)
    max_res = 0.0
    for i in 2:ix
        for j in 2:jx
            res = (ϕ[i+1, j] - 2*ϕ[i, j] + ϕ[i-1, j])/(Δx^2) + 
                  (ϕ[i, j+1] - 2*ϕ[i, j] + ϕ[i, j-1])/(Δy^2)
            
            max_res = max(max_res, abs(res))
        end
    end
    return max_res
end

# Residual calculation for Gauss-Seidel method
function calculate_gauss_seidel_residual(ϕ, ix, jx, Δx, Δy)
    max_res = 0.0
    for i in 2:ix
        for j in 2:jx
            res = (ϕ[i+1, j] - 2*ϕ[i, j] + ϕ[i-1, j])/(Δx^2) + 
                  (ϕ[i, j+1] - 2*ϕ[i, j] + ϕ[i, j-1])/(Δy^2)
            max_res = max(max_res, abs(res))
        end
    end
    return max_res
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
        for i in 2:ix
            for j in 2:jx
                ϕ_next[i, j] = (
                    (ϕ[i+1, j] + ϕ[i-1, j]) / (Δx^2) + 
                    (ϕ[i, j+1] + ϕ[i, j-1]) / (Δy^2)
                ) / (2 * (1/(Δx^2) + 1/(Δy^2)))
            end
        end
        
        # Calculate residual specifically for Jacobi method
        max_res = calculate_jacobi_residual(ϕ_next, ix, jx, Δx, Δy)
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

function gauss_seidel(ix, jx, Δx, Δy, ϵ)
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
        for i in 2:ix
            for j in 2:jx
                ϕ_next[i, j] = (
                    (ϕ[i+1, j] + ϕ_next[i-1, j]) / (Δx^2) + 
                    (ϕ[i, j+1] + ϕ_next[i, j-1]) / (Δy^2)
                ) / (2 * (1/(Δx^2) + 1/(Δy^2)))
            end
        end
        
        # Calculate residual specifically for Gauss-Seidel method
        max_res = calculate_gauss_seidel_residual(ϕ_next, ix, jx, Δx, Δy)
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

# Plot Generation 
# Solve using two different grid sizes for Jacobi 
ϕ1_jacobi, residuals1_jacobi = jacobi_method(ix1, jx1, Δx1, Δy1, ϵ)
ϕ2_jacobi, residuals2_jacobi = jacobi_method(ix2, jx2, Δx2, Δy2, ϵ)

# Solve using two different grid sizes for Gauss-Seidel 
ϕ1_gauss, residuals1_gauss = gauss_seidel(ix1, jx1, Δx1, Δy1, ϵ)
ϕ2_gauss, residuals2_gauss = gauss_seidel(ix2, jx2, Δx2, Δy2, ϵ)

# Create x and y coordinates for both grids
x_values1 = range(0, L, length=ix1+1)
y_values1 = range(0, W, length=jx1+1)
x_values2 = range(0, L, length=ix2+1)
y_values2 = range(0, W, length=jx2+1)

# Calculate exact solution on both grids
ϕ_exact1 = [exact_solution(x, y, L, W, n_terms) for y in y_values1, x in x_values1]
ϕ_exact2 = [exact_solution(x, y, L, W, n_terms) for y in y_values2, x in x_values2]

# Get middle indices
mid_index1 = div(jx1, 2) + 1
mid_index2 = div(jx2, 2) + 1

# Plot the max residual vs iteration number for Jacobi 
plt1_jacobi = plot(1:length(residuals1_jacobi), residuals1_jacobi, 
           xlabel="Iteration", ylabel="Max Residual", 
           title="Convergence: 50×50 vs 100×100 Grid For Jacobi Method",
           label="50×50 Grid", yscale=:log10, legend=:topright)
plot!(plt1_jacobi, 1:length(residuals2_jacobi), residuals2_jacobi, 
     label="100×100 Grid")

# Plot ϕ(x,W/2) for Jacobi
plt2_jacobi = plot(x_values1, ϕ_exact1[mid_index1, :], 
           xlabel="x", ylabel="ϕ(x,W/2)", 
           title="Solution at y=W/2 (Jacobi Method)",
           label="Exact Solution", legend=:topright)
plot!(plt2_jacobi, x_values1, ϕ1_jacobi[mid_index1, :], 
     linestyle=:dash, label="Numerical (50×50)")
plot!(plt2_jacobi, x_values2, ϕ2_jacobi[mid_index2, :], 
     linestyle=:dot, label="Numerical (100×100)")

# Plot the max residual vs iteration number for Gauss-Seidel
plt1_gauss = plot(1:length(residuals1_gauss), residuals1_gauss, 
           xlabel="Iteration", ylabel="Max Residual", 
           title="Convergence: 50×50 vs 100×100 Grid For Gauss-Seidel Method",
           label="50×50 Grid", yscale=:log10, legend=:topright)
plot!(plt1_gauss, 1:length(residuals2_gauss), residuals2_gauss, 
     label="100×100 Grid")

# Plot ϕ(x,W/2) for Gauss-Seidel
plt2_gauss = plot(x_values1, ϕ_exact1[mid_index1, :], 
           xlabel="x", ylabel="ϕ(x,W/2)", 
           title="Solution at y=W/2 (Gauss-Seidel Method)",
           label="Exact Solution", legend=:topright)
plot!(plt2_gauss, x_values1, ϕ1_gauss[mid_index1, :], 
     linestyle=:dash, label="Numerical (50×50)")
plot!(plt2_gauss, x_values2, ϕ2_gauss[mid_index2, :], 
     linestyle=:dot, label="Numerical (100×100)")

# Contour plot of the exact solution
plt3 = contour(x_values2, y_values2, ϕ_exact2, 
              title="Contour Plot of Exact Solution",
              xlabel="x", ylabel="y",
              color=:turbo, levels=20)

plot(plt1_jacobi, plt2_jacobi, plt1_gauss, plt2_gauss, plt3, layout=(5,1), size=(800, 1200))

# Option 2: Display Jacobi and Gauss-Seidel plots separately
# plot(plt1_jacobi, plt2_jacobi, plt3, layout=(3,1), size=(800, 900))
# plot(plt1_gauss, plt2_gauss, layout=(2,1), size=(800, 600))









