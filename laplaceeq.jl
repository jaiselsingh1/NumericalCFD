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
        term = ((-1)^(n+1) + 1) / n 
        term *= sin(n * π * x / L)
        term *= sinh(n * π * y / L) / sinh(n * π * W / L)
        result += term
    end 
    result *= (2.0 / π)  
    return result 
end 

# Residual calculation
function calculate_residual(ϕ, ix, jx, Δx, Δy)
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

function thomas_algorithm(a::Vector{Float64}, b::Vector{Float64}, c::Vector{Float64},
    d::Vector{Float64}, n::Int64)

    x = copy(d)
    c_prime = copy(c)

    # Setting initial elements
    c_prime[1] /= b[1]
    x[1] /= b[1]

    for i = 2:n
    # Scale factor is for c_prime and x
        scale = 1.0 / (b[i] - c_prime[i-1]*a[i])
        c_prime[i] *= scale
        x[i] = (x[i] - a[i] * x[i-1]) * scale
    end

    # Back-substitution
    for i = n-1:-1:1
        x[i] -= (c_prime[i] * x[i+1])
    end

    return x

end 

function jacobi_method(ix, jx, Δx, Δy, ϵ)
    ϕ = zeros(ix+1, jx+1)
    
    # Apply boundary conditions
    ϕ[:, 1] .= 0.0         # ϕ(x,0) = 0
    ϕ[:, jx+1] .= 1.0      # ϕ(x,W) = 1
    ϕ[1, :] .= 0.0         # ϕ(0,y) = 0
    ϕ[ix+1, :] .= 0.0      # ϕ(L,y) = 0
    
    ϕ_prev = copy(ϕ)
    ϕ_next = copy(ϕ)
    
    residuals = Float64[]
    converged = false
    
    while !converged
        ϕ_prev .= ϕ
        for i in 2:ix
            for j in 2:jx
                ϕ_next[i, j] = (
                    (ϕ[i+1, j] + ϕ[i-1, j]) / (Δx^2) + 
                    (ϕ[i, j+1] + ϕ[i, j-1]) / (Δy^2)
                ) / (2 * (1/(Δx^2) + 1/(Δy^2)))
            end
        end

        ϕ .= ϕ_next
        
        diff = ϕ - ϕ_prev
        change_norm = norm(diff)
        
        max_res = calculate_residual(ϕ, ix, jx, Δx, Δy)
        push!(residuals, max_res)
        
        if change_norm < ϵ
            converged = true
        end
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
    
    ϕ_prev = copy(ϕ)
    ϕ_next = copy(ϕ)
    
    residuals = Float64[]
    converged = false
    
    while !converged
        ϕ_prev .= ϕ
        
        for i in 2:ix
            for j in 2:jx
                ϕ_next[i, j] = (
                    (ϕ[i+1, j] + ϕ_next[i-1, j]) / (Δx^2) + 
                    (ϕ[i, j+1] + ϕ_next[i, j-1]) / (Δy^2)
                ) / (2 * (1/(Δx^2) + 1/(Δy^2)))
            end
        end
        
        ϕ .= ϕ_next
        
        diff = ϕ - ϕ_prev
        change_norm = norm(diff)
        
        max_res = calculate_residual(ϕ, ix, jx, Δx, Δy)
        push!(residuals, max_res)
        
        if change_norm < ϵ
            converged = true
        end
    end
    
    return ϕ, residuals
end

function sor_method(ix, jx, Δx, Δy, ϵ, ω=1.8)
    ϕ = zeros(ix+1, jx+1)
    
    # Apply boundary conditions
    ϕ[:, 1] .= 0.0         # ϕ(x,0) = 0
    ϕ[:, jx+1] .= 1.0      # ϕ(x,W) = 1
    ϕ[1, :] .= 0.0         # ϕ(0,y) = 0
    ϕ[ix+1, :] .= 0.0      # ϕ(L,y) = 0
    
    ϕ_prev = copy(ϕ)
    ϕ_next = copy(ϕ)
    
    residuals = Float64[]
    converged = false
    
    while !converged
        ϕ_prev .= ϕ
        
        for i in 2:ix
            for j in 2:jx
                gs_update = (
                    (ϕ[i+1, j] + ϕ_next[i-1, j]) / (Δx^2) + 
                    (ϕ[i, j+1] + ϕ_next[i, j-1]) / (Δy^2)
                ) / (2 * (1/(Δx^2) + 1/(Δy^2)))
                
                ϕ_next[i, j] = (1.0 - ω) * ϕ[i, j] + ω * gs_update   # gs_update is gauss-seidel equivalent 
            end
        end
        
        ϕ .= ϕ_next
        
        diff = ϕ - ϕ_prev
        change_norm = norm(diff)
        
        max_res = calculate_residual(ϕ, ix, jx, Δx, Δy)
        push!(residuals, max_res)
        
        if change_norm < ϵ
            converged = true
        end
    end
    
    return ϕ, residuals
end

function SLOR(ix, jx, Δx, Δy, ϵ, ω=1.8)
    ϕ = zeros(ix+1, jx+1)
    ϕ_next = copy(ϕ)
    
    # Apply boundary conditions
    ϕ[:, 1] .= 0.0          # Bottom boundary
    ϕ[:, jx+1] .= 1.0       # Top boundary
    ϕ[1, :] .= 0.0          # Left boundary
    ϕ[ix+1, :] .= 0.0       # Right boundary
    
    ϕ_next .= ϕ
    
    residuals = Float64[]
    converged = false
    
    while !converged 
        ϕ_prev = copy(ϕ)
        
        for i = 2:ix 
            a = zeros(jx-1)  # Lower diagonal
            b = zeros(jx-1)  # Main diagonal
            c = zeros(jx-1)  # Upper diagonal
            d = zeros(jx-1)  # Right-hand side
            
            for j = 1:(jx-1)
                # Main diagonal coefficient
                b[j] = -2.0 * ((1.0/(Δx^2)) + (1.0/(Δy^2)))
                
                # Lower diagonal (j-1 term)
                if j > 1
                    a[j] = 1.0/(Δy^2)
                end
                
                # Upper diagonal (j+1 term)
                if j < (jx-1)
                    c[j] = 1.0/(Δy^2)
                end
                
                # Right-hand side from x-direction terms
                d[j] = -(1.0/(Δx^2)) * (ϕ[i+1,j+1] + ϕ_next[i-1,j+1])
                
                # Boundary condition adjustments
                if j == 1
                    # Bottom neighbor (j=0) is a boundary point
                    d[j] -= (1.0/(Δy^2)) * ϕ[i,j]
                end
                
                if j == (jx-1)
                    # Top neighbor (j=jx+1) is a boundary point
                    d[j] -= (1.0/(Δy^2)) * ϕ[i,j+2]
                end
            end
            
            # Solve the tridiagonal system
            solution = thomas_algorithm(a, b, c, d, jx-1)
            
            # Apply SOR update
            for j = 1:(jx-1)
                ϕ_next[i,j+1] = ϕ[i,j+1] + ω * (solution[j] - ϕ[i,j+1])
            end
        end
        
        # Update the solution
        ϕ .= ϕ_next
        
        # Check convergence
        diff = ϕ - ϕ_prev
        change_norm = norm(diff)
        
        max_res = calculate_residual(ϕ, ix, jx, Δx, Δy)
        push!(residuals, max_res)
        
        if change_norm < ϵ
            converged = true
        end
    end
    
    return ϕ, residuals
end


# Generate grid data
x_values1 = range(0, L, length=ix1+1)
y_values1 = range(0, W, length=jx1+1)
x_values2 = range(0, L, length=ix2+1)
y_values2 = range(0, W, length=jx2+1)

# Calculate exact solution
ϕ_exact1 = [exact_solution(x, y, L, W, n_terms) for y in y_values1, x in x_values1]
ϕ_exact2 = [exact_solution(x, y, L, W, n_terms) for y in y_values2, x in x_values2]

# Get middle indices
mid_index1 = div(jx1, 2) + 1
mid_index2 = div(jx2, 2) + 1

# a. Exact solution plots
# a.i: Contour plot
plt_contour = contour(x_values2, y_values2, ϕ_exact2, 
                     title="Contour Plot of Exact Solution",
                     xlabel="x", ylabel="y",
                     color=:turbo, levels=20)

# a.ii: Solution at y=W/2
plt_exact_mid = plot(x_values2, ϕ_exact2[mid_index2, :],
                    title="Exact Solution at y=W/2",
                    xlabel="x", ylabel="ϕ(x,W/2)",
                    label="y1")

# b. Jacobi method
ϕ1_jacobi, residuals1_jacobi = jacobi_method(ix1, jx1, Δx1, Δy1, ϵ)
ϕ2_jacobi, residuals2_jacobi = jacobi_method(ix2, jx2, Δx2, Δy2, ϵ)

# b.a: Residual plot
plt_jacobi_residual = plot(1:length(residuals1_jacobi), residuals1_jacobi, 
                          xlabel="Iteration", ylabel="Max Residual", 
                          title="Jacobi Method: Residual vs Iteration",
                          label="50×50 Grid", yscale=:log10, legend=:topright)
plot!(plt_jacobi_residual, 1:length(residuals2_jacobi), residuals2_jacobi, 
     label="100×100 Grid")

# b.b: Solution plot
plt_jacobi_solution = plot(x_values2, ϕ_exact2[mid_index2, :], 
                          xlabel="x", ylabel="ϕ(x,W/2)", 
                          title="Jacobi Method: Solution at y=W/2",
                          label="Exact Solution", legend=:topright)
plot!(plt_jacobi_solution, x_values1, ϕ1_jacobi[:, mid_index1], 
     linestyle=:dash, label="50×50 Grid")
plot!(plt_jacobi_solution, x_values2, ϕ2_jacobi[:, mid_index2], 
     linestyle=:dot, label="100×100 Grid")

# c. Gauss-Seidel method
ϕ1_gauss, residuals1_gauss = gauss_seidel(ix1, jx1, Δx1, Δy1, ϵ)
ϕ2_gauss, residuals2_gauss = gauss_seidel(ix2, jx2, Δx2, Δy2, ϵ)

# c.a: Residual plot
plt_gauss_residual = plot(1:length(residuals1_gauss), residuals1_gauss, 
                         xlabel="Iteration", ylabel="Max Residual", 
                         title="Gauss-Seidel Method: Residual vs Iteration",
                         label="50×50 Grid", yscale=:log10, legend=:topright)
plot!(plt_gauss_residual, 1:length(residuals2_gauss), residuals2_gauss, 
     label="100×100 Grid")

# c.b: Solution plot
plt_gauss_solution = plot(x_values2, ϕ_exact2[mid_index2, :], 
                         xlabel="x", ylabel="ϕ(x,W/2)", 
                         title="Gauss-Seidel Method: Solution at y=W/2",
                         label="Exact Solution", legend=:topright)
plot!(plt_gauss_solution, x_values1, ϕ1_gauss[:, mid_index1], 
     linestyle=:dash, label="50×50 Grid")
plot!(plt_gauss_solution, x_values2, ϕ2_gauss[:, mid_index2], 
     linestyle=:dot, label="100×100 Grid")

# d. SOR method
ω = 1.8
ϕ1_sor, residuals1_sor = sor_method(ix1, jx1, Δx1, Δy1, ϵ, ω)
ϕ2_sor, residuals2_sor = sor_method(ix2, jx2, Δx2, Δy2, ϵ, ω)

# d.a: Residual plot
plt_sor_residual = plot(1:length(residuals1_sor), residuals1_sor, 
                       xlabel="Iteration", ylabel="Max Residual", 
                       title="SOR Method (ω=1.8): Residual vs Iteration",
                       label="50×50 Grid", yscale=:log10, legend=:topright)
plot!(plt_sor_residual, 1:length(residuals2_sor), residuals2_sor, 
     label="100×100 Grid")

# d.b: Solution plot
plt_sor_solution = plot(x_values2, ϕ_exact2[mid_index2, :], 
                       xlabel="x", ylabel="ϕ(x,W/2)", 
                       title="SOR Method (ω=1.8): Solution at y=W/2",
                       label="Exact Solution", legend=:topright)
plot!(plt_sor_solution, x_values1, ϕ1_sor[:, mid_index1], 
     linestyle=:dash, label="50×50 Grid")
plot!(plt_sor_solution, x_values2, ϕ2_sor[:, mid_index2], 
     linestyle=:dot, label="100×100 Grid")

# e. SLOR method
ω = 1.8
ϕ1_slor, residuals1_slor = SLOR(ix1, jx1, Δx1, Δy1, ϵ, ω)
ϕ2_slor, residuals2_slor = SLOR(ix2, jx2, Δx2, Δy2, ϵ, ω)

# e.a: Residual plot
plt_slor_residual = plot(1:length(residuals1_slor), residuals1_slor, 
                        xlabel="Iteration", ylabel="Max Residual", 
                        title="SLOR Method (ω=1.8): Residual vs Iteration",
                        label="50×50 Grid", yscale=:log10, legend=:topright)
plot!(plt_slor_residual, 1:length(residuals2_slor), residuals2_slor, 
      label="100×100 Grid")

# e.b: Solution plot
plt_slor_solution = plot(x_values2, ϕ_exact2[mid_index2, :], 
                        xlabel="x", ylabel="ϕ(x,W/2)", 
                        title="SLOR Method (ω=1.8): Solution at y=W/2",
                        label="Exact Solution", legend=:topright)
plot!(plt_slor_solution, x_values1, ϕ1_slor[:, mid_index1], 
      linestyle=:dash, label="50×50 Grid")
plot!(plt_slor_solution, x_values2, ϕ2_slor[:, mid_index2], 
      linestyle=:dot, label="100×100 Grid")

# Updated combined plot including SLOR method
plot(
    plt_contour, plt_exact_mid,
    plt_jacobi_residual, plt_jacobi_solution,
    plt_gauss_residual, plt_gauss_solution,
    plt_sor_residual, plt_sor_solution,
    plt_slor_residual, plt_slor_solution,
    layout=(5, 2), size=(1200, 2000)
)