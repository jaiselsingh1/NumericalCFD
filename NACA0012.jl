using LinearAlgebra
using Plots

Δx = Δy = 0.02
L = W = 5.0

ix = Int(2 * (L / Δx)) + 1
jx = Int(2 * (W / Δy)) + 1

h_sq = 1 / (2 * ((1 / Δx^2) + (1 / Δy^2)))


function dNACA_0012(x)
    return 0.12 / 0.2 * (0.2969 * 0.5 * x^(-0.5) - 0.1260 - 2 * 0.3516 * x + 3 * 0.2843 * x^2 - 4 * 0.1015 * x^3)
end

# Residual calculation
function calculate_residual(ϕ, ix, jx, Δx, Δy)
    max_res = 0.0
    for i in 2:ix-1
        for j in 2:jx-1
            res = (ϕ[i+1, j] - 2 * ϕ[i, j] + ϕ[i-1, j]) / (Δx^2) +
                  (ϕ[i, j+1] - 2 * ϕ[i, j] + ϕ[i, j-1]) / (Δy^2)

            max_res = max(max_res, abs(res))
        end
    end
    return max_res
end

# suggested method is SOR
function incompress_flow()
    Φ = zeros(ix, jx)
    ω = 1.8

    # BCs
    Φ[1, :] .= 0.0 # left
    Φ[ix, :] .= 0.0 # right
    Φ[:, jx] .= 0.0 # top
    # bottom BC is floating

    Φ_next = copy(Φ)
    Φ_prev = copy(Φ)
    ϵ = 1e-5

    residuals = Float64[]
    converged = false

    while !converged
        Φ_prev .= Φ
        for i = 2:ix-1 # 1 and ix are set by the BC
            x = -L + (i - 1) * Δx
            for j = 1:jx-1
                if j == 1
                    if x >= 0 && x <= 1 # airfoil update
                        f_prime = dNACA_0012(x)

                        Φ_next[i, j] .= (Φ[i, j+1] / Δy^2 + (Φ[i+1, j] + Φ[i-1, j]) / Δx - f_prime / Δy) / (2 / Δx + 1 / Δy^2)
                    else # x < 0 || x > 1, symmetry plane update
                        Φ_next[i, j] .= ((Φ[i+1, j] + Φ[i-1, j]) / Δx^2 + (2 * Φ[i, j+1] / Δy^2)) / (1 / h_sq)
                        # have to use the j+1 instead of j-1
                        # Φ_next[i, j] .= ((Φ[i+1, j] + Φ[i-1, j]) / Δx^2 + (2 * Φ[i, j-1] / Δy^2)) / (1 / h_sq)
                    end
                else
                    # general update using the laplace equation
                    Φ_next[i, j] .= ω * h_sq * ((1 / Δx^2) * (Φ[i+1, j] + Φ_next[i-1, j]) + (1 / Δy^2) * (Φ[i, j+1] + Φ_next[i, j-1]) + 1 / h_sq * (1 - 1 / ω) * Φ[i, j])
                end
            end
        end
        Φ .= Φ_next
        diff = Φ .- Φ_prev
        change_norm = norm(diff)

        max_res = calculate_residual(Φ, ix, jx, Δx, Δy)
        push!(residuals, max_res)

        if change_norm < ϵ
            converged = true
        end
    end

    return Φ, residuals
end

Φ, residuals = incompress_flow()
plt_residuals = plot(residuals, yscale=:log10,
    xlabel="Iteration", ylabel="Residual",
    title="Convergence of SOR Method",
    legend=false)
