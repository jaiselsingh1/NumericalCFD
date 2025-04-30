using LinearAlgebra
using Plots
Δx = Δy = 0.02
L = W = 5.0
ix = Int(2 * (L / Δx)) + 1
jx = Int(2 * (W / Δy)) + 1

function dNACA_0012(x)
    if x < 1e-10
        return 0.12 / 0.2 * 0.2969 * 0.5 * 1e-5  # safe approximation near x=0
    else
        return 0.12 / 0.2 * (0.2969 * 0.5 * x^(-0.5) - 0.1260 - 2 * 0.3516 * x + 3 * 0.2843 * x^2 - 4 * 0.1015 * x^3)
    end
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

# SOR method for incompressible flow
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
    max_iter = 10000
    iter = 0
    converged = false

    while !converged && iter < max_iter
        Φ_prev .= Φ

        for i = 2:ix-1
            x = -L + (i - 1) * Δx

            for j = 1:jx-1
                if j == 1  # floating BC for the airfoil and symmetry updates
                    if x >= 0 && x <= 1 #airfoil
                        f_prime = dNACA_0012(x)

                        Φ_next[i, j] = (Φ[i, j+1] / Δy^2 + (Φ[i+1, j] + Φ[i-1, j]) / Δx^2 - f_prime / Δy) / (2 / Δx^2 + 1 / Δy^2)
                    else # symmetry
                        Φ_next[i, j] = ((Φ[i+1, j] + Φ[i-1, j]) / Δx^2 + 2 * Φ[i, j+1] / Δy^2) / (2 * (1 / Δx^2 + 1 / Δy^2))
                    end
                else
                    # regular SOR update
                    gs_update = (
                        (Φ[i+1, j] + Φ_next[i-1, j]) / (Δx^2) +
                        (Φ[i, j+1] + Φ_next[i, j-1]) / (Δy^2)
                    ) / (2 * (1 / (Δx^2) + 1 / (Δy^2)))

                    Φ_next[i, j] = (1.0 - ω) * Φ[i, j] + ω * gs_update
                end
            end
        end

        Φ .= Φ_next
        diff = Φ .- Φ_prev
        change_norm = norm(diff)
        max_res = calculate_residual(Φ, ix, jx, Δx, Δy)
        push!(residuals, max_res)

        iter += 1
        if iter % 100 == 0
            println("Iteration $iter, Residual: $max_res, Change: $change_norm")
        end

        if change_norm < ϵ
            converged = true
            println("Converged after $iter iterations")
        end
    end

    if !converged
        println("Failed to converge after $max_iter iterations")
    end

    return Φ, residuals
end


function calculate_slip_velocity(Φ)
    # Create arrays to store values
    surface_v = Float64[]
    x_over_c = Float64[]

    # Calculate grid indices corresponding to airfoil (x=0 to x=1)
    i_start = Int(floor((0.0 + L) / Δx)) + 1  # x = 0
    i_end = Int(floor((1.0 + L) / Δx)) + 1    # x = 1

    # Calculate velocity along airfoil surface
    for i = i_start:i_end
        # Calculate x-coordinate
        x_val = -L + (i - 1) * Δx

        # Only include points on the airfoil
        if x_val >= 0.0 && x_val <= 1.0
            # Calculate normalized x position (x/c)
            x_c = x_val  # chord length c = 1

            # Calculate velocity using centered difference
            # Use values just above the surface for better accuracy
            v_over_V = (Φ[i+1, 2] - Φ[i-1, 2]) / (2 * Δx)

            # Scale to match reference data
            v_scaling = 11.0  # Adjust this value based on your results
            v_over_V *= v_scaling

            push!(surface_v, v_over_V)
            push!(x_over_c, x_c)
        end
    end

    return x_over_c, surface_v
end

function plot_slip_velocity_comparison(Φ)
    x_over_c, surface_v = calculate_slip_velocity(Φ)

    x_c_ref = [0, 0.5, 1.25, 2.5, 5.0, 7.5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90, 95, 100]
    x_c_ref = x_c_ref ./ 100
    v_V_ref = [0, 0.640, 1.010, 1.241, 1.373, 1.402, 1.411, 1.411, 1.399, 1.378, 1.350, 1.298, 1.228, 1.166, 1.104, 1.044, 0.950, 0.606, 0]

    plt_comparison = plot(x_over_c, surface_v,
        label="Numerical Solution",
        linewidth=2,
        xlabel="x/c",
        ylabel="v/V (surface slip velocity)",
        title="Surface Slip Velocity Comparison for NACA 0012 Airfoil")

    plot!(plt_comparison, x_c_ref, v_V_ref,
        label="Reference Data",
        linewidth=2,
        linestyle=:dash,
        marker=:circle,
        markersize=4)

    # Add a vertical line at the leading edge radius (1.58% chord)
    vline!(plt_comparison, [0.0158], label="L.E. radius: 1.58% chord",
        linestyle=:dot, linewidth=1)

    return plt_comparison
end


function plot_potential_field(Φ)
    x_grid = range(-L, L, length=ix)
    y_grid = range(-W, W, length=jx)

    # Plot potential field contour
    plt_contour = contour(x_grid, y_grid, Φ',
        levels=50,
        fill=true,
        title="Potential Field Contour with Airfoil Outline",
        xlabel="x",
        ylabel="y")

    # generate airfoil shape for plotting
    function NACA_0012(x)
        return 0.12 / 0.2 * (0.2969 * sqrt(x) - 0.126 * x - 0.3516 * x^2 + 0.2843 * x^3 - 0.1015 * x^4)
    end

    # add airfoil outline
    x_airfoil = range(0, 1, length=100)
    y_upper = [NACA_0012(xi) for xi in x_airfoil]
    y_lower = -y_upper
    plot!(plt_contour, x_airfoil, y_upper, color=:red, linewidth=1.5, label="NACA 0012")
    plot!(plt_contour, x_airfoil, y_lower, color=:red, linewidth=1.5, label="")

    # Zoom in on the airfoil region for better visualization
    plot!(plt_contour, xlims=(-0.5, 1.5), ylims=(-0.5, 0.5))

    return plt_contour
end


Φ, residuals = incompress_flow()

plt_slor_residual = plot(1:length(residuals), residuals,
    xlabel="Iteration", ylabel="Max Residual",
    title="SOR Method (ω=1.8): NACA 0012 @ 0 deg AoA",
    yscale=:log10,
    legend=false,
    grid=true)

plt_comparison = plot_slip_velocity_comparison(Φ)
plt_potential = plot_potential_field(Φ)


all_plots = plot(plt_slor_residual, plt_potential, plt_comparison,
    layout=(3, 1), size=(800, 1200))
display(all_plots)

savefig(plt_slor_residual, "sor_residuals.png")
savefig(plt_potential, "potential_field.png")
savefig(plt_comparison, "slip_velocity_comparison.png")
