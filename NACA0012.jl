using LinearAlgebra
using Plots

Δx = Δy = 0.02
L = W = 5.0

ix = 2 * (L / Δx)
jx = 2 * (W / Δy)

function NACA_0012(x)
    return 0.12 / 0.2 * (0.2969 * sqrt(x) - 0.1260 * x - 0.3516 * x^2 + 0.2843 * x^3 - 0.1015 * x^4)
end

# suggested method is SOR
function incompress_flow()
    Φ = zeros(ix, jx)
    # BCs
    Φ[1, :] = 0.0 # left
    Φ[ix, :] = 0.0 # right
    Φ[:, jx] = 0.0 # top
    # bottom BC is floating

    Φ_next = copy(Φ)

    for i = 2:ix-1 # 1 and ix are set by the BC
        x = -L
        for j = 1:jx
            x += Δx
            if j == 1
                if x >= 0 && x <= 1
                    Φ_next
                else
                    x <= 0 && x >= 1
                    Φ_next
                end
            else

                # general update using the laplace equation

            end

        end

    end

end
