using LinearAlgebra 
using Plots 

# computing the wave drag for a parabolic profile 

function wave_drag(M)
    β = sqrt(M^2 - 1)

    #based on the computational domain 
    y_min = 0 
    y_max = 2 
    x_min = 0 
    x_max = 1 
    
    e_m = 0.1 #percent thickness of airfoil 

    jx = 5001 #number of grid points (y-direction)
    dy = (y_max - y_min)/(jx)

    σ = 1 # convergence condition (courant number)
    dx = σ * β * dy # x is "time like" 
    ix = Int(floor((x_max - x_min)/dx)) + 1

    cd_num = 0.0  # coefficient of drag numerical
    cd_exact = 0.0 # exact soluton for the drag coefficient 

    ϕ = zeros(ix+1, jx) # pertubation potential 

    for i in 2:ix
        for j in 2:(jx-1)  # the jx = 1 falls in the boundary condition 
            # update form 
            ϕ[i+1, j] = (((ϕ[i,j+1] - 2*ϕ[i,j] + ϕ[i,j-1])/dy^2)*(dx^2/β^2)) + 2*ϕ[i,j] - ϕ[i-1,j]
        end 

        #add to the numerical cd value 
        x = x_min + (i-1)*dx

        # bottom boundary condition 
        ϕ[i+1,1] = ϕ[i+1,2] - (2*e_m - 4*e_m*x)*dy
        #top boundary condition 
        ϕ[i+1,jx] = ϕ[i+1,jx-1] 

        cd_num += -4*((ϕ[i+1,1] - ϕ[i,1])) * (2*e_m - 4*e_m*x)  # the second part is the tied to the parabolic profile's derivative 

    end 

    cd_exact = (16/3)*(e_m^2/β) #based on definition given to us 


    return (cd_num = cd_num , cd_exact = cd_exact)
end 

mach_numbers = [1.5, 2.0, 2.5, 3.0, 3.5, 4.0]
cd_numerical = []
cd_exact = []
for M in mach_numbers
    result = wave_drag(M)
    push!(cd_numerical, result.cd_num)
    push!(cd_exact, result.cd_exact)
end 

p = plot(mach_numbers, [cd_numerical cd_exact], 
    label=["Numerical" "Exact"],
    marker=[:circle :square],
    linestyle=[:solid :dash],
    linewidth=2,
    xlabel="Mach Number (M₀)",
    ylabel="Wave Drag Coefficient (cₐ)",
    title="Wave Drag vs Mach Number for Parabolic Profile",
    legend=:topright)

display(p)

