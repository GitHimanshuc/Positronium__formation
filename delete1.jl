using Distributions
using Plots


function final_scattering_energy_and_direction!( Ei , sigma,temp )


    vm = rand(Normal(0.0,sqrt(1.38e-23*temp/1.67e-27)),3)

    mm = 1.67e-27    # Kg
    m = 9.1e-31   # Kg
    M = mm + m
    Ei = Ei*1.6e-19
    v1 = sqrt(2*Ei/m)*sigma  # m/s

    vc = (m*v1 + mm*vm)/M

    g1 = v1 - vm

    θ = acos(1 -2*rand())
    ϕ = 2*pi*rand()

    g2 = norm(g1)*[sin(θ)*cos(ϕ),sin(ϕ)*sin(θ),cos(θ)]

    v2 = (vc + mm/M*g2)
    sigma = v2/norm(v2)

    return dot(v2,v2)*.5*m/1.6e-19


end

function enrgy_loss!(E_init, temp ,sigma )

    vm = rand(Normal(0.0,sqrt(1.38e-23*temp/1.67e-27)),3)

    mm = 1.67e-27    # Kg
    m = 9.1e-31   # Kg
    M = mm + m
    Eₗ = E_init*1.6e-19  # J
    Eₘ = .5*mm*dot(vm,vm)
    vi = sqrt(2*Eₗ/m)*sigma  # m/s

    pᶜ = (m*vi + mm*vm)              #Vector

    # DO not change the order of the next three sentences
    β = sqrt(2m)*1/M*dot(pᶜ,sigma)  # m/s
    final_scattering_angle!(sigma , pi*(2*rand()-1)/2, 2*pi*rand()) # Right now we are using the isotropic scattering case.
    α = sqrt(2m)*1/M*dot(pᶜ,sigma)  # m/s

    γ = dot(vi,vm)



    return 1/2*( α^2 -2β * sqrt(Eₗ) + 2Eₗ + sqrt(abs(α^4 - 4*α^2 * β * sqrt(Eₗ) + 4*α^2*Eₗ)))/(1.6e-19) # Final energy to be returned in eV

    # if γ < sqrt(dot(vi,vi))
    #     return 1/2*( α^2 -2β * sqrt(Eₗ) + 2Eₗ - sqrt(abs(α^4 - 4*α^2 * β * sqrt(Eₗ) + 4*α^2*Eₗ)))/(1.6e-19)
    # else
    #
    # end
end



N= 100000
diff_array = zeros(N)
i = 1
currene = 7.5e3
temp = 7.5e3
sigma = [0.0,1.0,0.0]

while (currene > 0.0) && (i<N)
    diff_array[i] = currene
    i = i+1
    println(currene)

    currene = final_scattering_energy_and_direction!(currene,sigma , temp)

end


plott = plot(diff_array[(diff_array .<1000.0).&(diff_array .!=0.0)],title = "Temperature :  " * string(temp) *" K", xlabel = "collisions", ylabel = "Energy (eV)")




histogram(diff_array[diff_array.!=0],nbins=100)
# savefig(plott,"Energy loss "*string(temp))
