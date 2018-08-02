using Distributions
using Plots
using StaticArrays

function final_scattering_energy_and_direction!( Ei , sigma_temp,tempa ,vma,v1,vc,v2,g1,g2)
    # println("Inside here")

    vma[1] = rand(Normal(0.0,sqrt(1.38e-23*tempa/1.67e-27)))
    vma[2] = rand(Normal(0.0,sqrt(1.38e-23*tempa/1.67e-27)))
    vma[3] = rand(Normal(0.0,sqrt(1.38e-23*tempa/1.67e-27)))

    mm = 1.67e-27    # Kg
    m = 9.1e-31   # Kg
    M = mm + m
    Ei = Ei*1.6e-19
    v1 = sqrt(2*Ei/m)*sigma_temp  # m/s

    vc = (m*v1 + mm*vma)/M

    g1 = v1 - vma

    θ = acos(1 -2*rand())
    ϕ = 2*pi*rand()

    g2[1] = norm(g1)*sin(θ)*cos(ϕ)
    g2[2] = norm(g1)*sin(ϕ)*sin(θ)
    g2[3] = norm(g1)*cos(θ)

    v2 = (vc + mm/M*g2)
    sigma_temp = v2/norm(v2)
    # println(dot(v2,v2)*.5*m/1.6e-19)
    return dot(v2,v2)*.5*m/1.6e-19
end
function   func(N,currene,temp)
    i = 1
    sigma = [0.0,0.0,1.0]
    vma=v1=vc=v2=g1=g2= @MVector [0.0,0.0,0.0]
    while (currene > 0.0) && (i<N)
        diff_array[i] = currene
        i = i+1
        # println(currene)

        currene = final_scattering_energy_and_direction!(currene,sigma , temp,vma,v1,vc,v2,g1,g2)

    end
    return diff_array
end

N= 50000*3
diff_array = zeros(N)
currene = 7.5e3
temp = 7.5e1

func(1,currene,temp)
@time diff_array=func(N,currene,temp)



plott = plot(diff_array[(diff_array .<0.007).&(diff_array .!=0.0)],title = "Temperature :  " * string(75.0) *" K", xlabel = "collisions", ylabel = "Energy (eV)",label = "")




histogram(diff_array[diff_array.!=0],nbins=100)
# savefig(plott,"Energy loss "*string(temp))
