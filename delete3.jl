using Distributions

# const tempstorage3 = @MVector zeros(3)
# function asd()
#     #tempstorage .= zeros(3)
#     tempstorage3 .= @MVector rand(3)
#     # a=rand(3)
#
# end
#
# function surge_exp(x,xth,varlambda)
#
#     if x < xth
#         return 0.0
#     else
#         return (exp(1)/varlambda*(x-xth)*exp(-(x-xth)/varlambda))
#     end
# end
# let
#     vma=v1=vc=v2=g1=g2=zeros(3)
# global function final_scattering_energy_and_direction!( Ei , sigma_temp,tempa )
#
#     #vma=v1=vc=v2=g1=g2=zeros(3)
#
#     vma .= (@SVector randn(3)) * sqrt(1.38e-23*tempa/1.67e-27)
#     #vma .= rand(Normal(0.0,sqrt(1.38e-23*tempa/1.67e-27)),3)
#     #vma .= ntuple
#
#     mm = 1.67e-27    # Kg
#     m = 9.1e-31   # Kg
#     M = mm + m
#     Ei = Ei*1.6e-19
#     v1 .= sqrt(2*Ei/m).*sigma_temp  # m/s
#
#     vc .= (m.*v1 + mm.*vma)./M
#
#     g1 .= v1 .- vma
#
#     θ = acos(1 -2*rand())
#     ϕ = 2*pi*rand()
#
#     g2 .= norm(g1)* @SVector [sin(θ)*cos(ϕ),sin(ϕ)*sin(θ),cos(θ)]
#
#     v2 .= (vc + mm/M*g2)
#     sigma_temp = v2/norm(v2)
#
#     return dot(v2,v2)*.5*m/1.6e-19
#
# end
# end


function final_scattering_energy_and_direction!( Ei , sigma_temp,tempa ,vma,v1,vc,v2,g1,g2)

    # vma=v1=vc=v2=g1=g2=zeros(3)

    vma[1] = rand(Normal(0.0,sqrt(1.38e-23*tempa/1.67e-27)))
    vma[2] = rand(Normal(0.0,sqrt(1.38e-23*tempa/1.67e-27)))
    vma[3] = rand(Normal(0.0,sqrt(1.38e-23*tempa/1.67e-27)))
    # vma .= rand(Normal(0.0,sqrt(1.38e-23*tempa/1.67e-27)),3)
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
    sigma_temp[1] = v2[1]/norm(v2)
    sigma_temp[2] = v2[2]/norm(v2)
    sigma_temp[3] = v2[3]/norm(v2)

    return dot(v2,v2)*.5*m/1.6e-19

end

Ei = 200.0
sigma_temp = [1.0,0.0,0.0]
tempa = 75.0
vma=v1a=vca=v2a=g1a=g2a=  zeros(3)

function qwe(N)
    for i in 1:N
        final_scattering_energy_and_direction!( Ei , sigma_temp,tempa,vma,v1a,vca,v2a,g1a,g2a )
    end
end


qwe(1)
N =100000
@time qwe(N)
# @time qwe(1)
