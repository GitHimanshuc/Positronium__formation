__precompile__()
module OtherFunctions
export enrgy_loss!
export make_array
export make_array2
export final_scattering_angle!
export surge_exp
export surge_poly
export elastic_cross_section
export bring_down_the_energy!

using DifferentialEquations
using Dierckx
using Distributions
# p[1]:= dens, number density
# p[2]:= elastic crosssections
# p[3]:= mm , ratio of masses
# p[4]:= vm , averaged velocity of the ISM particles

function make_array(dens , mass_ratio , avg_vel_ism , initial_vel , Elas , Ener; max_time::Float64 = 1e8 , start_time::Float64 = 0.0 , interporder::Int = 1 , time_step::Float64 = 1e4)

    u0 = [initial_vel]
    tspan = (0.0,max_time)
    p = Array{Float64}(3)





    temp = sqrt.(Ener.*(1.6e-19*2/(9.1e-31)))
    interp = Spline1D(temp,Elas,k=interporder)

    p[1] = dens
    p[2] = mass_ratio
    p[3] = avg_vel_ism

    function ext(x, inte)
        interp(x)*1e-20
    end


    function par(du,u,p,t)
     du[1] = -p[1]*ext(u[1],interp)*u[1]^2*(p[2]-(p[3]/u[1])^2)/(1+p[2])^2
    end


    prob = ODEProblem(par,u0,tspan,p)
    sol = solve(prob,Tsit5())
    return sol(collect(linspace(0,max_time,floor(Int,max_time/time_step))))


end

#make_array(dens, mass_ratio , avg_vel_ism , initial_vel , Elas , Ener; max_time::Float64 = 1e8 , start_time::Float64 = 0.0)
function make_array2(dens , mass_ratio , avg_vel_ism , initial_vel , Elas , Ener; max_time::Float64 = 1e8 , start_time::Float64 = 0.0 , interporder::Int = 1 , time_step::Float64 = 1e4)

    u0 = [initial_vel]
    tspan = (0.0,max_time)
    p = Array{Float64}(3)





    temp = sqrt.(Ener.*(1.6e-19*2/(9.1e-31)))
    interp = Spline1D(temp,Elas,k=interporder)

    p[1] = dens
    p[2] = mass_ratio
    p[3] = avg_vel_ism

    function ext(x, inte)
        interp(x)*1e-20
    end


    function par(du,u,p,t)
     du[1] = -p[1]*ext(u[1],interp)*u[1]^2*(p[2]-(p[3]/u[1])^2)/(1+p[2])^2
    end


    prob = ODEProblem(par,u0,tspan,p)
    sol = solve(prob,Tsit5())


end

function final_scattering_angle!(sigma ,theta, phi)
    u=sigma[1]
    v=sigma[2]
    w=sigma[3]

    s = sqrt(1-w*w)
    sigma[1] = ((u*w)*cos(phi) - v*sin(phi))*sin(theta)/s + u*cos(theta)
    sigma[2] = ((v*w)*cos(phi) + u*sin(phi))*sin(theta)/s + v*cos(theta)
    sigma[3] =  -s*(sin(theta)*cos(phi)) + w*cos(theta)


end


function surge_exp(x,xth,varlambda)

    if x < xth
        return 0.0
    else
        return (exp(1)/varlambda*(x-xth)*exp(-(x-xth)/varlambda))
    end
end

function surge_poly(x,xth,varlambda)

    if x < xth
        return 0.0
    else
        return (4*varlambda*(x-xth)/(x-xth+varlambda)^2)
    end
end

function elastic_cross_section(x,xth=200.0,below_xth=0.8,fact=25.0,power=0.7)

    if x < xth
        return below_xth
    else
        return (fact/(x^power))
    end
end






function bring_down_the_energy!(currene , dens  , temperature  , sigma  , varposition_vector  , para , energy_shared_with_inonized_elctron ,time_step = 1e5 ; photon_density::Float64 = 0.26 , magnetic_field_strength::Float64 = 6.3 , helium_fraction::Float64 = 0.0)



    elas = para[1]   #Armstrong squared
    Aion = para[2]
    eion = para[3]
    lion = para[4]
    Apsf = para[5]
    epsf = para[6]
    lpsf = para[7]
    Aexh = para[8]
    eexh = para[9]
    lexh = para[10]



    c = 299792458.0
    mₑ = 9.1e-31


    dist = 0.0
    time = 0.0
    direct_ion_cs = 0.0
    elas = 0.0
    total_cross = 0.0
    temp_dist = 0.0
    count_elas = 0
    count_ioni = 0


    varratio = 1835/(1+1835)^2 # Some factor required in the elastic energy loss formula

    mean_energy_ISM = 3/2 * 1.3807e-23 * temperature/1.6e-19            #some energy to be used in the thermalization formula (eV)


    density_hydrogen = dens*(1.0-helium_fraction)
    density_helium = dens*helium_fraction



    while currene > 5e6



        direct_ion_cs = Aion*surge_poly(currene,eion,lion)
        elas = elastic_cross_section(currene , 200.0, 0.8 , 25.0 , 0.7)
        total_cross = elas + direct_ion_cs

        temp_dist = -log(rand())/(total_cross*dens*1e-20)

        γ = currene*(1.6e-19)/((mₑ * c^2))
        println("γ  ",γ)
        β = sqrt(1.0-1.0/γ^2)

        Δt = temp_dist/(β * c)
        println("Time   ",Δt)

        time = time +Δt

        varposition_vector = varposition_vector + sigma * temp_dist

        final_scattering_angle!(sigma , asin((2*rand()-1)), 2*pi*rand()) # Right now we are using the isotropic scattering case.





        ΔE_elas = -2.0*varratio*(1.0 - mean_energy_ISM/currene)*currene  # Energy loss due to elastic collision
        count_elas = count_elas + 1



        number_of_small_steps = Δt/time_step
        left_time = Δt  - floor(number_of_small_steps)*time_step
        Δt = time_step




        while number_of_small_steps >0.0
            γ = currene*(1.6e-19)/((mₑ * c^2))
            β = sqrt(1.0-1.0/γ^2)
            dE_inco = -2.6e-14 * photon_density * (γ * β)^2     #Inverse compton scattering with photons
            dE_sync = -9.9e-16*(magnetic_field_strength * β * γ * 2/3)^2    #Synchrotron radiation
            dE_brem = (-4.1e-10 * density_hydrogen + -1.1e-9 * density_helium ) * γ*1e-6    # Bremsstrahlung radiation with neutral H and He gas


            currene = currene + (dE_brem + dE_sync + dE_inco)*Δt

            number_of_small_steps = number_of_small_steps -1.0

            if number_of_small_steps <  1.0
                Δt = left_time
                println("Sync   ",dE_sync*number_of_small_steps*time_step)
                println("Inco   ",dE_inco*number_of_small_steps*time_step)
                println("brem   ",dE_brem*number_of_small_steps*time_step)

            end

        end
        println("Elas   ",ΔE_elas)
        println("--------------------------------------------------------------------------------------------------------------")


        currene = currene + ΔE_elas



        if rand() < direct_ion_cs/(total_cross)

            currene = currene*(energy_shared_with_inonized_elctron)

            count_ioni = count_ioni +1

        end


    end

        if (currene > 3.5e6)
            currene = currene - 0.511e6
        end

        return currene ,  time ,  count_elas , count_ioni
    end




    function enrgy_loss!(E_init, temp ,sigma )

        vm = rand(Normal(0.0,sqrt(3/2*1.38e-23*temp/1.67e-27)),3)

        mm = 1.67e-27    # Kg
        m = 9.1e-31   # Kg
        M = mm + m
        Eₗ = E_init*1.6e-19  # J
        vi = sqrt(2*Eₗ/m)*sigma  # m/s

        pᶜ = (m*vi + mm*vm)              #Vector

        # DO not change the order of the next three sentences
        β = sqrt(2m)*1/M*dot(pᶜ,sigma)  # m/s
        final_scattering_angle!(sigma , asin((2*rand()-1)), 2*pi*rand()) # Right now we are using the isotropic scattering case.
        α = sqrt(2m)*1/M*dot(pᶜ,sigma)  # m/s

        1/2*( α^2 -2β * sqrt(Eₗ) + 2Eₗ - sqrt(abs(α^4 - 4*α^2 * β * sqrt(Eₗ) + 4*α^2*Eₗ)))/(1.6e-19) # Final energy to be returned in eV
    end



end
