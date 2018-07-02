__precompile__()
module OtherFunctions

export make_array
export make_array2
export final_scattering_angle!
export surge_exp
export surge_poly


using DifferentialEquations
using Dierckx
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




end
