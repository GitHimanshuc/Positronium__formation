using Plots
using Distributions
function surge_exp(x::Float64,xth::Float64,varp,varlambda::Float64)

    if x < xth
        return 0.0
    else
        return ((exp(1)/varlambda*(x-xth))^varp*exp(-(x-xth)*varp/varlambda))
    end
end

function surge_poly(x::Float64,xth::Float64,varp,varlambda::Float64)

    if x < xth
        return 0.0
    else
        return (varlambda^(varp)*(varp+1)^(varp+1)*(x-xth)/(x-xth+varlambda*varp)^(varp+1))
    end
end





N = 1000
first = 0.0
last= 500.0/2
currene = collect(first:(last-first)/N:last)[2:end]
psarr = zeros(N)
elasarr = zeros(N)
exarr = zeros(N)
ionarr = zeros(N)


@time for i in 1:N
    psarr[i] = surge_exp(currene[i],6.8,0.5,7.0)
end
@time for i in 1:N
    exarr[i] = .5*surge_poly(currene[i],10.0,1.5,12.0)
end
@time for i in 1:N
    ionarr[i] = 0.3*surge_poly(currene[i],13.6,1,30.0)
end
@time for i in 1:N
    elasarr[i] = 0.4*surge_poly(currene[i],-20.0,2,20.0)
end


plot(currene,psarr,label = "Ps",xlabel = "Energy (eV)",ylabel = "Cross-section")
plot!(currene,exarr,label = "Excitation")
plot!(currene,ionarr,label = "Ion")
plot!(currene,elasarr,label = "Elas")
# savefig("C:\\Users\\Himanshu\\Desktop\\Report\\cross_section_new")
