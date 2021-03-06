using Plots
using OtherFunctions




vardensity  = 90
varnumber = 10000

graphdata = Array{Float64}(vardensity,2)
step1 = 4500/Int(vardensity/2)
nextener = 5000.0 - Int(vardensity/2)*step1
step2 = (nextener-50.0)/Int(vardensity/2)


vfile ="Hydrogen"
vMMM = 1

vdens = 3.5e7
vtemp = 75.0



for i in 1:Int(vardensity/2)
    graphdata[i,1]=5000.0 -step1*i
    graphdata[i,2]=simulate(5000.0 -step1*i, varnumber, file = vfile, MMM = vMMM, dens = vdens, temp = vtemp)
end



for i in (Int(vardensity/2)+1):Int(vardensity)
    graphdata[i,1]= nextener -step2*(i-Int(vardensity/2))
    graphdata[i,2]= simulate(nextener -step2*(i-Int(vardensity/2)), varnumber, file = vfile, MMM = vMMM, dens = vdens, temp = vtemp)
end
plot(graphdata[1:Int(vardensity),1],graphdata[1:Int(vardensity),2],size=(1200,700),label = vfile,xlabel="Mean Energy (eV)",ylabel="Percentage Positronium Formed",title="Percentage Positronium formation vs Mean Energy(Normal Distribution with Std. deviation = Mean Energy/5)",shape=:circle)
savefig("VariationHydrogen.png")














#########################################################################################################################################
#########################################################################################################################################
#########################################################################################################################################


using Plots
using OtherFunctions


elas = 1.0   #Armstrong squared
Aion = .33
eion = 15.0
lion = 5.0
Apsf = 1.0
epsf = 10.0
lpsf = 1.0
Aexh = .33
eexh = 8.0
lexh = 1.0

x_max = 3000.0

currene = collect(20:.1:x_max)

thresholdps = Apsf*surge_exp.(currene,epsf,lpsf)
thresholddi = Aion*surge_poly.(currene,eion,lion)   #Storing the relevant cross sections in temporary variables
thresholdex = Aexh*surge_exp.(currene,eexh,lexh)
elass = rand(length(currene))
elass[:] = elas
atotal=rand(length(currene))
atotal[:] = elas + thresholdps + thresholddi + thresholdex

plot(currene,thresholdps,label="PS",xlabel = "Energy(eV)")
plot!(currene,thresholddi,label="ION")
plot!(currene,thresholdex,label="EX")
plot!(currene,elass,label="ELAS")
plot!(currene,atotal,label="TOT")
