graphdata = Array{Float64}(90,2)

for i in 1:45
    graphdata[i,1]=5000.0 -100.0*i
    graphdata[i,2]=simulate(5000.0 -100.0*i, 20000, file ="Hydrogen", MMM = 1 , dens = 3.5e7, temp = 75.0)
end



for i in 46:90
    graphdata[i,1]=500.0 -10.0*(i-45)
    graphdata[i,2]=simulate(500.0 -10.0*(i-45), 20000, file ="Hydrogen", MMM = 1 , dens = 3.5e7, temp = 75.0)
end
plot(graphdata[1:90,1]./2,graphdata[1:90,2],size=(1200,700),label = "Helium",xlabel="Mean Energy (eV)",ylabel="Percentage Positronium Formed",title="Percentage Positronium formation vs Mean Energy(Normal Distribution with Std. deviation = Mean Energy/5)",shape=:circle)
savefig("VariationHydrogen.png")
