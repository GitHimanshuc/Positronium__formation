using DataFrames
using Distributions
using Plots
using StaticArrays
# using OtherFunctions

function surge_exp(x::Float64,xth::Float64,varp,varlambda::Float64)

    if x < xth
        return 0.0
    else
        return ((exp(1)/varlambda*(x-xth))^varp*exp(-(x-xth)*varp/varlambda))
    end
end

function surge_poly(x,xth,varp,varlambda)

    if x < xth
        return 0.0
    else
        return (varlambda^(varp)*(varp+1)^(varp+1)*(x-xth)/(x-xth+varlambda*varp)^(varp+1))
    end
end


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

    return dot(v2,v2)*.5*m/1.6e-19

end



function simulate( energy = 5000.0, N = 10000, para = rand(10); MMM::Int = 1 , dens::Float64 = 3.5e7, temp::Float64 = 75.0, Q::Float64 = 1.0,elastic_present=true)



    Aelas = para[1]   #Armstrong squared
    Aion = para[2]
    eion = para[3]
    lion = para[4]
    Apsf = para[5]
    epsf = para[6]
    lpsf = para[7]
    Aexh = para[8]
    eexh = para[9]
    lexh = para[10]

    eelas = -20.0
    lelas = 20.0
    pelas = 2
    pion = 1
    pps = 0.5
    pex = 1.5



    Q = Q*1.0  # Energy share fraction with electron during ionization
    #Not included as of now because at a time we only work with either high energy or low energy excitations
    Aexl = 0.0
    eexl = eexh
    lexl = 1.0




    # println("------------------------------------------------------------------------")
    # #Astronomical data
    # println("The density of ISM is : ",dens," Particles per meter cubed")  #particle per m^3           #density of ISM                     #right now it is CNM
    # println("The temperature of ISM is : ",temp," Kelvin")   #in kelvin                    #temperature of ISM



    mass_ratio = 1/1835 #a ratio                      #ratio of mass of the positron and molecule/atom


    rng = MersenneTwister(124)   #""" fixing the seed for some unknown reason """







    #Every energy data unless stated otherwise is in eV



    psthresh = epsf   #""" Positronium formation threshold """
    dirthresh = eion   #""" Direct ionization threshold """
    exthresh = eexh   #""" Excitation threshold """



    # Varibles for time and distance measurements

    tempdist = 0.0
    temptime = 0.0
    avgdist = 0.0
    avgtime = 0.0

    #Defining quantities for future use
    v0 = 0.0  #will store Initial velocity of the particle
    tempvar = 0.0 # will store the collision time or collision distance





    posfor = 0   #""" To count the number of positroniums formed """
    collcount = 0   #""" Number of collision, will be divided by the number of particles to get the average number of collisions """
    excount = 0   #"""counting total number of excitations """
    dirioncount = 0   #""" Counting total number of ionizations """
    othcount = 0   #""" There were cases where nothing happened, i.e. all the cross sections at that energy ended up being 0, but this will lead to an infinite loop as our medium is infinite for now. These cases were checked for and this variable keeps track of these events """

    a = 0.0   #""" Just a variable that will be used to temporarily store a randomly generated number """
    currene = energy   #""" A temporary variable to store the current energy  """


    thresholdps = 0.0   #""" Positronium threshold """
    thresholdex = 0.0   #""" Excitation threshold """
    thresholddi = 0.0   #""" Direct ionization threshold """





    # println("\n\nThe number of particles simulated is:   ",N)
    # println("The mean energy is:  ",energy,"  eV")





    vm = sqrt(temp*3*8.314*1000/MMM)     #velocity of molecules/atoms of the ISM in meters/second   (Change according to the dimension the simulation is running in)
    em = 0.5*MMM*1.67e-27*vm^2/(1.6e-19)                #some weird energy to be used in the thermalization formula (eV)
    varratio = mass_ratio/(1+mass_ratio)^2



    arrcurrene = rand(Normal(energy, psthresh*2), N)    #This array stores the energy distribution / The positron energy will be sampled form this array
    arrinitial_sigma = 2*rand(N,3)-1 # normalization is Required
    # Normalizing the sigma's norm to 1
    norm  = arrinitial_sigma[:,1].*arrinitial_sigma[:,1]+arrinitial_sigma[:,2].*arrinitial_sigma[:,2]+arrinitial_sigma[:,3].*arrinitial_sigma[:,3]
    arrinitial_sigma = arrinitial_sigma[:,:]./sqrt.(norm[:]) # somehow it works

    sigma = [0.0,0.0,0.0] # Will store the direction cosine
    velocity_vector = [.0,.0,.0] # Will hold the current velocity.
    varalpha = 0.0 # Will store the factor multiplied by which we will get the velocity vector from the direction cosine

    position_vector = Array{Float64}(N,3) # Will be used to hold the final 3D position of the positrons before annhilation
    arrtime_low = Array{Float64}(N)
    varposition_vector = zeros(3) # Will hold the current position vector because I am not sure if accesing array every time inside the while loop is efficient
    vma = @MVector zeros(3)
    v1 = @MVector zeros(3)
    vc = @MVector zeros(3)
    v2 = @MVector zeros(3)
    g1 = @MVector zeros(3)
    g2 = @MVector zeros(3)
    #arrcurrene = rand(N)*energy

    thermal_temp_energy = temp/1e4   # Approximately the thermal energy.

    stopping_energy = thermal_temp_energy

    if !elastic_present  # If elastic collisions are absent set the thershold to be the positron formation threshold otherwise the code may never stop.
        stopping_energy = psthresh
    end

    flag = true
    ps_formed_flag = false
    arr_timelow = zeros(N)
    arr_energylow = zeros(N)




    avgtime_num = 0
    temp_elascount = 0
    only_non_ps_elascount = 0.0

    for i in 1:N   #""" THE loop """

        currene = arrcurrene[i]   #""" A value of energy choosen from the distribution """
        sigma = arrinitial_sigma[i,:]  # A value of initial velocity is choosen


        # Recursively calculating the average and storing information
        # position_vector[i,:] = varposition_vector
        # tempdist = sqrt(sum(varposition_vector.*varposition_vector))
        # avgdist = (i-1)*avgdist/i + tempdist/i

        flag = true
        ps_formed_flag = false
        temptime = 0.0
        temp_elascount = 0

        while currene >= stopping_energy   #""" Simulates life of a particle """





            # Without extrapolation
            ########################################################################################################
            ########################################################################################################
            thresholdps = Apsf*surge_exp(currene,epsf,pps,lpsf)
            thresholddi = Aion*surge_poly(currene,eion,pion,lion)   #Storing the relevant cross sections in temporary variables
            thresholdex = Aexh*surge_poly(currene,eexh,pex,lexh)
            elas = Aelas*surge_poly(currene,eelas,pelas,lelas)
            a = rand()   #""" Random number to be used in simulation """
            total_cross = elas + thresholdps + thresholddi + thresholdex
            a = a*total_cross # Normalizing


            #Warning:
            #Do not use the same random number in tempvar and for checking the type of interaction
            #else the correlation will give an additional 2% positron formation.
            v0 = sqrt((2*currene*1.6e-19)/(9.1e-31))  #velocity in meter/second

            tempvar = -log(rand())/(total_cross*dens*1e-20)  # The distance covered before the next interaction

            # varposition_vector = varposition_vector+sigma*tempvar
            temptime = temptime + tempvar/v0
            ########################################################################################################################






            if a < thresholdps
                posfor = posfor + 1
                arr_energylow[i] = currene
                arr_timelow[i]  = temptime
                ps_formed_flag = true
                # println(currene)
                break
            elseif a < (thresholdps + thresholddi)
                dirioncount = dirioncount + 1
                currene = currene - dirthresh
                currene = Q*currene           # The fraction of energy shared with the emitted electron
            elseif a < (thresholdps + thresholddi + thresholdex)
                excount = excount + 1
                currene = currene - exthresh
            else
                collcount = collcount + 1
                temp_elascount += 1
            end


            # print(currene,"  ->  ")
            if elastic_present
                currene =final_scattering_energy_and_direction!(currene,sigma,temp,vma,v1,vc,v2,g1,g2) # Right now we are using the isotropic scattering case.
            end


            if flag && (currene < psthresh)
                arr_energylow[i] = currene
                arr_timelow[i]  = temptime
                flag = false
            end







        end

        if !ps_formed_flag
            avgtime_num += 1
            avgtime = (avgtime_num-1)*avgtime/avgtime_num + temptime/avgtime_num
            only_non_ps_elascount = only_non_ps_elascount*(avgtime_num-1)/avgtime_num + temp_elascount/avgtime_num

        end


    end


    # println("\n\nThe results are:\n\n")
    # println("Postronium formed:")
    # print((posfor)/N*100)
    # println("%")
    # println("Average direct ionization count is:")
    # println((dirioncount)/N)
    # # println("Average number of excitations:")
    # # println((excount/N))
    # println("Avearage number of elastic collisions:")
    # println((collcount)/N)
    # println("Average distance travelled is:")
    # print(avgdist/(3.086e16))
    # println(" Pc")
    # println("Average time travelled for is:")
    # print(avgtime/31536000)
    # println(" Years")
    # println("Average number of others:")
    # println((othcount/N),"\n\n")
    # println(arrtime_low[1:10])
    # return (posfor)/N*100,dirioncount/N,collcount/N,mean(arrtime_low)/31536000,avgtime/31536000


    return (posfor)/N*100,dirioncount/N,only_non_ps_elascount,avgtime/31536000,mean(arr_timelow)/31536000,arr_energylow
end




# Calling the function



psformation = Array{Float64}(1) # To store ps formation percentage

q = 1.0
particles = 5000*3
varenergy = 10000
vardensity = 3.5e7
vartemp = 75.0



#para = [Aelas, Aion, eion, lion, Apsf, epsf, lpsf, Aexh, eexh, lexh]
para = [1.0,0.30,13.6,30.0,1.0,6.8,7.0,0.0,10.0,12.0]
varpsthresh = para[6]
varnum = 10  #NUmber of points the parmeter should be divided into



array_ps = Array{Float64}(varnum)
array_ion = Array{Float64}(varnum)
array_elas = Array{Float64}(varnum)
array_timelow = Array{Float64}(varnum)
array_thermalization_time = Array{Float64}(varnum)
array_energylow = Array{Float64}(varnum,particles+1)
array_energylow[varnum,particles + 1] = 0.3145

array_psfalse = Array{Float64}(varnum)
array_ionfalse = Array{Float64}(varnum)
array_elasfalse = Array{Float64}(varnum)
array_timelowfalse = Array{Float64}(varnum)
array_thermalization_timefalse = Array{Float64}(varnum)
array_energylowfalse = Array{Float64}(varnum,particles+1)
array_energylowfalse[varnum,particles + 1] = 0.3145


arrparaA = Array{Float64}(varnum,10)
arrparal = Array{Float64}(varnum,10)
arrparae = Array{Float64}(varnum,10)
paranow  = Array{Float64}(varnum,10)
variation_Q = collect(0:1/varnum:1)[2:end]
Qnow  = variation_Q
variation_paraA = collect(0:1.5/varnum:1.5)[2:end]   # Variation of the amplitude
variation_paral = collect(10:(45-10)/varnum:45)[2:end]   # Variation of the falling rate
variation_parae = collect(1:(40-1)/varnum:40)[2:end]   # Variation of the Excitation threshold


for i in 1:varnum
    arrparaA[i,:] = para
    arrparaA[i,2] = variation_paraA[i]  #Aion
    arrparal[i,:] = para
    arrparal[i,4] = variation_paral[i]  #lion
    arrparae[i,:] = para
    arrparae[i,3] = variation_parae[i]  #lion
end

for k in 1:10


paranow = arrparaA
name = "A"
q = k/20.0 + 0.5
Qnow[:] = q
name  = name*" with Q = "*string(q)

for i in 1:varnum
    @time array_ps[i],array_ion[i],array_elas[i], array_thermalization_time[i],array_timelow[i],array_energylow[i,1:end-1]=simulate(varenergy , particles, paranow[i,:], MMM = 1 , dens = vardensity, temp = vartemp, Q=Qnow[i],elastic_present=true)
    @time array_psfalse[i],array_ionfalse[i],array_elasfalse[i],array_thermalization_timefalse[i],array_timelowfalse[i],array_energylowfalse[i,1:end-1]=simulate(varenergy , particles, paranow[i,:], MMM = 1 , dens = vardensity, temp = vartemp, Q=Qnow[i],elastic_present=false)
    println(i/varnum*100,"% done\n")
end
println("--------------------------------------------------------------------------------------------------------")


varfactor = floor(Int,maximum(array_thermalization_time)/maximum(array_timelow))+1.0
array_timelow = array_timelow*varfactor
array_timelowfalse = array_timelowfalse*varfactor

varfactor2 = floor(Int,maximum(array_elas)/maximum(array_ion))+1.0
array_elas = array_elas/varfactor2
varfactor3 = floor(Int,maximum(array_elasfalse)/maximum(array_ion))+1.0
array_elasfalse = array_elasfalse/varfactor3




if paranow == arrparaA
    plot(variation_paraA, array_ps,xlabel = "Relative Magnitude of ionization cross section",ylabel = "Percentage Positronium formed",label = "With recoil")
    plot!(variation_paraA, array_psfalse,shape = :circle,label = "Without recoil")
    savefig("C:\\Users\\Himanshu\\Desktop\\Report\\graphs\\SetA\\SetA"*name)

    plot(variation_paraA, array_ion,xlabel = "Relative Magnitude of ionization cross section",ylabel = "Number of events",label = "Ionizations With recoil")
    plot!(variation_paraA, array_elas,label = "Elastic With recoil/"*string(varfactor2))
    plot!(variation_paraA, array_ionfalse,label = "Ionizations Without recoil",shape = :circle)
    plot!(variation_paraA, array_elasfalse,label = "Elastic Without recoi/"*string(varfactor3),shape = :circle)
    savefig("C:\\Users\\Himanshu\\Desktop\\Report\\graphs\\SetA\\SetA"*name*"_collsions")

    plot(variation_paraA, array_thermalization_time,xlabel = "Relative Magnitude of ionization cross section",ylabel = "Time (years)",label = "Thermalization time")
    plot!(variation_paraA, array_timelow,label = "Time_low*"*string(varfactor)*" with recoil")
    plot!(variation_paraA, array_timelowfalse,label = "Time_low*"*string(varfactor)*" Without recoil",shape = :circle)
    savefig("C:\\Users\\Himanshu\\Desktop\\Report\\graphs\\SetA\\SetA"*name*"_times")



    varsomething = Integer(length(array_energylow[:,1])/2)
    histogram(array_energylow[varsomething,:],label="with recoil",xlabel = "Energy(eV)",alpha = .5)
    histogram!(array_energylowfalse[varsomething,:],label="without recoil",alpha = .5)
    savefig("C:\\Users\\Himanshu\\Desktop\\Report\\graphs\\SetA\\SetA"*name*"_energy_distribution_on_forming_Ps")


    histogram(array_energylow[varsomething,array_energylow[varsomething,:].<varpsthresh],label="with recoil",xlabel = "Energy(eV)",alpha = .5)
    histogram!(array_energylowfalse[varsomething,array_energylowfalse[varsomething,:].<varpsthresh],label="without recoil",alpha = .5)
    savefig("C:\\Users\\Himanshu\\Desktop\\Report\\graphs\\SetA\\SetA"*name*"_energy_distribution_below_psthresh")

end

if paranow == arrparal
    plot(variation_paral, array_ps,xlabel = "Parameter lambda of ionization cross section",ylabel = "Percentage Positronium formed",label = "With recoil")
    plot!(variation_paral, array_psfalse,shape = :circle,label = "Without recoil")
    savefig("C:\\Users\\Himanshu\\Desktop\\Report\\graphs\\SetA\\SetA"*name)

    plot(variation_paral, array_ion,xlabel = "Parameter lambda of ionization cross section",ylabel = "Number of events",label = "Ionizations With recoil")
    plot!(variation_paral, array_elas,label = "Elastic With recoil/"*string(varfactor2))
    plot!(variation_paral, array_ionfalse,label = "Ionizations Without recoil",shape = :circle)
    plot!(variation_paral, array_elasfalse,label = "Elastic Without recoil/"*string(varfactor3),shape = :circle)
    savefig("C:\\Users\\Himanshu\\Desktop\\Report\\graphs\\SetA\\SetA"*name*"_collsions")

    plot(variation_paral, array_thermalization_time,xlabel = "Parameter lambda of ionization cross section",ylabel = "Time (years)",label = "Thermalization time")
    plot!(variation_paral, array_timelow,label = "Time_low*"*string(varfactor)*" with recoil")
    plot!(variation_paral, array_timelowfalse,label = "Time_low*"*string(varfactor)*" Without recoil",shape = :circle)
    savefig("C:\\Users\\Himanshu\\Desktop\\Report\\graphs\\SetA\\SetA"*name*"_times")

    varsomething = Integer(length(array_energylow[:,1])/2)
    histogram(array_energylow[varsomething,:],label="with recoil",xlabel = "Energy(eV)",alpha = .5)
    histogram!(array_energylowfalse[varsomething,:],label="without recoil",alpha = .5)
    savefig("C:\\Users\\Himanshu\\Desktop\\Report\\graphs\\SetA\\SetA"*name*"_energy_distribution_on_forming_Ps")


    histogram(array_energylow[varsomething,array_energylow[varsomething,:].<varpsthresh],label="with recoil",xlabel = "Energy(eV)",alpha = .5)
    histogram!(array_energylowfalse[varsomething,array_energylowfalse[varsomething,:].<varpsthresh],label="without recoil",alpha = .5)
    savefig("C:\\Users\\Himanshu\\Desktop\\Report\\graphs\\SetA\\SetA"*name*"_energy_distribution_below_psthresh")

end

if paranow ==arrparae
    plot(variation_parae, array_ps,xlabel = "Ionization threshold",ylabel = "Percentage Positronium formed",label = "With recoil")
    plot!(variation_parae, array_psfalse,shape = :circle,label = "Without recoil")
    savefig("C:\\Users\\Himanshu\\Desktop\\Report\\graphs\\SetA\\SetA"*name)

    plot(variation_parae, array_ion,xlabel = "Ionization threshold",ylabel = "Number of events",label = "Ionizations With recoil")
    plot!(variation_parae, array_elas,label = "Elastic With recoil/"*string(varfactor2))
    plot!(variation_parae, array_ionfalse,label = "Ionizations Without recoil",shape = :circle)
    plot!(variation_parae, array_elasfalse,label = "Elastic Without recoil/"*string(varfactor3),shape = :circle)
    savefig("C:\\Users\\Himanshu\\Desktop\\Report\\graphs\\SetA\\SetA"*name*"_collsions")


    plot(variation_parae, array_thermalization_time,xlabel = "Ionization threshold",ylabel = "Time (years)",label = "Thermalization time")
    plot!(variation_parae, array_timelow,label = "Time_low*"*string(varfactor)*" with recoil")
    plot!(variation_parae, array_timelowfalse,label = "Time_low*"*string(varfactor)*" Without recoil",shape = :circle)
    savefig("C:\\Users\\Himanshu\\Desktop\\Report\\graphs\\SetA\\SetA"*name*"_times")

    varsomething = Integer(length(array_energylow[:,1])/2)
    histogram(array_energylow[varsomething,:],label="with recoil",xlabel = "Energy(eV)",alpha = .5)
    histogram!(array_energylowfalse[varsomething,:],label="without recoil",alpha = .5)
    savefig("C:\\Users\\Himanshu\\Desktop\\Report\\graphs\\SetA\\SetA"*name*"_energy_distribution_on_forming_Ps")


    histogram(array_energylow[varsomething,array_energylow[varsomething,:].<varpsthresh],label="with recoil",xlabel = "Energy(eV)",alpha = .5)
    histogram!(array_energylowfalse[varsomething,array_energylowfalse[varsomething,:].<varpsthresh],label="without recoil",alpha = .5)
    savefig("C:\\Users\\Himanshu\\Desktop\\Report\\graphs\\SetA\\SetA"*name*"_energy_distribution_below_psthresh")




end


end
