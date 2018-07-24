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
    @. v1 = sqrt(2*Ei/m)*sigma_temp  # m/s

    vc = (m*v1 + mm*vma)/M

    g1 = v1 - vma

    θ = acos(1 -2*rand())
    ϕ = 2*pi*rand()

    g2[1] = norm(g1)*sin(θ)*cos(ϕ)
    g2[2] = norm(g1)*sin(ϕ)*sin(θ)
    g2[3] = norm(g1)*cos(θ)

    v2 = (vc + mm/M*g2)
    sigma_temp = v2/norm(v2)
    # sigma_temp[2] = v2[2]/norm(v2)
    # sigma_temp[3] = v2[3]/norm(v2)

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
    varposition_vector = zeros(3) # Will hold the current position vector because I am not sure if accesing array every time inside the while loop is efficient
    vma =v1=vc=v2=g1=g2= @MVector zeros(3)
    # v1 = @MVector [0.0,0.0,0.0]
    # vc = @MVector [0.0,0.0,0.0]
    # v2 = @MVector [0.0,0.0,0.0]
    # g1 = @MVector [0.0,0.0,0.0]
    # g2 = @MVector [0.0,0.0,0.0]
    #arrcurrene = rand(N)*energy



    for i in 1:N   #""" THE loop """

        currene = arrcurrene[i]   #""" A value of energy choosen from the distribution """
        sigma = @views arrinitial_sigma[i,:]  # A value of initial velocity is choosen



        # Recursively calculating the average and storing information
        # position_vector[i,:] = varposition_vector
        # tempdist = sqrt(sum(varposition_vector.*varposition_vector))
        # avgdist = (i-1)*avgdist/i + tempdist/i
        # avgtime = (i-1)*avgtime/i + temptime/i
        # temptime = 0.0
        # varposition_vector = [.0,.0,.0]


        while currene >= psthresh   #""" Simulates life of a particle """





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
            # v0 = sqrt((2*currene*1.6e-19)/(9.1e-31))  #velocity in meter/second
            #
            # tempvar = -log(rand())/(total_cross*dens*1e-20)  # The distance covered before the next interaction
            #
            # varposition_vector = varposition_vector+sigma*tempvar
            # temptime = temptime + tempvar/v0
            ########################################################################################################################






            if a < thresholdps
                posfor = posfor + 1
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
            end


            # print(currene,"  ->  ")
            if elastic_present
                currene = final_scattering_energy_and_direction!(currene,sigma,temp,vma,v1,vc,v2,g1,g2) # Right now we are using the isotropic scattering case.
            end
            # println(currene)

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
    return (posfor)/N*100,dirioncount/N,collcount/N
end




# Calling the function



psformation = Array{Float64}(1) # To store ps formation percentage

q = 1.0
particles = 50000
varenergy = 1000
vardensity = 3.5e7
vartemp = 75.0



#para = [Aelas, Aion, eion, lion, Apsf, epsf, lpsf, Aexh, eexh, lexh]
para = [1.0,0.30,13.6,30.0,1.0,6.8,7.0,0.5,10.0,12.0]
varnum = 10  #NUmber of points the parmeter should be divided into



array_ps = Array{Float64}(varnum)
array_ion = Array{Float64}(varnum)
array_elas = Array{Float64}(varnum)
array_psfalse = Array{Float64}(varnum)
array_ionfalse = Array{Float64}(varnum)
array_elasfalse = Array{Float64}(varnum)
arrparaA = Array{Float64}(varnum,10)
arrparal = Array{Float64}(varnum,10)
arrparae = Array{Float64}(varnum,10)
paranow  = Array{Float64}(varnum,10)
variation_Q = collect(0:1/varnum:1)[2:end]
Qnow  = variation_Q
variation_paraA = collect(0:1.5/varnum:1.5)[2:end]   # Variation of the amplitude
variation_paral = collect(1:(30-1)/varnum:30)[2:end]   # Variation of the falling rate
variation_parae = collect(8:(40-8)/varnum:40)[2:end]   # Variation of the Excitation threshold


for i in 1:varnum
    arrparaA[i,:] = para
    arrparaA[i,8] = variation_paraA[i]  #Aexh
    arrparal[i,:] = para
    arrparal[i,10] = variation_paral[i]  #lexh
    arrparae[i,:] = para
    arrparae[i,9] = variation_parae[i]  #eexh
end


paranow = arrparal
name = "l"
Qnow[:] = q

for i in 1:varnum
    @time array_ps[i],array_ion[i],array_elas[i]=simulate(varenergy , particles, paranow[i,:], MMM = 1 , dens = vardensity, temp = vartemp, Q=Qnow[i],elastic_present=true)
    @time array_psfalse[i],array_ionfalse[i],array_elasfalse[i]=simulate(varenergy , particles, paranow[i,:], MMM = 1 , dens = vardensity, temp = vartemp, Q=Qnow[i],elastic_present=false)
    println(i/varnum*100,"% done\n")
end
println("--------------------------------------------------------------------------------------------------------")
if paranow == arrparaA
    plot(variation_paraA, array_ps,xlabel = "Relative Magnitude of Excitation (SetB) cross section",ylabel = "Percentage Positronium formed",label = "With recoil")
    plot!(variation_paraA, array_psfalse,shape = :circle,label = "Without recoil")
end
# plot(variation_paraA, array_ion)
# plot(variation_paraA, array_elas)
if paranow == arrparal
    plot(variation_paral, array_ps,xlabel = "Parameter lambda of Excitation (SetB) cross section",ylabel = "Percentage Positronium formed",label = "With recoil")
    plot!(variation_paral, array_psfalse,shape = :circle,label = "Without recoil")
end
# plot(variation_paral, array_ion)
# plot(variation_paral, array_elas)

# plot(variation_Q, array_ps,xlabel = "Q",ylabel = "Percentage Positronium formed")
# plot(variation_Q, array_ion)
# plot(variation_Q, array_elas)
if paranow ==arrparae
    plot(variation_parae, array_ps,xlabel = "Excitation (SetB) threshold",ylabel = "Percentage Positronium formed",label = "With recoil")
    plot!(variation_parae, array_psfalse,shape = :circle,label = "Without recoil")
end
# plot(variation_parae, array_ion)
# plot(variation_parae, array_elas)


savefig("C:\\Users\\Himanshu\\Desktop\\Report\\graphs\\SetB"*name)
