using Distributions
using OtherFunctions




rng = MersenneTwister(1234)


function simulate( energy = 5000.0, N = 10000, para = rand(10); MMM::Int = 1 , dens::Float64 = 3.5e7, temp::Float64 = 75.0, Q::Float64 = 1.0)



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


    Q = Q*1.0  # Energy share fraction with electron during ionization




    println("------------------------------------------------------------------------")
    #Astronomical data
    println("The density of ISM is : ",dens," Particles per meter cubed")  #particle per m^3           #density of ISM                     #right now it is CNM
    println("The temperature of ISM is : ",temp," Kelvin")   #in kelvin                    #temperature of ISM



    const mass_ratio = 1/1835 #a ratio                      #ratio of mass of the positron and molecule/atom


    rng = MersenneTwister(124)   #""" fixing the seed for some unknown reason """







    #Every energy data unless stated otherwise is in eV



    const psthresh = 6.8   #""" Positronium formation threshold """
    const dirthresh = 13.6   #""" Direct ionization threshold """
    const exthresh = 10.2   #""" Excitation threshold """



    # Varibles for time and distance measurements

    const tempdist = 0.0
    const temptime = 0.0
    const avgdist = 0.0
    const avgtime = 0.0

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





    println("\n\nThe number of particles simulated is:   ",N)
    println("The mean energy is:  ",energy,"  eV")





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
    varposition_vector = [.0,.0,.0] # Will hold the current position vector because I am not sure if accesing array every time inside the while loop is efficient

    return_value_of_function = rand(4)
    #arrcurrene = rand(N)*energy



    for i in 1:N   #""" THE loop """

        currene = arrcurrene[i]   #""" A value of energy choosen from the distribution """
        sigma = arrinitial_sigma[i,:]  # A value of initial velocity is choosen
        # Recursively calculating the average and storing information
        # position_vector[i,:] = varposition_vector
        # tempdist = sqrt(sum(varposition_vector.*varposition_vector))
        # avgdist = (i-1)*avgdist/i + tempdist/i
        # avgtime = (i-1)*avgtime/i + temptime/i
        # temptime = 0.0
        # varposition_vector = [.0,.0,.0]


        # return_value_of_function = bring_down_the_energy!(currene , dens  , temp  , sigma  , varposition_vector  , para, Q )

        # currene = return_value_of_function[1]
        # temptime = retu
        # rn_value_of_function[2]
        # dirioncount = dirioncount + return_value_of_function[4]
        # collcount = collcount + return_value_of_function[3]

        while currene >= psthresh/10   #""" Simulates life of a particle """





            # Without extrapolation
            ########################################################################################################
            ########################################################################################################
            thresholdps = Apsf*surge_exp(currene,epsf,lpsf)
            thresholddi = Aion*surge_poly(currene,eion,lion)   #Storing the relevant cross sections in temporary variables
            thresholdex = Aexh*surge_exp(currene,eexh,lexh)
            elas = elastic_cross_section(currene , 200.0, 0.8 , 25.0 , 0.7)
            a = rand()   #""" Random number to be used in simulation """
            total_cross = elas + thresholdps + thresholddi + thresholdex
            a = a*total_cross # Normalizing


            #Warning:
            #Do not use the same random number in tempvar and for checking the type of interaction
            #else the correlation will give an additional 2% positron formation.
            # v0 = sqrt((2*currene*1.6e-19)/(9.1e-31))  #velocity in meter/second
            #
            # tempvar = -log(rand())/(total_cross*dens*1e-20)  # The distance covered before the next interaction

            # varposition_vector = varposition_vector+sigma*tempvar
            ########################################################################################################################
            # temptime = temptime + tempvar/v0

             # final_scattering_angle!(sigma , asin((2*rand()-1)), 2*pi*rand()) # Right now we are using the isotropic scattering case.



             # If not very important keep this elastic energy loss above all losses, otherwise at low energies domain error may occur (IN the function we require velocity of the particle thus we take squareroot of energy and this requires energy to be positive)
             print(currene)
             currene = enrgy_loss!(currene, temp ,sigma )

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
             println("   ",currene)
"""

                if (1-2*varratio*(1 - em/currene))<1
                    println(currene, "  gained   ", -currene*(2*varratio*(1 - em/currene)))
                end
"""






        end
    end


    println("\n\nThe results are:\n\n")
    println("Postronium formed:")
    print((posfor)/N*100)
    println("%")
    println("Average direct ionization count is:")
    println((dirioncount)/N)
    println("Average number of excitations:")
    println((excount/N))
    println("Avearage number of elastic collisions:")
    println((collcount)/N)
    # println("Average distance travelled is:")
    # print(avgdist/(3.086e16))
    # println(" Pc")
    # println("Average time travelled for is:")
    # print(avgtime/31536000)
    # println(" Years")
    # println("Average number of others:")
    # println((othcount/N),"\n\n")
    (posfor)/N*100
end




# Calling the function

#simulate( energy = 5000.0, N = 10000;, MMM::Int = 1 , dens::Float64 = 3.5e7, temp::Float64 = 75.0)



#@time a=simulate(5000.0 , 5000, MMM = 1 , dens = 3.5e7, temp = 75.0)



"""
points = 2000     # Try to limit the number of points to 10000 or so otherwise the ploting may take time
scatter3d(a[1:points,2],a[1:points,1],a[1:points,3])
writedlm("data.txt",a)

"""


#para = [elas, Aion, eion, lion, Apsf, epsf, lpsf, Aexh, eexh, lexh]

paraA8 = [1.0,1/3,15.0,5.0,1.0,10.0,1.0,1/3,8.0,1.0]
paraA10 = [1.0,1/3,15.0,5.0,1.0,10.0,1.0,1/3,10.0,1.0]
paraA12 = [1.0,1/3,15.0,5.0,1.0,10.0,1.0,1/3,13.0,1.0]
paraB6 = [1.0,1/3,15.0,5.0,1.0,10.0,1.0,1/6,2.0,1.0]
paraB3 = [1.0,1/3,15.0,5.0,1.0,10.0,1.0,1/3,2.0,1.0]
paraB1 = [1.0,1/3,15.0,5.0,1.0,10.0,1.0,1.0,2.0,1.0]
paraC6 = [1.0,1/6,10.0,5.0,1.0,10.0,1.0,0.0,8.0,1.0]
paraC3 = [1.0,1/3,10.0,5.0,1.0,10.0,1.0,0.0,8.0,1.0]
paraC1 = [1.0,1.0,10.0,5.0,1.0,10.0,1.0,0.0,8.0,1.0]


psformation = Array{Float64}(9) # To store ps formation percentage

q = 1.0
particles = 5000
varenergy = 1e4
vardensity = 3.5e7
vartemp = 75.0

"""
@time psformation[1]=simulate(varenergy , particles, paraA8, MMM = 1 , dens = vardensity, temp = vartemp, Q=q)
@time psformation[2]=simulate(varenergy , particles, paraA10, MMM = 1 , dens = vardensity, temp = vartemp, Q=q)
@time psformation[3]=simulate(varenergy , particles, paraA12, MMM = 1 , dens = vardensity, temp = vartemp, Q=q)
@time psformation[4]=simulate(varenergy , particles, paraB6, MMM = 1 , dens = vardensity, temp = vartemp, Q=q)
@time psformation[5]=simulate(varenergy , particles, paraB3, MMM = 1 , dens = vardensity, temp = vartemp, Q=q)
@time psformation[6]=simulate(varenergy , particles, paraB1, MMM = 1 , dens = vardensity, temp = vartemp, Q=q)
@time psformation[7]=simulate(varenergy , particles, paraC6, MMM = 1 , dens = vardensity, temp = vartemp, Q=q)
@time psformation[8]=simulate(varenergy , particles, paraC3, MMM = 1 , dens = vardensity, temp = vartemp, Q=q)
@time psformation[9]=simulate(varenergy , particles, paraC1, MMM = 1 , dens = vardensity, temp = vartemp, Q=q)
"""




q = 1.0
particles = 5000
varenergy = 5e3
vardensity = 3.5e7
vartemp = 75.0


para = [1.0, 4.0, 18.0, 18.0, 8.0, 6.0, 10.0, 3.0, 10.0, 20.0]

@time simulate(varenergy , particles, para, MMM = 1 , dens = vardensity, temp = vartemp, Q=q)



#println(psformation)
