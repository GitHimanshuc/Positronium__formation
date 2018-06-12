using DataFrames
using Distributions
using Plots



using OtherFunctions



rng = MersenneTwister(1234)

function simulate( energy = 5000.0, N = 10000; file::String ="Hydrogen", MMM::Int = 1 , dens::Float64 = 3.5e7, temp::Float64 = 75.0)



    filename = joinpath("C:\\Users\\Himanshu\\Desktop\\Julia\\Base code for positron transport",file*"XS.csv")
    datada = readtable(filename)   #""" The name of the file containing the data, set the working directory accordingly """
    println("------------------------------------------------------------------------")
    #Astronomical data
    println("\n\nThe ISM medium is : ",file)
    println("The density of ISM is : ",dens," Particles per meter cubed")  #particle per m^3           #density of ISM                     #right now it is CNM
    println("The temperature of ISM is : ",temp," Kelvin")   #in kelvin                    #temperature of ISM



    mass_ratio = 1/1835 #a ratio                      #ratio of mass of the positron and molecule/atom


    rng = MersenneTwister(124)   #""" fixing the seed for some unknown reason """

    normalization = datada[2]   #""" Introduced so that overall normalization can be changed for all in one go """
    threshps = datada[3]./normalization   #""" Positronium formation cross section """
    threshdi = datada[4]./normalization   #""" Direct ionization cross section """
    threshex = datada[6]./normalization  # """ Excitation cross section """
    sizeofdata = size(threshps)[1]    #""" Just calculating the size of the above defined arrays for future use/ Required do not change """


    #Every energy data unless stated otherwise is in eV



    psthresh = datada[8][1]   #""" Positronium formation threshold """
    dirthresh = datada[8][2]   #""" Direct ionization threshold """
    exthresh = datada[8][3]   #""" Excitation threshold """



    # Varibles for time and distance measurements

    tempdist = 0.0
    temptime = 0.0
    avgdist = 0.0
    avgtime = 0.0

    #Defining quantities for future use
    v0 = 0.0  #will store Initial velocity of the particle
    tempvar = 0.0 # will store the collision time or collision distance


     # temporary varibles that will be used to take into account the energy loss due to thermalization and due to other events
    dcurrene = 0.0  # Will keep track of change in energy due to non-thermalization processos
    tcurrene = 0.0  # Will keep track of continuous energy change due to thermalization



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


    fraction = 0.0 # Will be used to lineary extrapolate the values of current crosssections



    println("\n\nThe number of particles simulated is:   ",N)
    println("The mean energy is:  ",energy/2,"  eV")





    vm = sqrt(temp*1*8.314*1000/MMM)     #velocity of molecules/atoms of the ISM in meters/second   (Change according to the dimension the simulation is running in)
    em = 0.5*MMM*1.67e-27*vm^2/(1.6e-19)                #some weird energy to be used in the thermalization formula (eV)


    arrcurrene = rand(Normal(energy, psthresh*2), N)    #This array stores the energy distribution / The positron energy will be sampled form this array

    Elas = datada[2] - datada[3]-datada[4]-datada[5]-datada[6]
    thresholdel = mean(Elas)*1e-20   #Averaged elastic scattering cross section in meters square
    varalpha = sqrt(mass_ratio)*(dens*thresholdel*vm)/(1+mass_ratio)^2  #Parameter to be used in the formula for thermalization, refer William C. Sauder

    #arrcurrene = rand(N)*energy



    for i in 1:N   #""" THE loop """

        currene = arrcurrene[i]   #""" A value of energy choosen from the distribution """

        j = sizeofdata   #""" Stores the index in datada[2] of current energy of the positron  """

        v0 = sqrt((2*currene*1.6e-19)/(9.1e-31))  #velocity in meter/second
        varmax_time = 1e9
        vartime_step = 1e5

        vel_time = make_array2(dens , mass_ratio , vm , v0 , Elas , datada[1]; max_time = varmax_time , start_time= 0.0 , interporder = 1 , time_step = vartime_step)

        # Recursively calculating the average
        avgdist = (i-1)*avgdist/i + tempdist/i
        avgtime = (i-1)*avgtime/i + temptime/i
        temptime = 0.0
        tempdist = 0.0


        dcurrene = 0.0
        tcurrene = currene


        varbeta = acoth(sqrt((currene)/(em)))   #Parameter to be used in the formula for thermalization, refer William C. Sauder



        while currene >= psthresh   #""" Simulates life of a particle """


            a = rand()   #""" Random number to be used in simulation """


            while datada[1][j]>currene   #""" Finds the index corresponding to the current energy/////"""
                j=j-1                      #""" Some improvements required """
            end

            # Without extrapolation
            ########################################################################################################
            ########################################################################################################
            thresholdps = threshps[j]
            thresholddi = threshdi[j]   #Storing the relevant cross sections in temporary variables
            thresholdex = threshex[j]


            # With extrapolation
            ########################################################################################################
            ########################################################################################################
"""
            #The simple linear extrapolation used below can break the code if j becomes 1.
            #For the data I am working with this will not happen as the loop ends at much higher index.

            fraction = (currene - datada[1][j-1])/(datada[1][j] - datada[1][j-1])


            thresholdps = threshps[j-1] + (threshps[j]-threshps[j-1])*fraction
            thresholddi = threshdi[j-1] + (threshdi[j]-threshdi[j-1])*fraction #Storing the relevant cross sections in temporary variables
            thresholdex = threshex[j-1] + (threshex[j]-threshex[j-1])*fraction
"""

            # *********** The formula was derived using the assumption that cross section is independent of the velocity, but that is not the case here.



            # The time estimate is used to get distance estimate
            ########################################################################################################
            ########################################################################################################
"""
            v0 = sqrt((2*currene*1.6e-19)/(9.1e-31))  #velocity in meter/second
            tempvar = 1/(dens*datada[2][j]*1e-20*v0)
            temptime =  temptime + tempvar
            tempdist = tempdist + tempvar*v0
            tcurrene = em*(coth(varbeta + varalpha*temptime))^2
"""
            # distance estimate is used to get time estimate


            #Warning:
            #Do not use the same random number in tempvar and for checking the type of interaction
            #else the correlation will give an additional 2% positron formation.
            v0 = sqrt((2*currene*1.6e-19)/(9.1e-31))  #velocity in meter/second

            tempvar = -log(rand())/(datada[2][j]*1.0e-20*dens)  # The distance covered before the next interaction
            tempdist = tempdist + tempvar
            temptime = temptime + tempvar/v0
            #v0 = vel_time[floor(Int,temptime/varmax_time*length(vel_time))+1][1]
            v0 = vel_time(temptime)[1]
            tcurrene = .5 * 9.1e-31 * v0^2/1.6e-19

            #tcurrene = em*(coth(varbeta + varalpha*temptime))^2



            if a < thresholdps
                posfor = posfor + 1
                break
            elseif a < (thresholdps + thresholddi)
                dirioncount = dirioncount + 1
                dcurrene = dcurrene + dirthresh
            elseif a < (thresholdps + thresholddi + thresholdex)
                excount = excount + 1
                dcurrene = dcurrene + exthresh
            elseif (thresholdps + thresholddi + thresholdex) == 0
                othcount = othcount + 1
                break
            else
                collcount = collcount +1
            end



            currene = tcurrene - dcurrene


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
    println("Average distance travelled is:")
    print(avgdist/3.086e16)
    println(" Pc")
    println("Average time travelled for is:")
    print(avgtime/31536000)
    println(" Years")
    println("Average number of others:")
    println((othcount/N),"\n\n")
    (posfor)/N*100
end




# Calling the function

#simulate( energy = 5000.0, N = 10000; file::String ="Hydrogen", MMM::Int = 1 , dens::Float64 = 3.5e7, temp::Float64 = 75.0)

@time simulate(5000.0 , 10000, file ="Helium", MMM = 4 , dens = 3.5e7, temp = 75.0)
