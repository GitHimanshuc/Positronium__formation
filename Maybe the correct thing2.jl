using DataFrames
using Distributions

rng = MersenneTwister(1234)

function simulate( energy = 5000.0, N = 10000; file::String ="Hydrogen", MMM::Int = 1 , dens::Float64 = 3.5e7, temp::Float64 = 75.0)



    filename = joinpath("C:\\Users\\Himanshu\\Desktop\\Julia\\Base code for positron transport",file*"XS.csv")
    datada = readtable(filename)   #""" The name of the file containing the data, set the working directory accordingly """
    println("------------------------------------------------------------------------")
    #Astronomical data
    println("\n\nThe ISM medium is : ",file)
    println("The density of ISM is : ",dens," Particles per meter cubed")  #particle per m^3           #density of ISM                     #right now it is CNM
    println("The temperature of ISM is : ",temp," Kelvin")   #in kelvin                    #temperature of ISM



    mm = 1/1835 #a ratio                      #ratio of mass of the positron and molecule/atom


    rng = MersenneTwister(124)   #""" fixing the seed for some unknown reason """

    normalization = datada[2]   #""" Introduced so that overall normalization can be changed for all in one go """   """ This contains the total cross section """
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
    tm = 0.0 # will store the collision time
     # temporary varibles that will be used to take into account the energy loss due to thermalization and due to other events
    dcurrene = 0.0  #
    tcurrene = 0.0



    posfor = 0   #""" To count the number of positroniums formed """
    collcount = 0   #""" Number of collision, will be divided by the number of particles to get the average number of collisions """
    excount = 0   #""" Same as before, counting total number of excitations """
    dirioncount = 0   #""" Counting total number of ionizations """
    othcount = 0   #""" There were cases where nothing happened, i.e. all the cross sections at that energy ended up being 0, but this will lead to an infinite loop as our medium is infinite for now. These cases were checked for and this variable keeps track of these events """

    a = 0.0   #""" Just a variable that will be used to temporarily store a randomly generated number """
    currene = energy   #""" A temporary variable to store the current energy  """
    thresholdps = 0.0   #""" Positronium threshold """
    thresholdex = 0.0   #""" Excitation threshold """
    thresholddi = 0.0   #""" Direct ionization threshold """




    println("\n\nThe number of particles simulated is:   ",N)
    println("The mean energy is:  ",energy/2,"  eV")





    vm = sqrt(temp*1*8.314*1000/MMM)     #velocity of molecules/atoms of the ISM in meters/second
    em = 0.5*MMM*1.67e-27*vm^2/(1.6e-19)                #some weird energy to be used in the thermalization formula eV


    arrcurrene = rand(Normal(energy/2, energy/10), N)    #""" This array stores the energy distribution """
    thresholdel = mean((1 .-(threshex[200:sizeofdata] + threshdi[200:sizeofdata] + threshps[200:sizeofdata])))*1e-20      #Averaged elastic scattering cross section in meters square
    varalpha = sqrt(mm)*(dens*thresholdel*vm)/(1+mm)^2

    #arrcurrene = rand(N)*energy



    for i in 1:N   #""" THE loop """

        currene = arrcurrene[i]   #""" A value of energy choosen from the distribution """

        j = sizeofdata   #""" Stores the index in datada[2] of current energy of the positron  """

        avgdist = (i-1)*avgdist/i + tempdist/i # Recursively calculating the average
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



            thresholdps = threshps[j]
            thresholddi = threshdi[j]   #""" Storing the relevant cross sections in temporary variables """
            thresholdex = threshex[j]



            #  Model 4, at every loop current energy will be decreasd by the amount dictated by the formula.
            # *********** The formula was derived using the assumption that cross section is independent of the velocity, but that is not the case here.
            # *********** Alternate approach could be to use average elatics cross section that will also make the code faster

            #vm = abs(rand(Normal(0,sqrt(8.314*temp*1000/MMM))))
            v0 = sqrt((2*currene*1.6e-19)/(9.1e-31))  #velocity in meter/second
            tm = 1/(dens*datada[2][j]*1e-20*v0)
            temptime =  temptime + tm
            tempdist = tempdist + tm*v0
            tcurrene = em*(coth(varbeta + varalpha*temptime))^2



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

graphdata = Array{Float64}(49,2)

for i in 1:49
    graphdata[i,1]=5000.0 -100.0*i
    graphdata[i,2]=simulate(5000.0 -100.0*i, 1000, file ="Hydrogen", MMM = 1 , dens = 3.5e7, temp = 75.0)
end
