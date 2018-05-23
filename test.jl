

using DataFrames #""" To read data from .csv files """


tic()#""" To calculate timing of the code """

datada = readtable("MolecularHydrogenXS.csv")   #""" The name of the file containing the data, set the working directory accordingly """



rng = MersenneTwister(124)   #""" fixing the seed for some unknown reason """

normalization = datada[2]   #""" Introduced so that overall normalization can be changed for all in one go """   """ This contains the total cross section """
threshps = datada[3]./normalization   #""" Positronium formation cross section """
threshdi = datada[4]./normalization   #""" Direct ionization cross section """
threshex = datada[6]./normalization  # """ Excitation cross section """
sizeofdata = size(thresh)[1]    #""" Just calculating the size of the above defined arrays for future use/ Required do not change """



psthresh = datada[8][1]   #""" Positronium formation threshold """
dirthresh = datada[8][2]   #""" Direct ionization threshold """
exthresh = datada[8][3]   #""" Excitation threshold """


posfor = 0   #""" To count the number of positroniums formed """
collcount = 0   #""" Number of collision, will be divided by the number of particles to get the average number of collisions """
excount = 0   #""" Same as before, counting total number of excitations """
dirioncount = 0   #""" Counting total number of ionizations """
othcount = 0   #""" There were cases where nothing happened, i.e. all the cross sections at that energy ended up being 0, but this will lead to an infinite loop as our medium is infinite for now. These cases were checked for and this variable keeps track of these events """

a = 0.0   #""" Just a variable that will be used to temporarily store a randomly generated number """
energy = 5000-100   #""" Initial energy """
currene = energy   #""" A temporary variable to store the current energy  """
thresholdps = 0.0   #""" Positronium threshold """
thresholdex = 0.0   #""" Excitation threshold """
thresholddi = 0.0   #""" Direct ionization threshold """
N = 4000   #""" Number of particles to be simulated """
i=1   #""" Variable used in for loop/ redundant """


arrcurrene = rand(Normal(energy, sqrt(13.6*2)), N)    #""" This array stores the energy distribution """



for i in 1:N   #""" THE loop """

    currene = arrcurrene[i]   #""" A value of energy choosen from the distribution """

    j = sizeofdata   #""" Stores the index in datada[2] of current energy of the positron  """


    while currene >= psthresh   #""" Simulates life of a particle """


        a = rand()   #""" Random number to be used in simulation """


        while datada[1][j]>currene   #""" Finds the index corresponding to the current energy/////"""
            j=j-1                      #""" Some improvements required """
        end


        thresholdps = threshps[j]
        thresholddi = threshdi[j]   #""" Storing the relevant cross sections in temporary variables """
        thresholdex = threshex[j]




        if a < thresholdps
            posfor = posfor + 1
            break
        elseif a < (thresholdps + thresholddi)
            dirioncount = dirioncount + 1
            currene = currene - dirthresh
        elseif a < (thresholdps + thresholddi + thresholdex)
            excount = excount + 1
            currene = currene - exthresh
        elseif (thresholdps + thresholddi + thresholdex) == 0
            othcount = othcount + 1
            break
        else
            collcount = collcount +1

        end
    end
end


println("\n\n\n\n\n\nThe results are:\n")
println("Postronium formed:")
print((posfor)/N*100)
println("%")
println("Average direct ionization count is:")
println((dirioncount)/N)
println("Average number of excitations:")
println((excount/N))
println("Average number of others:")
println((othcount/N))
println("Avearage number of elastic collisions:")
println((collcount)/N)


toc()
