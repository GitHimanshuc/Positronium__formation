

using DataFrames


tic()

datada = readtable("HeliumXS.csv")



rng = MersenneTwister(124)

normalization = datada[2]
threshps = datada[3]./normalization
threshdi = datada[4]./normalization
threshex = datada[6]./normalization
sizeofdata = size(thresh)[1]



psthresh = datada[8][1]
dirthresh = datada[8][2]
exthresh = datada[8][3]


posfor = 0
collcount = 0
excount = 0
dirioncount = 0
othcount = 0

a = 0.0
energy = 5000-100
currene = energy
thresholdps = 0.0
thresholdex = 0.0
thresholddi = 0.0
N = 1000
i=1


arrcurrene = rand(Normal(energy, sqrt(13.6*2)), N)



for i in 1:N

    currene = arrcurrene[i]

    j = sizeofdata


    while currene >= psthresh


        a = rand()


        while datada[1][j]>currene
            j=j-1
        end


        thresholdps = threshps[j]
        thresholddi = threshdi[j]
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
