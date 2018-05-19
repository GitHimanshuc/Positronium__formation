rng = MersenneTwister(1234)
randarr = rand(rng,1000)
using DataFrames

datada = readtable("HydrogenXS.csv")

thresh = datada[3]./datada[2]
sizeofdata = size(thresh)[1]


posfor = 0
collcount = 0

a = 0.0
energy = 1000
currene = energy
threshold = 0.0
N = 1000
i=1







for i in 1:N

    currene = energy -13.6+ 2*13.6/N*i
    while currene >= 6.8
        a = rand()


        i=1
        while datada[1][i]<currene
            i=i+1
            if i==382
                break
            end
        end
        threshold = thresh[i]




        if a < threshold
            posfor=posfor +1
            break
        else
            collcount = collcount +1
            currene = currene -13.6
        end
    end
end


println("\n\n\n\n\n\nThe results are:\n")
println("Postronium formed:")
print((posfor)/N*100)
println("%")

println("Avearage number of collisions:")
print((collcount)/N)
