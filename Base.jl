rng = MersenneTwister(1234)
randarr = rand(rng,1000)
using DataFrames

iris = readtable("iris.csv")


posfor = 0
collcount = 0
ann = 0
a = 0.0
energy = 10.0
N = 1e6



formation = 0.3
elastic = 0.8





for i = 1:N

    currene = energy
    while currene >= 6.8
        a = rand()

        if a < formation
            posfor=posfor +1
            break
        elseif a < elastic
            collcount = collcount +1
            currene = currene -13.6
        else
            ann = ann +1
            break

        end
    end
end


println("\n\n\n\n\n\nThe results are:\n")
println("Postronium formed:")
print((posfor)/N*100)
println("%")
println("Annihalated:")
print((ann)/N*100)
println("%")
println("Avearage number of collisions:")
print((collcount)/N)
