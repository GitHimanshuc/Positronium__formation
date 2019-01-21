using DataFrames
using Distributions
using Plots
using StaticArrays
using LinearAlgebra
using Dates
using DelimitedFiles

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

    return dot(v2,v2)*0.5*m/1.6e-19

end



function simulate( energy = 5000.0, N = 10000, para = rand(10); MMM::Int = 1 , dens::Float64 = 3.5e7, temp::Float64 = 75.0, Q::Float64 = 1.0,elastic_present=true)


    # These variables store various parameters of the surge functions

    #e.g.   Aion -> A_{ion}     ||      eion -> ϵ   ||  lion ->  λ
    # ion -> ionization   ||  psf -> positron formation   || elas -> elastic   || exh -> Excitation

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



    Q = Q*1.0  # Energy share fraction with electron during ionization  Q = 1.0 no energy is taken away by the electron after the ionization




    # Not included as of now because at a time we only work with either high energy or low energy excitations
    # Aexl = 0.0
    # eexl = eexh
    # lexl = 1.0




    # println("------------------------------------------------------------------------")
    # Astronomical data
    # println("The density of ISM is : ",dens," Particles per meter cubed")  #particle per m^3           #density of ISM                     #right now it is CNM
    # println("The temperature of ISM is : ",temp," Kelvin")   #in kelvin                    #temperature of ISM



    mass_ratio = 1/1835      #ratio of mass of the positron and molecule/atom


    # rng = MersenneTwister(124)   #""" fixing the seed for some unknown reason """




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
    collcount = 0   #Number of elastic collisions
    excount = 0   #"""counting total number of excitations """
    dirioncount = 0   #""" Counting total number of ionizations """

    # Not relevant for the modelled data
    #othcount = 0   #""" There were cases where nothing happened, i.e. all the cross sections at that energy ended up being 0, but this will lead to an infinite loop as our medium is infinite for now. These cases were checked for and this variable keeps track of these events """

    a = 0.0   #"""it will be used to temporarily store a randomly generated number """
    currene = energy   #""" A temporary variable to store the current energy  """


    thresholdps = 0.0   #""" Positronium cross_section """
    thresholdex = 0.0   #""" Excitation cross_section """
    thresholddi = 0.0   #""" Direct ionization cross_section """





    # println("\n\nThe number of particles simulated is:   ",N)
    # println("The mean energy is:  ",energy,"  eV")




    arrcurrene = rand(Normal(energy, 1000), N)    #This array_ stores the energy distribution / The positron energy will be sampled form this array_
    arrinitial_sigma = 2*rand(N,3) .-1 # normalization is Required
    # Normalizing the sigma's norm to 1
    normalization  = arrinitial_sigma[:,1].*arrinitial_sigma[:,1]+arrinitial_sigma[:,2].*arrinitial_sigma[:,2]+arrinitial_sigma[:,3].*arrinitial_sigma[:,3]
    arrinitial_sigma = arrinitial_sigma[:,:]./sqrt.(normalization[:]) # somehow it works

    sigma = [0.0,0.0,0.0] # Will store the direction cosine
    velocity_vector = [.0,.0,.0] # Will hold the current velocity.

    # position_vector = zeros(N,3) # Will be used to hold the final 3D position of the positrons before annhilation
    # varposition_vector = zeros(3) # Will hold the current position vector because I am not sure if accesing array_ every time inside the while loop is efficient


    # The following parameters will be used in the calculation of the final scattering angle
    vma = @MVector zeros(3)
    v1 = @MVector zeros(3)
    vc = @MVector zeros(3)
    v2 = @MVector zeros(3)
    g1 = @MVector zeros(3)
    g2 = @MVector zeros(3)

    stopping_energy = temp/1e4 # Approximately the thermal energy of the ISM particles. If energy falls below this it is assumed that the positron is thermalized.
    # stopping_energy = temp/1e4   # Approximately the thermal energy of the ISM particles. If energy falls below this it is assumed that the positron is thermalized.


    if !elastic_present  # If elastic collisions are absent set the thershold energy to be the positron formation threshold otherwise the code may never stop.
        stopping_energy = psthresh
    end

    flag = true
    ps_formed_flag = false
    dir_happened_flag = false
    arr_timelow = zeros(N) # Will store the time at which a positron either formed positronium or fell below the ps formation threshold.
    arr_energylow = zeros(N)  # Will store the energy at which a positron either formed positronium or fell below the ps formation threshold.




    only_non_ps_elascount = 0.0 # Will store the average number of elastic collision a positron goes through before thermalizing.
    only_ps_elascount  = 0.0 # Will store the average number of elastic collision a positron goes through before forming ps
    avgtime_num = 0 # Stores the number of positrons that did not form positronium. Used in calculating thermalization time and average number of elastic collision.
    temp_elascount = 0 # A temporary variable that will store the number of elastic collisions the last positron went through.

    for i in 1:N   #""" THE loop """

        currene = arrcurrene[i]   #""" A value of energy choosen from the distribution """
        sigma = arrinitial_sigma[i,:]  # A value of initial velocity is choosen


        # Recursively calculating the average and storing information
        # position_vector[i,:] = varposition_vector
        # tempdist = sqrt(sum(varposition_vector.*varposition_vector))
        # avgdist = (i-1)*avgdist/i + tempdist/i

        if ps_formed_flag == true
            only_ps_elascount += temp_elascount  # This is not the average, average will be found using combination of number of particles and ps formation %
        end

        flag = true
        ps_formed_flag = false
        dir_happened_flag = false
        temptime = 0.0
        temp_elascount = 0


        while currene >= stopping_energy   #""" Simulates life of a particle """



            dir_happened_flag = false



            thresholdps = Apsf*surge_exp(currene,epsf,pps,lpsf)
            thresholddi = Aion*surge_poly(currene,eion,pion,lion)   #Storing the relevant cross sections in temporary variables
            thresholdex = Aexh*surge_poly(currene,eexh,pex,lexh)
            elas = Aelas*surge_poly(currene,eelas,pelas,lelas)
            total_cross = elas + thresholdps + thresholddi + thresholdex
            a = rand()   #""" Random number to be used in simulation """
            a = a*total_cross # Normalizing



            #Warning:
            #Do not use the same random number in tempvar and for checking the type of interaction else the correlation will give an additional 2% positron formation.
            v0 = sqrt((2*currene*1.6e-19)/(9.1e-31))  #velocity of positron in meter/second

            tempvar = -log(rand())/(total_cross*dens*1e-20)  # The distance covered before the next interaction

            # varposition_vector = varposition_vector+sigma*tempvar     # Final 3D position

            temptime = temptime + tempvar/v0         # Time positron travelled for before the next _collision
            ########################################################################################################################






            if a < thresholdps
                posfor = posfor + 1
                arr_energylow[i] = currene
                arr_timelow[i]  = temptime
                ps_formed_flag = true    # True means that the positron formed and False implies that the positron fell below the stopping_energy = thermalization energy
                break
            elseif a < (thresholdps + thresholddi)
                dirioncount = dirioncount + 1
                qqq = rand()
                ΔE = ((dirthresh/2)^(-1.1)*(1-qqq) + (currene - (dirthresh/2))^(-1.1)*qqq)^(1/(-1.1)) - (dirthresh/2)
                currene = currene - dirthresh - ΔE
                # currene = currene - dirthresh
                dir_happened_flag = true
                if currene  < 0.0
                    println("SHOUT!!!!!!!!!!!!!!!!!!!!!!!!!!!")
                end

            elseif a < (thresholdps + thresholddi + thresholdex)
                excount = excount + 1
                currene = currene - exthresh
            else
                collcount = collcount + 1    # This counts the number of elastic collision
                # There is no energy loss due to elastic collision until elastic_present is set to True
                temp_elascount += 1
            end


            if elastic_present && !dir_happened_flag
                # If elastic_present is set to True then every collision event (ionization/excitations) is also treated as an elastic collision and energy loss takes place. But these collision are not counted as elastic collision/collcount is not incremented.
                currene = final_scattering_energy_and_direction!(currene,sigma,temp,vma,v1,vc,v2,g1,g2) # Right now we are using the isotropic scattering case.
            else
                # If elastic_present is False but we want to track 3D movement then the following line can be uncommented. This will not change the energy but the direction cosine will be updated.
                # final_scattering_energy_and_direction!(currene,sigma,temp,vma,v1,vc,v2,g1,g2)
            end



            # The following segment stores the distribution of positron energy as they fall below ps formation threshold (The same array also stores the energy at which positrons formed positronium (in the if section of ps formation).). A part of this information is used to generate the histograms.
            # This is not storing the final energy of thermalized positrons.
            if flag && (currene < psthresh)
                arr_energylow[i] = currene
                arr_timelow[i]  = temptime  # This is the time it takes for the positron to fall below ps thershold (The same array also stores the time at which the positron formed positronium (in the if section of ps formation).). (is good approximation to the median life time of a positron.)
                flag = false
            end



        end


        # The following segment calculates the average time it takes for a positron to fall below the stopping_energy/thermalization energy, which a least give some idea about the thermalization time.
        # It also counts the number of elastic collisions that a positron undergoes before falling below the stopping_energy/thermalization energy.
        # In calculating the average we are only considering the positrons that did not form positronium. If we average over all the positrons then the results obtained will be very small, because positron spend a considerable time before thermalizing after falling below ps formation threshold.
        # If elastic_present is set to False then we only calculate the number of collisions and time it take to fall below ps threshold, which are not that meaningful.
        if !ps_formed_flag
            avgtime_num += 1         # Number of positrons that did not form positronium
            avgtime = (avgtime_num-1)*avgtime/avgtime_num + temptime/avgtime_num   # Average thermalization time
            only_non_ps_elascount = only_non_ps_elascount*(avgtime_num-1)/avgtime_num + temp_elascount/avgtime_num   # Average number of elastic collision before thermalizing

        end


    end

    # Uncomment the required lines if individual runs are to be monitered
    # println("\n\nThe results are:\n\n")
    # println("Postronium formed:")
    # print((posfor)/N*100,"%\n")
    # println("Average direct ionization count is:")
    # println((dirioncount)/N)
    # println("Average number of excitations:")
    # println((excount/N))
    # println("Avearage number of elastic collisions:")
    # println((collcount)/N)
    # println("Average distance travelled is:")
    # print(avgdist/(3.086e16), " Pc\n")
    # println("Average time travelled for is:")
    # print(avgtime/31536000,"Years\n")
    # println("Average number of others:")
    # println((othcount/N),"\n\n")


    # return (posfor)/N*100,dirioncount/N,only_non_ps_elascount,excount/N, only_ps_elascount/posfor,avgtime,mean(arr_timelow),arr_energylow
    return (posfor)/N*100,dirioncount/N,only_non_ps_elascount,excount/N, only_ps_elascount/posfor,avgtime,mean(arr_timelow),mean(arr_energylow[arr_energylow .>psthresh])
end




# Calling the function



psformation = zeros(1) # To store ps formation percentage
dir_path = "/home/himanshu/Desktop/graphdel/"
q = 1.0
particles = 5000*15
varenergy = 10000
vardensity = 3.5e7
vartemp = 75.0



#para = [Aelas, Aion, eion, lion, Apsf, epsf, lpsf, Aexh, eexh, lexh]
para = [0.40,0.30,13.6,30.0,1.0,6.8,7.0,0.60,2.0,10.0]
varpsthresh = para[6]


varps = 0.0
varion = 0.0
varelas = 0.0
varexc = 0.0
varpselas = 0.0
varthermalization_time = 0.0
vartime_low = 0.0
varenergy_ps = 0.0

@time varps,varion,varelas,varexc,varpselas, varthermalization_time,vartime_low,varenergy_ps=simulate(varenergy , particles, para, MMM = 1 , dens = vardensity, temp = vartemp, Q=q,elastic_present=true)

folder_path = "/home/himanshu/Desktop/julia_data/setB/"*string(Dates.now())
mkdir(folder_path)
file = open(folder_path*"/with_elastic_Q=0.5.txt","a")


scaled_time1 = varthermalization_time*vardensity/1e15
scaled_time2 = vartime_low*vardensity/1e15


writedlm(file,hcat(100-varps,scaled_time2,varion,varpselas,varexc,varenergy_ps,varelas,scaled_time1))

println(100-varps," ",scaled_time2," ",varion," ",varpselas," ",varexc," ",varenergy_ps)




close(file)
