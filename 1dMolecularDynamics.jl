### Imports
using Plots
using Dates
using Random

using DataFrames
using CSV

println("Start Time:  ", now())
t1 = time()

### User Variables
NAtoms = 2
L = 5
δr = 0
Vo = 0.0001
Δt = 1e-3
rCut = 3
iters = 100

InterAtomDist = 1.1

### Heating Parameters
Heatbath = false
InterHeatingIters = 100
HeatFactor = 1.5

PBCBound = 0.9
PBCDist = PBCBound * L

### Seeding random number generator - Comment out if unneeded
#Random.seed!(2142)



### Setting Up Simulation Parameters
xAtoms = zeros(NAtoms)
xVelocity = zeros(NAtoms)
xForce = zeros(NAtoms,NAtoms)

yAxis = ones(1,NAtoms)

r = zeros(NAtoms,NAtoms)

Force = zeros(NAtoms,NAtoms)

KE = zeros(NAtoms)
PE = zeros(NAtoms,NAtoms)
KineticEnergy = zeros(iters)
PotentialEnergy = zeros(iters)
TotalEnergy = zeros(iters)

atoms = zeros(1,NAtoms)
atoms_prev = zeros(1,NAtoms)
atoms_new = zeros(1,NAtoms)

atomsOrigins = zeros(1,NAtoms)
MSD = zeros(NAtoms)

AvgMSD = zeros(iters)

today = Dates.today()
mylocale = pwd()
hour = Dates.hour(now())
minute = Dates.minute(now())

mkdir("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations")
mkdir("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/Screens")



### Creating Grid of Atoms
Offset = (((L-((NAtoms-1)*InterAtomDist))/2)) - InterAtomDist

Gx = zeros(NAtoms)
for i in 1:NAtoms
    Gx[i] = (i*InterAtomDist)+Offset
end

### Setting Initial Positions and Velocities
for i in 1:NAtoms
    atoms[1,i] = Gx[i] +(2*(rand()-0.5)*δr)

    xVelocity[i] = 2*(rand()-0.5)*Vo

    atoms_prev[1,i] = atoms[1,i] - (xVelocity[i]*Δt)

    atomsOrigins[1,i] = atoms[1,i]
end

Plots.pyplot()

##### Force Animation
anim = @animate for count ∈ 1:iters

    ### Calculating Distances
    Threads.@threads for x in 1:NAtoms
        for y in 1:NAtoms
            #### Checking Original Boundary
            r[x,y] = atoms[1,x] - atoms[1,y]

            #### Checking Periodic Boundary Atoms
            if r[x,y] > (PBCDist) #x + L
                r[x,y] =  (atoms[1,x]+L) - atoms[1,y]
            end

            if r[x,y] > (PBCDist) #x - L
                r[x,y] =  (atoms[1,x]-L) - atoms[1,y]
            end
         
        end
    end

    ### Calculating Forces
    Threads.@threads for x in 1:NAtoms
        for y in 1:NAtoms
            if x != y && x <= y
                if  abs(r[x,y]) > rCut
                    xForce[x,y] = 0
                else
                    xForce[x,y] = ( 24*(((2/((r[x,y])^13)))-((1/((r[x,y])^7)))) )
                    
                end
                xForce[y,x] = -xForce[x,y]
                PE[x,y] = xForce[x,y]
            end
        end
    end


    ### Heating
    if Heatbath == true
        if (count % InterHeatingIters) == 0
            for x in 1:NAtoms      
                atoms_prev[1,x] = atoms[1,x] - (HeatFactor*(atoms[1,x] - atoms_prev[1,x]))
            end
        end
    end

    ### Applying Force
    Threads.@threads for x in 1:NAtoms      
        atoms_new[1,x] = (2 * atoms[1,x]) - (atoms_prev[1,x]) + ((sum(xForce[x,:]))*Δt*Δt)
    end


    ### Periodic Boundary Correction
    for j in 1:NAtoms
        
        if atoms_new[1,j] < 0
            #println("X Went Under 0, Iter: ", count, "    Atom: ", j)
            while atoms_new[1,j] < 0          
                atoms_new[1,j] = mod(atoms_new[1,j],L)

                atoms[1,j] = atoms[1,j] + L
                atoms_prev[1,j] = atoms_prev[1,j] + L
                atomsOrigins[1,j] = atomsOrigins[1,j] + L
        
            end
            #println("ExitLoop")
        end

        if atoms_new[1,j] > L
            #println("X Went Over L, Iter: ", count, "    Atom: ", j)
            
            atoms_new[1,j] = atoms_new[1,j] % L

            atoms[1,j] = atoms[1,j] - L
            atoms_prev[1,j] = atoms_prev[1,j] - L
            atomsOrigins[1,j] = atomsOrigins[1,j] - L
        end

    end

    #### Plotting Output
    Plots.scatter(atoms_new[1,:], yAxis[1,:], label = "Step: $count" ,xlim = (0,L),ylim = (0,2),  title = "1-Dimensional Molecular Dynamics", background_color = :white, seriescolor = [:green1], size=(1280,720))    
   
    #### Saving a screen shot every n timesteps
    if (count % 10) == 0 || count == 1
        png("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/Screens/$count")
    end

    ### MSD Calculations
    for i in 1:NAtoms
        MSD[i] = (((atoms_new[1,i]-atomsOrigins[1,i])^2))
    end
    AvgMSD[count] = (sum(MSD)/NAtoms)

    ### Velocity Calculations
    for i in 1:NAtoms
        xVelocity[i] = (atoms_new[1,i] - atoms_prev[1,i])/(2*Δt)
    end

    ### Per Atom Energy Calculations
    for i in 1:NAtoms
        KE[i] = 0.5 * (xVelocity[i]^2)
        PE[i] = sum(xForce[i,:])
    end

    ### Average Energy Calculations
    KineticEnergy[count] = (sum(KE)/NAtoms)
    PotentialEnergy[count] = (sum(PE)/NAtoms)
    TotalEnergy[count] = KineticEnergy[count] + PotentialEnergy[count]

    ### Setting Atom Position for Next Timestep
    for i in 1:NAtoms
        atoms_prev[1,i] = 1*atoms[1,i]

        atoms[1,i] = 1*atoms_new[1,i]
    end

    ### Completion Indicator
    if mod(count,iters/10) == 0
        println("Iterations Done: ",count, "    Time:  ", now())
    end



end

t2 = time()

println("Seconds Taken: ", t2-t1)
println("In Minutes: ", (t2-t1)/60)
println("In Hours: ", ((t2-t1)/60)/60)

### Saving All Relevant Simulation Data
mp4(anim, "$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/MolDynVideo.mp4", fps = 24)

vels = [xVelocity,xVelocity]
Es = [KineticEnergy,PotentialEnergy,TotalEnergy]
MSDs = [AvgMSD,AvgMSD]

df = DataFrames.DataFrame(vels)
df1 = DataFrames.DataFrame(xForce)
df3 = DataFrames.DataFrame(atoms)
df4 = DataFrames.DataFrame(atoms_prev)
df5 = DataFrames.DataFrame(MSDs)
df6 = DataFrames.DataFrame(Es)

CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/Velocity.csv",df)
CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/XForce.csv",df1)
CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/Atoms.csv",df3)
CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/AtomsPrev.csv",df4)
CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/MSDs.csv",df5)
CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/Energies.csv",df6)


open("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/Details.txt","w") do io
    println(io,"Number Of Atoms: ", NAtoms)
    println(io,"Size of Box: ", L)
    println(io,"δr: ", δr)
    println(io,"Initial Velocity: ", Vo)
    println(io,"Timestep: ", Δt)
    println(io,"Cutoff Distance: ", rCut)
    println(io,"Iterations: ", iters)
    println(io,"Starting Inter Atom Distance: ", InterAtomDist)
    println(io,"    ")
    if Heatbath == true
        println(io,"Heating Settings", )
        println(io,"Iterations Between Heating: ", InterHeatingIters)
        println(io,"Heating Factor: ", HeatFactor)
    end
    println(io,"    ")
    println(io,"Time Taken", )
    println(io,"Seconds Taken: ", t2-t1)
    println(io,"In Minutes: ", (t2-t1)/60)
    println(io,"In Hours: ", ((t2-t1)/60)/60)
end

plot(KineticEnergy,title = "KineticEnergy", label = false, size = (1920,1080))
png("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/KineticEnergy")
plot(PotentialEnergy,title = "PotentialEnergy", label = false, size = (1920,1080))
png("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/PotentialEnergy")
plot(TotalEnergy,title = "TotalEnergy", label = false, size = (1920,1080))
png("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/TotalEnergy")

#gif(anim, "$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/MolDynSlow.gif", fps = 24)
