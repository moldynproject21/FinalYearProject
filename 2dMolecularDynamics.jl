### Imports
using Plots
using Dates
using Random
using DataFrames
using CSV

println("Start Time:  ", now())
t1 = time()

### User Variables
NAtomsIsh = 16
L = 4.5
δr = 0.1
Vo = 0
Δt = 1e-4
rCut = 3
iters = 10

InterAtomDist = 1.0

### Heating Parameters
Heatbath = true
InterHeatingIters = 500
HeatFactor = 1.05

PBCBound = 0.9
PBCDist = PBCBound * L

### Seeding random number generator - Comment out if unneeded
#Random.seed!(2142)

### Setting Up Simulation Parameters
SQAtoms = trunc(Int,ceil(sqrt(NAtomsIsh)))
NAtoms = SQAtoms * SQAtoms

xVelocity = zeros(NAtoms)
xForce = zeros(NAtoms,NAtoms)

yVelocity = zeros(NAtoms)
yForce = zeros(NAtoms,NAtoms)

r = zeros(NAtoms,NAtoms)
angles = zeros(NAtoms,NAtoms)

KE = zeros(NAtoms)
PE = zeros(NAtoms)

KineticEnergy = zeros(iters)
PotentialEnergy = zeros(iters)
TotalEnergy = zeros(iters)

atoms = zeros(2,NAtoms)
atoms_prev = zeros(2,NAtoms)
atoms_new = zeros(2,NAtoms)

atomsOrigins = zeros(2,NAtoms)
MSD = zeros(NAtoms)
AvgMSD = zeros(iters)




today = Dates.today()
mylocale = pwd()
hour = Dates.hour(now())
minute = Dates.minute(now())

mkdir("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations")
mkdir("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/Screens")


### Setting Colours for Plot
colours = [:red, :dodgerblue, :magenta3, :yellow, :green2, :darkorange3, :mediumorchid1, :goldenrod1]

### Creating Grid of Atoms
Offset = ((L - ((SQAtoms-1)*(InterAtomDist)))/2)-InterAtomDist

Gx = zeros(SQAtoms)
for i in 1:SQAtoms
    Gx[i] = (i*InterAtomDist)+Offset
end
Gx = repeat(Gx; outer=[SQAtoms])

Gy = (((1:SQAtoms)*InterAtomDist).+Offset)' .* (ones(SQAtoms))
Gy = reshape(Gy, (1, :))


### Setting Initial Positions and Velocities
for i in 1:NAtoms
    xVelocity[i] = 2*(rand()-0.5)*Vo
    yVelocity[i] = 2*(rand()-0.5)*Vo

    atoms[1,i] = Gx[i] + (2*(rand()-0.5)*δr)
    atoms[2,i] = Gy[i] + (2*(rand()-0.5)*δr)

    atoms_prev[1,i] = atoms[1,i] - (xVelocity[i]*Δt)
    atoms_prev[2,i] = atoms[2,i] - (yVelocity[i]*Δt)

    atomsOrigins[1,i] = atoms[1,i]
    atomsOrigins[2,i] = atoms[2,i]
end

Gx = nothing
Gy = nothing

### Choosing Plotting Backend
Plots.pyplot()

##### Force Animation
anim = @animate for count ∈ 1:iters

    ### Calculating Distances
    Threads.@threads for x in 1:NAtoms
        for y in 1:NAtoms
            ####Original Box
            r[x,y] = sqrt(((atoms[1,x]-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2))
            angles[x,y] = atan((atoms[2,x]-atoms[2,y]),(atoms[1,x]-atoms[1,y]))

            ####Checking Periodic Boundary Atoms

            ####Cardinals
            if r[x,y] > (PBCDist) #x + L
                r[x,y] = sqrt((((atoms[1,x]+L)-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2))
                angles[x,y] = atan((atoms[2,x]-atoms[2,y]),((atoms[1,x]+L)-atoms[1,y]))
            end

            if r[x,y] > (PBCDist) #x - L
                r[x,y] = sqrt((((atoms[1,x]-L)-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2))
                angles[x,y] = atan((atoms[2,x]-atoms[2,y]),((atoms[1,x]-L)-atoms[1,y]))
            end

            if r[x,y] > (PBCDist) #y + L
                r[x,y] = sqrt(((atoms[1,x]-atoms[1,y])^2)+(((atoms[2,x]+L)-atoms[2,y])^2))
                angles[x,y] = atan(((atoms[2,x]+L)-atoms[2,y]),(atoms[1,x]-atoms[1,y]))
            end

            if r[x,y] > (PBCDist) #y - L
                r[x,y] = sqrt(((atoms[1,x]-atoms[1,y])^2)+(((atoms[2,x]-L)-atoms[2,y])^2))
                angles[x,y] = atan(((atoms[2,x]-L)-atoms[2,y]),(atoms[1,x]-atoms[1,y]))
            end
            
            #### Corners
            if r[x,y] > (PBCDist) #x + L, y + L
                r[x,y] = sqrt((((atoms[1,x]+L)-atoms[1,y])^2)+(((atoms[2,x]+L)-atoms[2,y])^2))
                angles[x,y] = atan(((atoms[2,x]+L)-atoms[2,y]),((atoms[1,x]+L)-atoms[1,y]))
            end

            if r[x,y] > (PBCDist) #x + L, y - L
                r[x,y] = sqrt((((atoms[1,x]+L)-atoms[1,y])^2)+(((atoms[2,x]-L)-atoms[2,y])^2))
                angles[x,y] = atan(((atoms[2,x]-L)-atoms[2,y]),((atoms[1,x]+L)-atoms[1,y]))
            end

            if r[x,y] > (PBCDist) #x - L, y + L
                r[x,y] = sqrt((((atoms[1,x]-L)-atoms[1,y])^2)+(((atoms[2,x]+L)-atoms[2,y])^2))
                angles[x,y] = atan(((atoms[2,x]+L)-atoms[2,y]),((atoms[1,x]-L)-atoms[1,y]))
            end

            if r[x,y] > (PBCDist) #x - L, y - L
                r[x,y] = sqrt((((atoms[1,x]-L)-atoms[1,y])^2)+(((atoms[2,x]-L)-atoms[2,y])^2))
                angles[x,y] = atan(((atoms[2,x]-L)-atoms[2,y]),((atoms[1,x]-L)-atoms[1,y]))
            end

            if r[x,y] > (PBCDist) # Back to Center
                r[x,y] = sqrt(((atoms[1,x]-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2))
                angles[x,y] = atan((atoms[2,x]-atoms[2,y]),(atoms[1,x]-atoms[1,y]))
            end
            
        end
    end

    ### Calculating Forces
    Threads.@threads for x in 1:NAtoms
        for y in 1:NAtoms
            if x != y && x <= y
                if  abs(r[x,y]) > rCut
                    xForce[x,y] = 0
                    yForce[x,y] = 0
                else
                    xForce[x,y] = ( 24*(((2/((r[x,y])^13)))-((1/((r[x,y])^7)))) ) * cos(angles[x,y])
                    yForce[x,y] = ( 24*(((2/((r[x,y]^13))))-((1/((r[x,y])^7)))) ) * sin(angles[x,y])
                end
                xForce[y,x] = -xForce[x,y]
                yForce[y,x] = -yForce[x,y]
            end
        end
    end


    ### Heating
    if Heatbath == true
        if (count % InterHeatingIters) == 0
            for x in 1:NAtoms      
                atoms_prev[1,x] = atoms[1,x] - (HeatFactor*(atoms[1,x] - atoms_prev[1,x]))
                atoms_prev[2,x] = atoms[2,x] - (HeatFactor*(atoms[2,x] - atoms_prev[2,x]))
            end
        end
    end


    ### Applying Force
    Threads.@threads for x in 1:NAtoms      
        atoms_new[1,x] = (2 * atoms[1,x]) - (atoms_prev[1,x]) + ((sum(xForce[x,:]))*Δt*Δt)
        atoms_new[2,x] = (2 * atoms[2,x]) - (atoms_prev[2,x]) + ((sum(yForce[x,:]))*Δt*Δt)
    end


    ### Periodic Boundary Correction
    for j in 1:NAtoms
        if atoms_new[1,j] < 0
            ###println("X Went Under 0, Iter: ", count, "    Atom: ", j)
            while atoms_new[1,j] < 0  
                atoms_new[1,j] = mod(atoms_new[1,j],L)

                atoms[1,j] = atoms[1,j] + L
                atoms_prev[1,j] = atoms_prev[1,j] + L
                atomsOrigins[1,j] = atomsOrigins[1,j] + L
                
            end
        end

        if atoms_new[1,j] > L
            ###println("X Went Over L, Iter: ", count, "    Atom: ", j)    
            atoms_new[1,j] = atoms_new[1,j] % L

            atoms[1,j] = atoms[1,j] - L
            atoms_prev[1,j] = atoms_prev[1,j] - L
            atomsOrigins[1,j] = atomsOrigins[1,j] - L         
        end

        if atoms_new[2,j] < 0
            ###println("Y Went Under 0, Iter:  ", count, "    Atom: ", j)
            while atoms_new[2,j] < 0   
                atoms_new[2,j] = mod(atoms_new[2,j],L)

                atoms[2,j] = atoms[2,j] + L
                atoms_prev[2,j] = atoms_prev[2,j] + L
                atomsOrigins[2,j] = atomsOrigins[2,j] + L
                
            end
        end

        if atoms_new[2,j] > L
            ###println("Y Went Over L, Iter: ", count, "    Atom: ", j)     
            atoms_new[2,j] = atoms_new[2,j] % L

            atoms[2,j] = atoms[2,j] - L
            atoms_prev[2,j] = atoms_prev[2,j] - L
            atomsOrigins[2,j] = atomsOrigins[2,j] - L
        end
    end
    
    ### Plot Output
    Plots.scatter(atoms_new[1,:],atoms_new[2,:], label = false ,xlim = (0,L),ylim = (0,L),  title = "2-Dimensional Molecular Dynamics", background_color = :white, seriescolor = colours, size=(1280,720))    

    #### Saving a screen shot every n timesteps
    if (count % 10) == 0 || count == 1
        png("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/Screens/$count")
    end
   
    ### MSD Calculations
    for i in 1:NAtoms
        MSD[i] = (((atoms_new[1,i]-atomsOrigins[1,i])^2)+((atoms_new[2,i]-atomsOrigins[2,i])^2))
    end
    AvgMSD[count] = (sum(MSD)/NAtoms)

    ### Velocity Calculations
    for i in 1:NAtoms
        xVelocity[i] = (atoms_new[1,i] - atoms_prev[1,i])/(2*Δt)
        yVelocity[i] = (atoms_new[2,i] - atoms_prev[2,i])/(2*Δt)
    end

    ### Per Atom Energy Calculations
    for i in 1:NAtoms
        KE[i] = 0.5 * ((xVelocity[i]^2) + (yVelocity[i]^2))
        PE[i] = sqrt(((sum(xForce[i,:]))^2)+((sum(yForce[i,:]))^2))
    end

    ### Average Energy Calculations
    KineticEnergy[count] = (sum(KE)/NAtoms)
    PotentialEnergy[count] = (sum(PE)/NAtoms)
    TotalEnergy[count] = KineticEnergy[count] + PotentialEnergy[count]
    
    ### Setting Atom Position for Next Timestep
    for i in 1:NAtoms
        atoms_prev[1,i] = 1*atoms[1,i]
        atoms_prev[2,i] = 1*atoms[2,i]
        atoms[1,i] = 1*atoms_new[1,i]
        atoms[2,i] = 1*atoms_new[2,i]
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

### Saving all Relevant Simulation Data

mp4(anim, "$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/MolDynVideo.mp4", fps = 24)

vels = [xVelocity,yVelocity]
Es = [KineticEnergy,PotentialEnergy,TotalEnergy]
MSDs = [AvgMSD,AvgMSD]

df = DataFrames.DataFrame(vels)
df1 = DataFrames.DataFrame(xForce)
df2 = DataFrames.DataFrame(yForce)
df3 = DataFrames.DataFrame(atoms)
df4 = DataFrames.DataFrame(atoms_prev)
df5 = DataFrames.DataFrame(MSDs)
df6 = DataFrames.DataFrame(Es)

CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/Velocity.csv",df)
CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/XForce.csv",df1)
CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/YForce.csv",df2)
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

plot(KineticEnergy,title = "KineticEnergy", label = false)
png("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/KineticEnergy")
plot(PotentialEnergy,title = "PotentialEnergy", label = false)
png("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/PotentialEnergy")
plot(TotalEnergy,title = "TotalEnergy", label = false)
png("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/TotalEnergy")

plot(AvgMSD,title = "Δr Squared", label = false)
png("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/AvgMSD")

#gif(anim, "$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/MolDynSlow.gif", fps = 24)
