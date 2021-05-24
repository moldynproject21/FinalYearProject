### Imports
using Plots
using Dates
using Random
using DataFrames
using CSV


println("Start Time:  ", now())
t1 = time()

### User Variables
NAtomsIsh = 400
L = 80
δr = 0.15
Vo = 0.5
Δt = 1e-3
rCut = 3
iters = 100

InterAtomDist = 1.2
#### Heating Settings
Heatbath = true
InterHeatingIters = 500
HeatFactor = 0.99

PBCBound = 0.75
PBCDist = PBCBound * L

### Ballistic Parameters
Ballistics = true
FireAfter = 50
M = 1
B2 = 4e-5
BM = B2/M
g = 0#9.8
FiringX = 30
FiringY = 0.0
Range = 1

numFiring = 10

#### If triangles = true, atoms start in triangle lattice instead of square
triangles = true

### Seeding random number generator - Comment out if unneeded
Random.seed!(2412)


### Setting Up Simulation Parameters
SQAtoms = trunc(Int,ceil(sqrt(NAtomsIsh)))
NAtoms = SQAtoms * SQAtoms


xVelocity = zeros(NAtoms+numFiring)
xForce = zeros(NAtoms+numFiring,NAtoms+numFiring)
yVelocity = zeros(NAtoms+numFiring)
yForce = zeros(NAtoms+numFiring,NAtoms+numFiring)

r = zeros(NAtoms+numFiring,NAtoms+numFiring)
angles = zeros(NAtoms+numFiring,NAtoms+numFiring)

KE = zeros(NAtoms+numFiring)
PE = zeros(NAtoms+numFiring)
KineticEnergy = zeros(iters)
PotentialEnergy = zeros(iters)
TotalEnergy = zeros(iters)

atoms = zeros(2,NAtoms+numFiring)
atoms_prev = zeros(2,NAtoms+numFiring)
atoms_new = zeros(2,NAtoms+numFiring)

atomsOrigins = zeros(2,NAtoms+numFiring)
MSD = zeros(NAtoms+numFiring)

AvgMSD = zeros(iters)
D = zeros(iters)

withinRange = false

#### Setting Colours for Target and Projectile
Colours = []
for i in 1:(NAtoms+numFiring)
    if i <= NAtoms
        push!(Colours,:red)
    else
        push!(Colours, :deepskyblue)
    end
end



today = Dates.today()
mylocale = pwd()
hour = Dates.hour(now())
minute = Dates.minute(now())




### Creating Grid of Atoms
Offset = ((L - ((SQAtoms-1)*(InterAtomDist)))/2)-InterAtomDist

Gx = zeros(SQAtoms)
for i in 1:SQAtoms
    Gx[i] = (i*InterAtomDist)+Offset
end
Gx = repeat(Gx; outer=[SQAtoms])

Gy = (((1:SQAtoms)*InterAtomDist).+Offset)' .* (ones(SQAtoms))

Gy = reshape(Gy, (1, :))



if triangles == true
    for i in 1:SQAtoms
        if i % 2 == 0
            for j in 1:SQAtoms
                Gx[((i-1)*SQAtoms)+j] = Gx[((i-1)*SQAtoms)+j] + (InterAtomDist/2)
            end   
        end
    end
end



### Setting Initial Positions and Velocities
for i in 1:NAtoms
    xVelocity[i] =2*(rand()-0.5)*Vo
    yVelocity[i] = 2*(rand()-0.5)*Vo

    atoms[1,i] = Gx[i] +(2*(rand()-0.5)*δr)
    atoms[2,i] = Gy[i] +(2*(rand()-0.5)*δr)

    atoms_prev[1,i] = atoms[1,i] - (xVelocity[i]*Δt)
    atoms_prev[2,i] = atoms[2,i] - (yVelocity[i]*Δt)

    atomsOrigins[1,i] = atoms[1,i]
    atomsOrigins[2,i] = atoms[2,i]
end

Gx = nothing
Gy = nothing




mkdir("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations")
mkdir("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/Screens")


#### Manually Placing Projectile Atoms Into Bullet Shape

## Bullet 1
atoms_prev[1,NAtoms+1] = 0.001
atoms[1,NAtoms+1] = 0.001
atoms_new[1,NAtoms+1] = 0.001

atoms_prev[2,NAtoms+1] = L/2 - 1.2
atoms[2,NAtoms+1] = L/2 - 1.2
atoms_new[2,NAtoms+1] = L/2 - 1.2

xVelocity[NAtoms+1] = FiringX
yVelocity[NAtoms+1] = FiringY
## Bullet 2
atoms_prev[1,NAtoms+2] = 0.001
atoms[1,NAtoms+2] = 0.001
atoms_new[1,NAtoms+2] = 0.001

atoms_prev[2,NAtoms+2] = L/2
atoms[2,NAtoms+2] = L/2
atoms_new[2,NAtoms+2] = L/2

xVelocity[NAtoms+2] = FiringX
yVelocity[NAtoms+2] = FiringY
## Bullet 3
atoms_prev[1,NAtoms+3] = 0.001
atoms[1,NAtoms+3] = 0.001
atoms_new[1,NAtoms+3] = 0.001

atoms_prev[2,NAtoms+3] = L/2 + 1.2
atoms[2,NAtoms+3] = L/2 + 1.2
atoms_new[2,NAtoms+3] = L/2 + 1.2

xVelocity[NAtoms+3] = FiringX
yVelocity[NAtoms+3] = FiringY
## Bullet 4
atoms_prev[1,NAtoms+4] = 0.001 + 1.2
atoms[1,NAtoms+4] = 0.001 + 1.2
atoms_new[1,NAtoms+4] = 0.001 + 1.2

atoms_prev[2,NAtoms+4] = L/2
atoms[2,NAtoms+4] = L/2 
atoms_new[2,NAtoms+4] = L/2 

xVelocity[NAtoms+4] = FiringX
yVelocity[NAtoms+4] = FiringY
## Bullet 5
atoms_prev[1,NAtoms+5] = 0.001 - 1.2
atoms[1,NAtoms+5] = 0.001 - 1.2
atoms_new[1,NAtoms+5] = 0.001 - 1.2

atoms_prev[2,NAtoms+5] = L/2
atoms[2,NAtoms+5] = L/2 
atoms_new[2,NAtoms+5] = L/2 

xVelocity[NAtoms+5] = FiringX
yVelocity[NAtoms+5] = FiringY

## Bullet 6
atoms_prev[1,NAtoms+6] = 0.001 - 1.2
atoms[1,NAtoms+6] = 0.001 - 1.2
atoms_new[1,NAtoms+6] = 0.001 - 1.2

atoms_prev[2,NAtoms+6] = L/2 - 1.2
atoms[2,NAtoms+6] = L/2 - 1.2
atoms_new[2,NAtoms+6] = L/2 - 1.2

xVelocity[NAtoms+6] = FiringX
yVelocity[NAtoms+6] = FiringY

## Bullet 7
atoms_prev[1,NAtoms+7] = 0.001 - 1.2
atoms[1,NAtoms+7] = 0.001 - 1.2
atoms_new[1,NAtoms+7] = 0.001 - 1.2

atoms_prev[2,NAtoms+7] = L/2 + 1.2
atoms[2,NAtoms+7] = L/2 + 1.2
atoms_new[2,NAtoms+7] = L/2 + 1.2

xVelocity[NAtoms+7] = FiringX
yVelocity[NAtoms+7] = FiringY

## Bullet 8
atoms_prev[1,NAtoms+8] = 0.001 + 1.2
atoms[1,NAtoms+8] = 0.001 + 1.2
atoms_new[1,NAtoms+8] = 0.001 + 1.2

atoms_prev[2,NAtoms+8] = L/2 + 1.2
atoms[2,NAtoms+8] = L/2 + 1.2
atoms_new[2,NAtoms+8] = L/2 + 1.2

xVelocity[NAtoms+8] = FiringX
yVelocity[NAtoms+8] = FiringY

## Bullet 9
atoms_prev[1,NAtoms+9] = 0.001 + 1.2
atoms[1,NAtoms+9] = 0.001 + 1.2
atoms_new[1,NAtoms+9] = 0.001 + 1.2

atoms_prev[2,NAtoms+9] = L/2 - 1.2
atoms[2,NAtoms+9] = L/2 - 1.2
atoms_new[2,NAtoms+9] = L/2 - 1.2

xVelocity[NAtoms+9] = FiringX
yVelocity[NAtoms+9] = FiringY

## Bullet 9
atoms_prev[1,NAtoms+10] = 0.001 + 2.4
atoms[1,NAtoms+10] = 0.001 + 2.4
atoms_new[1,NAtoms+10] = 0.001 + 2.4

atoms_prev[2,NAtoms+10] = L/2 
atoms[2,NAtoms+10] = L/2 
atoms_new[2,NAtoms+10] = L/2 

xVelocity[NAtoms+10] = FiringX
yVelocity[NAtoms+10] = FiringY

Plots.pyplot()

##### Force Animation
anim = @animate for count ∈ 1:iters

    ### Calculating Distances
    Threads.@threads for x in 1:NAtoms+numFiring
        for y in 1:NAtoms+numFiring
            ####Original Box
            r[x,y] = sqrt(((atoms[1,x]-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2))
            angles[x,y] = atan((atoms[2,x]-atoms[2,y]),(atoms[1,x]-atoms[1,y]))

            #### Checking Periodic Boundary Atoms
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

    
    ### Swapping From Ballistic to LJ
    if withinRange == false && minimum(r[(NAtoms+1:NAtoms+numFiring),(1:NAtoms)]) < Range
        global withinRange = true
        global Δt = 1e-5
        println("Switch to LJ in Iter:  ",  count)
        #global storer[:] .= (r[126,:])
        global Ballistics = false
    end

    ### Ballistic Calculations
    if Ballistics == true && withinRange == false

        if (count > FireAfter) 
            for numBull in 1:numFiring
                atoms_new[1,NAtoms+numBull] = atoms[1,NAtoms+numBull] + (xVelocity[NAtoms+numBull]*Δt)
                atoms_new[2,NAtoms+numBull] = atoms[2,NAtoms+numBull] + (yVelocity[NAtoms+numBull]*Δt)

                if atoms_new[2,NAtoms+numBull] < 0
                    atoms_new[2,NAtoms+numBull] = 0
                end

                v = sqrt((xVelocity[NAtoms+numBull]^2)+(yVelocity[NAtoms+numBull]^2))
        
                XDrag = - (BM)*(v)*(xVelocity[NAtoms+numBull])
                YDrag = - (BM)*(v)*(yVelocity[NAtoms+numBull])

                xVelocity[NAtoms+numBull] = xVelocity[NAtoms+numBull] + ((XDrag/M)*(Δt))
                yVelocity[NAtoms+numBull] = yVelocity[NAtoms+numBull] + ((-g + (YDrag/M))*(Δt))
            end
        end
        
    end






    ### Calculating Forces
    Threads.@threads for x in 1:NAtoms+numFiring
        for y in 1:NAtoms+numFiring
            if x != y
                if  abs(r[x,y]) > rCut
                    xForce[x,y] = 0
                    yForce[x,y] = 0
                else
                    xForce[x,y] = ( 24*(((2/((r[x,y])^13)))-((1/((r[x,y])^7)))) ) * cos(angles[x,y])
                    yForce[x,y] = ( 24*(((2/((r[x,y]^13))))-((1/((r[x,y])^7)))) ) * sin(angles[x,y])
                end
            end
        end
    end


    ### Heating
    if Heatbath == true && withinRange == false
        if (count % InterHeatingIters) == 0
            for x in 1:NAtoms      
                atoms_prev[1,x] = atoms[1,x] - (HeatFactor*(atoms[1,x] - atoms_prev[1,x]))
                atoms_prev[2,x] = atoms[2,x] - (HeatFactor*(atoms[2,x] - atoms_prev[2,x]))
            end
        end
    end




    ### Applying Force
    if withinRange == false
        for x in 1:NAtoms      
            atoms_new[1,x] = (2 * atoms[1,x]) - (atoms_prev[1,x]) + ((sum(xForce[x,:]))*Δt*Δt)
            atoms_new[2,x] = (2 * atoms[2,x]) - (atoms_prev[2,x]) + ((sum(yForce[x,:]))*Δt*Δt)
        end
    else
        for x in 1:NAtoms+numFiring      
            atoms_new[1,x] = (2 * atoms[1,x]) - (atoms_prev[1,x]) + ((sum(xForce[x,:]))*Δt*Δt)
            atoms_new[2,x] = (2 * atoms[2,x]) - (atoms_prev[2,x]) + ((sum(yForce[x,:]))*Δt*Δt)
        end
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
        end

        if atoms_new[1,j] > L
            #println("X Went Over L, Iter: ", count, "    Atom: ", j)   
            atoms_new[1,j] = atoms_new[1,j] % L

            atoms[1,j] = atoms[1,j] - L
            atoms_prev[1,j] = atoms_prev[1,j] - L
            atomsOrigins[1,j] = atomsOrigins[1,j] - L
        end

        if atoms_new[2,j] < 0
            #println("Y Went Under 0, Iter:  ", count, "    Atom: ", j)
            while atoms_new[2,j] < 0
                atoms_new[2,j] = mod(atoms_new[2,j],L)

                atoms[2,j] = atoms[2,j] + L
                atoms_prev[2,j] = atoms_prev[2,j] + L
                atomsOrigins[2,j] = atomsOrigins[2,j] + L
            end
        end

        if atoms_new[2,j] > L
            #println("Y Went Over L, Iter: ", count, "    Atom: ", j)  
            atoms_new[2,j] = atoms_new[2,j] % L

            atoms[2,j] = atoms[2,j] - L
            atoms_prev[2,j] = atoms_prev[2,j] - L
            atomsOrigins[2,j] = atomsOrigins[2,j] - L
            
        end
    end

    ### Plotting Output
    Plots.scatter(atoms_new[1,:],atoms_new[2,:], label = "Step: $count" ,xlim = (0,L),ylim = (0,L),  title = "2-Dimensional Molecular Dynamics", background_color = :white, seriescolor = Colours[:], size=(1920,1080))    
   
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

    ### Per Atom Energy Calculation
    for i in 1:NAtoms
        KE[i] = 0.5 * ((xVelocity[i]^2) + (yVelocity[i]^2))
        PE[i] = sqrt(((sum(xForce[i,:]))^2)+((sum(yForce[i,:]))^2))
    end

    ### Average Energy Calculation
    KineticEnergy[count] = (sum(KE)/NAtoms)
    PotentialEnergy[count] = (sum(PE)/NAtoms)
    TotalEnergy[count] = KineticEnergy[count] + PotentialEnergy[count]

    ### Setting Atom Position for Next Timestep
    for i in 1:NAtoms+numFiring
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

### Saving all Revelant Simulation Data

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

plot(AvgMSD,title = "Δr Squared", label = false,size = (1920,1080))
png("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/AvgMSD")

# gif(anim, "$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/MolDynSlow.gif", fps = 24)
