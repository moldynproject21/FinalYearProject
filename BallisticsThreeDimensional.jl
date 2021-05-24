### Imports
using Plots
using Dates
using Random
using DataFrames
using CSV

println("Start Time:  ", now())
t1 = time()

### User Variables
NAtomsIsh = 1000
L = 30
δr = 0.1
Vo = 1
Δt = 1e-3
rCut = 3
iters = 12500

InterAtomDist = 1.1

### Heating Parameters
Heatbath = true
InterHeatingIters = 500
HeatFactor = 0.99

### Ballistics Parameters
Ballistics = true
FireAfter = 1000
M = 1
B2 = 4e-5
BM = B2/M
g = 9.8
FiringX = 3.5
FiringY = 23
Range = 1

numFiring = 7

PBCBound = 0.9
PBCDist = PBCBound * L

### Seeding random number generator - Comment out if unneeded
#Random.seed!(2142)


### Setting Up Simulation Parameters
CBRAtoms = trunc(Int,ceil(cbrt(NAtomsIsh)))
NAtoms = (CBRAtoms * CBRAtoms * CBRAtoms)
SQAtoms = CBRAtoms * CBRAtoms

xAtoms = zeros(NAtoms+numFiring)
xVelocity = zeros(NAtoms+numFiring)
xForce = zeros(NAtoms+numFiring,NAtoms+numFiring)

yAtoms = zeros(NAtoms+numFiring)
yVelocity = zeros(NAtoms+numFiring)
yForce = zeros(NAtoms+numFiring,NAtoms+numFiring)

zAtoms = zeros(NAtoms+numFiring)
zVelocity = zeros(NAtoms+numFiring)
zForce = zeros(NAtoms+numFiring,NAtoms+numFiring)

r = zeros(NAtoms+numFiring,NAtoms+numFiring)
Theta = zeros(NAtoms+numFiring,NAtoms+numFiring)
Phi = zeros(NAtoms+numFiring,NAtoms+numFiring)

KineticEnergy = zeros(iters)
PotentialEnergy = zeros(iters)
TotalEnergy = zeros(iters)
KE = zeros(NAtoms+numFiring)
PE = zeros(NAtoms+numFiring)

atoms = zeros(3,NAtoms+numFiring)
atoms_prev = zeros(3,NAtoms+numFiring)
atoms_new = zeros(3,NAtoms+numFiring)

atomsOrigins = zeros(3,NAtoms+numFiring)
MSD = zeros(NAtoms+numFiring)

AvgMSD = zeros(iters)



today = Dates.today()
mylocale = pwd()
hour = Dates.hour(now())
minute = Dates.minute(now())



Colours = []
for i in 1:(NAtoms+numFiring)
    if i <= NAtoms
        push!(Colours,:magenta3)
    else
        push!(Colours, :yellow)
    end
end



### Creating Grid of Atoms
Offset = ((L - ((CBRAtoms-1)*(InterAtomDist)))/2)-InterAtomDist

Gx = zeros(CBRAtoms)
for i in 1:CBRAtoms
    Gx[i] = (i*InterAtomDist)+Offset
end
Gx = repeat(Gx; outer=[CBRAtoms])
Gx = repeat(Gx; outer=[CBRAtoms])

Gy = (((1:CBRAtoms)*InterAtomDist).+Offset)' .* (ones(CBRAtoms))
Gy = reshape(Gy, (:))
Gy = repeat(Gy; outer=[CBRAtoms])

Gz = (((1:CBRAtoms)*InterAtomDist).+Offset)' .* (ones(SQAtoms))
Gz = reshape(Gz, (1, :))
Gz = Gz .- 7

### Setting Initial Positions and Velocities
for i in 1:NAtoms
    xVelocity[i] = 2*(rand()-0.5)*Vo
    yVelocity[i] = 2*(rand()-0.5)*Vo
    zVelocity[i] = 2*(rand()-0.5)*Vo

    atoms[1,i] = Gx[i] +(2*(rand()-0.5)*δr)
    atoms[2,i] = Gy[i] +(2*(rand()-0.5)*δr)
    atoms[3,i] = Gz[i] +(2*(rand()-0.5)*δr)

    atoms_prev[1,i] = atoms[1,i] - (xVelocity[i]*Δt)
    atoms_prev[2,i] = atoms[2,i] - (yVelocity[i]*Δt)
    atoms_prev[3,i] = atoms[3,i] - (zVelocity[i]*Δt)

    atomsOrigins[1,i] = atoms[1,i]
    atomsOrigins[2,i] = atoms[2,i]
    atomsOrigins[3,i] = atoms[3,i]
end

### Ballistics Placement - Manually Placing Bullets In Spiked Shape
withinRange = false
## Bullet 1 - Center
atoms_prev[1,NAtoms+1] = 0.001
atoms[1,NAtoms+1] = 0.001
atoms_new[1,NAtoms+1] = 0.001

atoms_prev[2,NAtoms+1] = L/2 
atoms[2,NAtoms+1] = L/2 
atoms_new[2,NAtoms+1] = L/2

atoms_prev[3,NAtoms+1] = 0.001
atoms[3,NAtoms+1] = 0.001
atoms_new[3,NAtoms+1] = 0.001

xVelocity[NAtoms+1] = FiringX
zVelocity[NAtoms+1] = FiringY
## Bullet 2 Negative Y
atoms_prev[1,NAtoms+2] = 0.001
atoms[1,NAtoms+2] = 0.001
atoms_new[1,NAtoms+2] = 0.001

atoms_prev[2,NAtoms+2] = L/2 - 1.2
atoms[2,NAtoms+2] = L/2 - 1.2
atoms_new[2,NAtoms+2] = L/2 - 1.2

atoms_prev[3,NAtoms+2] = 0.001
atoms[3,NAtoms+2] = 0.001
atoms_new[3,NAtoms+2] = 0.001

xVelocity[NAtoms+2] = FiringX
zVelocity[NAtoms+2] = FiringY
## Bullet 3 Positive Y
atoms_prev[1,NAtoms+3] = 0.001
atoms[1,NAtoms+3] = 0.001
atoms_new[1,NAtoms+3] = 0.001

atoms_prev[2,NAtoms+3] = L/2 + 1.2
atoms[2,NAtoms+3] = L/2 + 1.2
atoms_new[2,NAtoms+3] = L/2 + 1.2

atoms_prev[3,NAtoms+3] = 0.001
atoms[3,NAtoms+3] = 0.001
atoms_new[3,NAtoms+3] = 0.001

xVelocity[NAtoms+3] = FiringX
zVelocity[NAtoms+3] = FiringY
## Bullet 4 Negative X
atoms_prev[1,NAtoms+4] = 0.001 - 1.2
atoms[1,NAtoms+4] = 0.001 - 1.2
atoms_new[1,NAtoms+4] = 0.001 - 1.2

atoms_prev[2,NAtoms+4] = L/2
atoms[2,NAtoms+4] = L/2
atoms_new[2,NAtoms+4] = L/2

atoms_prev[3,NAtoms+4] = 0.001
atoms[3,NAtoms+4] = 0.001
atoms_new[3,NAtoms+4] = 0.001

xVelocity[NAtoms+4] = FiringX
zVelocity[NAtoms+4] = FiringY
## Bullet 5 Positive X
atoms_prev[1,NAtoms+5] = 0.001 + 1.2
atoms[1,NAtoms+5] = 0.001 + 1.2
atoms_new[1,NAtoms+5] = 0.001 + 1.2

atoms_prev[2,NAtoms+5] = L/2
atoms[2,NAtoms+5] = L/2
atoms_new[2,NAtoms+5] = L/2

atoms_prev[3,NAtoms+5] = 0.001
atoms[3,NAtoms+5] = 0.001
atoms_new[3,NAtoms+5] = 0.001

xVelocity[NAtoms+5] = FiringX
zVelocity[NAtoms+5] = FiringY
## Bullet 6 Negative Z
atoms_prev[1,NAtoms+6] = 0.001
atoms[1,NAtoms+6] = 0.001
atoms_new[1,NAtoms+6] = 0.001

atoms_prev[2,NAtoms+6] = L/2
atoms[2,NAtoms+6] = L/2
atoms_new[2,NAtoms+6] = L/2

atoms_prev[3,NAtoms+6] = 0.001 - 1.2
atoms[3,NAtoms+6] = 0.001 - 1.2
atoms_new[3,NAtoms+6] = 0.001 - 1.2

xVelocity[NAtoms+6] = FiringX
zVelocity[NAtoms+6] = FiringY
## Bullet 7 Positive Z
atoms_prev[1,NAtoms+7] = 0.001
atoms[1,NAtoms+7] = 0.001
atoms_new[1,NAtoms+7] = 0.001

atoms_prev[2,NAtoms+7] = L/2 
atoms[2,NAtoms+7] = L/2
atoms_new[2,NAtoms+7] = L/2

atoms_prev[3,NAtoms+7] = 0.001 + 1.2
atoms[3,NAtoms+7] = 0.001 + 1.2
atoms_new[3,NAtoms+7] = 0.001 + 1.2

xVelocity[NAtoms+7] = FiringX
zVelocity[NAtoms+7] = FiringY







mkdir("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations")
mkdir("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/Screens")


Gx = nothing
Gy = nothing

Plots.pyplot()

##### Force Animation
anim = @animate for count ∈ 1:iters

    ### Calculating Distances
    Threads.@threads for x in 1:NAtoms+numFiring
        for y in 1:NAtoms+numFiring
            ####Original Box
            r[x,y] = sqrt(((atoms[1,x]-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2)+((atoms[3,x]-atoms[3,y])^2))
            Theta[x,y] = atan((atoms[2,x]-atoms[2,y]),(atoms[1,x]-atoms[1,y]))
            Phi[x,y] = atan(sqrt(((atoms[1,x]-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2)),(atoms[3,x]-atoms[3,y]))

            #### Checking Periodic Boundary Atoms

            if r[x,y] > PBCDist #x + L
                r[x,y] = sqrt((((atoms[1,x]+L)-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2)+((atoms[3,x]-atoms[3,y])^2))
                Theta[x,y] = atan((atoms[2,x]-atoms[2,y]),((atoms[1,x]+L)-atoms[1,y]))
                Phi[x,y] = atan(sqrt((((atoms[1,x]+L)-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2)),(atoms[3,x]-atoms[3,y]))
            end

            if r[x,y] > PBCDist #x - L
                r[x,y] = sqrt((((atoms[1,x]-L)-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2)+((atoms[3,x]-atoms[3,y])^2))
                Theta[x,y] = atan((atoms[2,x]-atoms[2,y]),((atoms[1,x]-L)-atoms[1,y]))
                Phi[x,y] = atan(sqrt((((atoms[1,x]-L)-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2)),(atoms[3,x]-atoms[3,y]))
            end

            if r[x,y] > PBCDist #y + L
                r[x,y] = sqrt(((atoms[1,x]-atoms[1,y])^2)+(((atoms[2,x]+L)-atoms[2,y])^2)+((atoms[3,x]-atoms[3,y])^2))
                Theta[x,y] = atan(((atoms[2,x]+L)-atoms[2,y]),(atoms[1,x]-atoms[1,y]))
                Phi[x,y] = atan(sqrt(((atoms[1,x]-atoms[1,y])^2)+(((atoms[2,x]+L)-atoms[2,y])^2)),(atoms[3,x]-atoms[3,y]))
            end

            if r[x,y] > PBCDist #y - L
                r[x,y] = sqrt(((atoms[1,x]-atoms[1,y])^2)+(((atoms[2,x]-L)-atoms[2,y])^2)+((atoms[3,x]-atoms[3,y])^2))
                Theta[x,y] = atan(((atoms[2,x]-L)-atoms[2,y]),(atoms[1,x]-atoms[1,y]))
                Phi[x,y] = atan(sqrt(((atoms[1,x]-atoms[1,y])^2)+(((atoms[2,x]-L)-atoms[2,y])^2)),(atoms[3,x]-atoms[3,y]))
            end

            if r[x,y] > PBCDist #z + L
                r[x,y] = sqrt(((atoms[1,x]-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2)+(((atoms[3,x]+L)-atoms[3,y])^2))
                Theta[x,y] = atan((atoms[2,x]-atoms[2,y]),(atoms[1,x]-atoms[1,y]))
                Phi[x,y] = atan(sqrt(((atoms[1,x]-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2)),((atoms[3,x]+L)-atoms[3,y]))
            end

            if r[x,y] > PBCDist #z - L
                r[x,y] = sqrt(((atoms[1,x]-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2)+(((atoms[3,x]-L)-atoms[3,y])^2))
                Theta[x,y] = atan((atoms[2,x]-atoms[2,y]),(atoms[1,x]-atoms[1,y]))
                Phi[x,y] = atan(sqrt(((atoms[1,x]-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2)),((atoms[3,x]-L)-atoms[3,y]))
            end

            if r[x,y] > PBCDist #Back to center
                r[x,y] = sqrt(((atoms[1,x]-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2)+((atoms[3,x]-atoms[3,y])^2))
                Theta[x,y] = atan((atoms[2,x]-atoms[2,y]),(atoms[1,x]-atoms[1,y]))
                Phi[x,y] = atan(sqrt(((atoms[1,x]-atoms[1,y])^2)+((atoms[2,x]-atoms[2,y])^2)),(atoms[3,x]-atoms[3,y]))
            end

        end
    end

    ### Swapping From Ballistic to LJ
    if withinRange == false && minimum(r[(NAtoms+1:NAtoms+numFiring),(1:NAtoms)]) < Range
        global withinRange = true
        println("Switch to LJ in Iter:  ",  count)
        global Ballistics = false

        ### Adaptive Timestepping
        global Δt = 1e-4
    end

    ### Ballistic Motion Calculation
    if Ballistics == true && withinRange == false

        if (count > FireAfter) 
            for numBull in 1:numFiring
                atoms_new[1,NAtoms+numBull] = atoms[1,NAtoms+numBull] + (xVelocity[NAtoms+numBull]*Δt)
                atoms_new[3,NAtoms+numBull] = atoms[3,NAtoms+numBull] + (zVelocity[NAtoms+numBull]*Δt)

                v = sqrt((xVelocity[NAtoms+numBull]^2)+(zVelocity[NAtoms+numBull]^2))
        
                XDrag = - (BM)*(v)*(xVelocity[NAtoms+numBull])
                YDrag = - (BM)*(v)*(zVelocity[NAtoms+numBull])

                xVelocity[NAtoms+numBull] = xVelocity[NAtoms+numBull] + ((XDrag/M)*(Δt))
                zVelocity[NAtoms+numBull] = zVelocity[NAtoms+numBull] + ((-g + (YDrag/M))*(Δt))
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
                    zForce[x,y] = 0
                else
                    xForce[x,y] = ( 24*(((2/((r[x,y])^13)))-((1/((r[x,y])^7)))) ) * sin(Phi[x,y]) * cos(Theta[x,y])
                    yForce[x,y] = ( 24*(((2/((r[x,y]^13))))-((1/((r[x,y])^7)))) ) * sin(Phi[x,y]) * sin(Theta[x,y])
                    zForce[x,y] = ( 24*(((2/((r[x,y]^13))))-((1/((r[x,y])^7)))) ) * cos(Phi[x,y])
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
                atoms_prev[3,x] = atoms[3,x] - (HeatFactor*(atoms[3,x] - atoms_prev[3,x]))
            end
        end
    end




    ### Applying Force
    if withinRange == false
        for x in 1:NAtoms      
            atoms_new[1,x] = (2 * atoms[1,x]) - (atoms_prev[1,x]) + ((sum(xForce[x,:]))*Δt*Δt)
            atoms_new[2,x] = (2 * atoms[2,x]) - (atoms_prev[2,x]) + ((sum(yForce[x,:]))*Δt*Δt)
            atoms_new[3,x] = (2 * atoms[3,x]) - (atoms_prev[3,x]) + ((sum(zForce[x,:]))*Δt*Δt)
        end
    else
        for x in 1:NAtoms+numFiring      
            atoms_new[1,x] = (2 * atoms[1,x]) - (atoms_prev[1,x]) + ((sum(xForce[x,:]))*Δt*Δt)
            atoms_new[2,x] = (2 * atoms[2,x]) - (atoms_prev[2,x]) + ((sum(yForce[x,:]))*Δt*Δt)
            atoms_new[3,x] = (2 * atoms[3,x]) - (atoms_prev[3,x]) + ((sum(zForce[x,:]))*Δt*Δt)
        end
    end


    ### Periodic Boundary Correction
    for j in 1:NAtoms

            if atoms_new[1,j] < 0
            println("X Went Under 0, Iter: ", count, "    Atom: ", j)
            while atoms_new[1,j] < 0          
                atoms_new[1,j] = mod(atoms_new[1,j],L)

                atoms[1,j] = atoms[1,j] + L
                atoms_prev[1,j] = atoms_prev[1,j] + L
                atomsOrigins[1,j] = atomsOrigins[1,j] + L
            end
        end

        if atoms_new[1,j] > L
            println("X Went Over L, Iter: ", count, "    Atom: ", j)
            atoms_new[1,j] = atoms_new[1,j] % L

            atoms[1,j] = atoms[1,j] - L
            atoms_prev[1,j] = atoms_prev[1,j] - L
            atomsOrigins[1,j] = atomsOrigins[1,j] - L
        end



        if atoms_new[2,j] < 0
            println("Y Went Under 0, Iter:  ", count, "    Atom: ", j)
            while atoms_new[2,j] < 0
                atoms_new[2,j] = mod(atoms_new[2,j],L)

                atoms[2,j] = atoms[2,j] + L
                atoms_prev[2,j] = atoms_prev[2,j] + L
                atomsOrigins[2,j] = atomsOrigins[2,j] + L
                
            end
        end

        if atoms_new[2,j] > L
            println("Y Went Over L, Iter: ", count, "    Atom: ", j)
            atoms_new[2,j] = atoms_new[2,j] % L

            atoms[2,j] = atoms[2,j] - L
            atoms_prev[2,j] = atoms_prev[2,j] - L
            atomsOrigins[2,j] = atomsOrigins[2,j] - L
        end

        

        if atoms_new[3,j] < 0
            println("Z Went Under 0, Iter:  ", count, "    Atom: ", j)
            while atoms_new[3,j] < 0
                atoms_new[3,j] = mod(atoms_new[3,j],L)

                atoms[3,j] = atoms[3,j] + L
                atoms_prev[3,j] = atoms_prev[3,j] + L
                atomsOrigins[3,j] = atomsOrigins[3,j] + L
            end
        end

        if atoms_new[3,j] > L
            println("Z Went Over L, Iter: ", count, "    Atom: ", j)
            atoms_new[3,j] = atoms_new[3,j] % L

            atoms[3,j] = atoms[3,j] - L
            atoms_prev[3,j] = atoms_prev[3,j] - L
            atomsOrigins[3,j] = atomsOrigins[3,j] - L
        end
    end

    ### Plotting Output
    Plots.scatter(atoms_new[1,:],atoms_new[2,:],atoms_new[3,:], label = "Step: $count" ,xlim = (0,L),ylim = (0,L), zlim = (0,L),  title = "3-Dimensional Molecular Dynamics", seriescolor = Colours[:], size=(1920,1080),camera=(((-30 - FireAfter) + (count/100)),15))    
    
    #### Saving a screen shot every n timesteps
    if (count % 10) == 0 || count == 1
        png("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/Screens/$count")
    end

    ### MSD Calculations
    for i in 1:NAtoms
        MSD[i] = (((atoms_new[1,i]-atomsOrigins[1,i])^2)+((atoms_new[2,i]-atomsOrigins[2,i])^2)+((atoms_new[3,i]-atomsOrigins[3,i])^2))
    end
    AvgMSD[count] = (sum(MSD)/NAtoms)

    ### Velocity Calculations
    for i in 1:NAtoms
        xVelocity[i] = (atoms_new[1,i] - atoms_prev[1,i])/(2*Δt)
        yVelocity[i] = (atoms_new[2,i] - atoms_prev[2,i])/(2*Δt)
        zVelocity[i] = (atoms_new[3,i] - atoms_prev[3,i])/(2*Δt)
    end

    ### Per Atom Energy Calculations
    for i in 1:NAtoms
        KE[i] = 0.5 * ((xVelocity[i]^2) + (yVelocity[i]^2) + (zVelocity[i]^2))
        PE[i] = sqrt(((sum(xForce[i,:]))^2)+((sum(yForce[i,:]))^2)+((sum(zForce[i,:]))^2))
    end

    ### Average Energy Calculations
    KineticEnergy[count] = (sum(KE)/NAtoms)
    PotentialEnergy[count] = (sum(PE)/NAtoms)
    TotalEnergy[count] = KineticEnergy[count]+PotentialEnergy[count]

    ### Setting Atom Position for Next Timestep
    for i in 1:NAtoms+numFiring
        atoms_prev[1,i] = 1*atoms[1,i]
        atoms_prev[2,i] = 1*atoms[2,i]
        atoms_prev[3,i] = 1*atoms[3,i]
        atoms[1,i] = 1*atoms_new[1,i]
        atoms[2,i] = 1*atoms_new[2,i]
        atoms[3,i] = 1*atoms_new[3,i]
    end

    ### Completion Indicator
    if mod(count,iters/10) == 0
        println("Iterations Done: ",count, "    Time:  ", now(), "    Minutes Taken: ", ((time()-t1)/60))
    end

end

t2 = time()

println("Seconds Taken: ", t2-t1)
println("In Minutes: ", (t2-t1)/60)
println("In Hours: ", ((t2-t1)/60)/60)

### Saving all Revelant Simulation Data

mp4(anim, "$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/MolDynVideo.mp4", fps = 24)

vels = [xVelocity,yVelocity,zVelocity]
Es = [KineticEnergy,PotentialEnergy,TotalEnergy]
MSDs = [AvgMSD,AvgMSD]

df = DataFrames.DataFrame(vels)
df1 = DataFrames.DataFrame(xForce)
df2 = DataFrames.DataFrame(yForce)
df3 = DataFrames.DataFrame(zForce)
df4 = DataFrames.DataFrame(atoms)
df5 = DataFrames.DataFrame(atoms_prev)
df6 = DataFrames.DataFrame(MSDs)
df7 = DataFrames.DataFrame(Es)

CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/Velocity.csv",df)
CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/XForce.csv",df1)
CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/YForce.csv",df2)
CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/ZForce.csv",df3)
CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/Atoms.csv",df4)
CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/AtomsPrev.csv",df5)
CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/MSDs.csv",df6)
CSV.write("$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/Energies.csv",df7)


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

gif(anim, "$mylocale/JuliaPlot/$today-$hour $minute - $NAtoms atoms - $iters iterations/MolDynSlow.gif", fps = 24)
