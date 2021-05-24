using CUDA
using Random
using Dates
using Test
using DataFrames
using CSV
using Plots

Plots.pyplot()

println("Start Time:  ", now())
t1 = time()

### GPU Threads that can be leveraged
num_threads = 1024

### User Variables
NAtomsIsh = 16
L = 5.0f0
δr = 0.1
Vo = 1
Δt = 1e-3
rCut = 3
iters = 5000

InterAtomDist = 1.2

PBCBound = 0.6
PBCDist = PBCBound * L

### Seeding random number generator - Comment out if unneeded
#Random.seed!(2412)


### Setting Up Simulation Parameters
SQAtoms = trunc(Int,ceil(sqrt(NAtomsIsh)))
NAtoms = (SQAtoms * SQAtoms ::Int64)


xAtoms = CUDA.zeros(NAtoms)
xVelocity = CUDA.zeros(NAtoms)
ForceX = CUDA.zeros(NAtoms,NAtoms)

yAtoms = CUDA.zeros(NAtoms)
yVelocity = CUDA.zeros(NAtoms)
ForceY = CUDA.zeros(NAtoms,NAtoms)

dist = CUDA.zeros(NAtoms,NAtoms)
angles = CUDA.zeros(NAtoms,NAtoms)
xDif = CUDA.zeros(NAtoms,NAtoms)
yDif = CUDA.zeros(NAtoms,NAtoms)



distOG = CUDA.zeros(NAtoms,NAtoms)
anglesOG = CUDA.zeros(NAtoms,NAtoms)
dist1 = CUDA.zeros(NAtoms,NAtoms)
angles1 = CUDA.zeros(NAtoms,NAtoms)
dist2 = CUDA.zeros(NAtoms,NAtoms)
angles2 = CUDA.zeros(NAtoms,NAtoms)
dist3 = CUDA.zeros(NAtoms,NAtoms)
angles3 = CUDA.zeros(NAtoms,NAtoms)
dist4 = CUDA.zeros(NAtoms,NAtoms)
angles4 = CUDA.zeros(NAtoms,NAtoms)
dist5 = CUDA.zeros(NAtoms,NAtoms)
angles5 = CUDA.zeros(NAtoms,NAtoms)
dist6 = CUDA.zeros(NAtoms,NAtoms)
angles6 = CUDA.zeros(NAtoms,NAtoms)
dist7 = CUDA.zeros(NAtoms,NAtoms)
angles7 = CUDA.zeros(NAtoms,NAtoms)
dist8 = CUDA.zeros(NAtoms,NAtoms)
angles8 = CUDA.zeros(NAtoms,NAtoms)


KineticEnergy = zeros(iters)
PotentialEnergy = zeros(iters)
TotalEnergy = zeros(iters)
KE = CUDA.zeros(NAtoms)
PE = CUDA.zeros(NAtoms)



atomsX = CUDA.zeros(NAtoms)
atomsY = CUDA.zeros(NAtoms)
atoms_prevX = CUDA.zeros(NAtoms)
atoms_prevY = CUDA.zeros(NAtoms)
atoms_newX = CUDA.zeros(NAtoms)
atoms_newY = CUDA.zeros(NAtoms)


atomsOriginsX = CUDA.zeros(NAtoms)
atomsOriginsY = CUDA.zeros(NAtoms)
MSD = CUDA.zeros(NAtoms)
AvgMSD = zeros(iters)




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

### Setting Atoms Initial Positions
xVelocity[:] .= 2*(rand()-0.5)*Vo
yVelocity[:] .= 2*(rand()-0.5)*Vo

atomsX .= Gx[:] .+ (2*(CUDA.rand(length(xAtoms)).-0.5)*δr)
atomsY .= Gy[:] .+ (2*(CUDA.rand(length(yAtoms)).-0.5)*δr)

atoms_prevX = atomsX - (xVelocity[:] * Δt) 
atoms_prevY = atomsY - (yVelocity[:] * Δt)

atomsOriginsX = atomsX
atomsOriginsY = atomsY


### Functions To Run Simulations on GPU

function gpu_dif!(r, y, x)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(y)
        @inbounds r[i] = y[i] - x[i]
    end
    return
end


function gpu_dist!(r, y, x)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(y)

        @inbounds r[i] = CUDA.sqrt((CUDA.pow(x[i],2))+(CUDA.pow(y[i],2)))

    end
    return
end

function gpu_PBC!(r, ang, y, x, L, PBCDist)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(r)

        if (r[i] > PBCDist) ### X + L
            r[i] = CUDA.sqrt((CUDA.pow((x[i]+L),2))+(CUDA.pow(y[i],2)))
            ang[i] = CUDA.atan(y[i],(x[i]+L))
        end

        if (r[i] > PBCDist) ### X - L
            r[i] = CUDA.sqrt((CUDA.pow((x[i]-L),2))+(CUDA.pow(y[i],2)))
            ang[i] = CUDA.atan(y[i],(x[i]-L))
        end

        if (r[i] > PBCDist) ### Y + L
            r[i] = CUDA.sqrt((CUDA.pow(x[i],2))+(CUDA.pow((y[i]+L),2)))
            ang[i] = CUDA.atan((y[i]+L),x[i])
        end

        if (r[i] > PBCDist) ### Y - L
            r[i] = CUDA.sqrt((CUDA.pow(x[i],2))+(CUDA.pow((y[i]-L),2)))
            ang[i] = CUDA.atan((y[i]-L),x[i])
        end



        if (r[i] > PBCDist) ### X + L, Y + L
            r[i] = CUDA.sqrt((CUDA.pow((x[i]+L),2))+(CUDA.pow((y[i]+L),2)))
            ang[i] = CUDA.atan((y[i]+L),(x[i]+L))
        end

        if (r[i] > PBCDist) ### X + L, Y - L
            r[i] = CUDA.sqrt((CUDA.pow((x[i]+L),2))+(CUDA.pow((y[i]-L),2)))
            ang[i] = CUDA.atan((y[i]-L),(x[i]+L))
        end

        if (r[i] > PBCDist) ### X - L, Y + L
            r[i] = CUDA.sqrt((CUDA.pow((x[i]-L),2))+(CUDA.pow((y[i]+L),2)))
            ang[i] = CUDA.atan((y[i]+L),(x[i]-L))
        end

        if (r[i] > PBCDist) ### X - L, Y - L
            r[i] = CUDA.sqrt((CUDA.pow((x[i]-L),2))+(CUDA.pow((y[i]-L),2)))
            ang[i] = CUDA.atan((y[i]-L),(x[i]-L))
        end

        if (r[i] > PBCDist) ### Back to Center
            r[i] = CUDA.sqrt((CUDA.pow(x[i],2))+(CUDA.pow(y[i],2)))
            ang[i] = CUDA.atan(y[i],x[i])
        end

    end
    return
end


function gpu_angle!(r, y, x)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(y)

        @inbounds r[i] = CUDA.atan(y[i],x[i])

    end
    r = r'
    return
end




function gpu_forceXCalc!(f, d, a)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(d)

        @inbounds f[i] = (24 * (((2 / ((CUDA.pow(d[i],13)))) - ((1 / ((CUDA.pow.(d[i],7))))) ))) * CUDA.cos(a[i])

    end
    return
end


function gpu_forceYCalc!(f, d, a)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(d)

        @inbounds f[i] = (24 * (((2 / ((CUDA.pow(d[i],13)))) - ((1 / ((CUDA.pow.(d[i],7))))) ))) * CUDA.sin(a[i])

    end
    return
end



function gpu_move!(an, ac, ap, f)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(an)

        @inbounds an[i] = (2 * ac[i]) - (ap[i]) + ((f[i]))

    end
    return
end



function gpu_VelCalc!(v, an, ap)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(v)

        @inbounds v[i] = ((an[i]) - (ap[i]))

    end
    return
end

function gpu_KeCalc!(k, vx, vy)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(k)

        @inbounds k[i] =  0.5 * ((CUDA.pow(vx[i],2)) + (CUDA.pow(vy[i],2)))

    end
    return
end

function gpu_PeCalc!(p, fx, fy)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(p)

        @inbounds p[i] =  sqrt((CUDA.pow(fx[i],2)) + (CUDA.pow(fy[i],2)))

    end
    return
end


function gpu_MSDCalc!(m, xn, xo, yn, yo)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(m)

        @inbounds m[i] =  CUDA.pow((xn[i] - xo[i]),2) + CUDA.pow((yn[i] - yo[i]),2)

    end
    return
end



function gpu_OverCheck!(an, a, ap, ao, L)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(an)

        @inbounds if (an[i] > L)
            an[i] %= L
            a[i] += (-L) 
            ap[i] += (-L) 
            ao[i] += (-L) 
        end
                    
    end
    return
end


function gpu_UnderCheck!(an, a, ap, ao, L)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(an)

        @inbounds if (an[i] < 0)
            while (an[i] < 0)
                j = an[i]
                an[i] = L - (CUDA.mod(j,L))
            end     
            a[i] += L
            ap[i] += L
            ao[i] += L
        end

    end
    return
end


function gpu_rCheck!(r, rcut, fx, fy)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(r)
        @inbounds if (r[i] > rcut)
            fx[i] = 0
            fy[i] = 0
        end
                    
    end
    return
end

function gpu_divZeroCheck!(r, rcut, fx, fy)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(r)
        @inbounds if (r[i] == 0)
            fx[i] = 0
            fy[i] = 0
        end
                    
    end
    return
end

function gpu_setAtoms!(anew, aold)
    index = (blockIdx().x - 1) * blockDim().x + threadIdx().x
    stride = blockDim().x * gridDim().x
    for i = index:stride:length(anew)

        @inbounds anew[i] = aold[i]

                    
    end
    return
end


numblocks = ceil(Int, NAtoms/num_threads)


## THE MAIN LOOP
anim = @animate for iterator in 1:iters
    
    ###Increase dimensions of atoms array - Used for Calculating Distance
    atomsDifX = (repeat((atomsX)', NAtoms))
    atomsDifY = (repeat((atomsY)', NAtoms))

    ### Find Distance in 2 Planes 
    @cuda threads=num_threads blocks=numblocks gpu_dif!(xDif, atomsDifX', atomsDifX)
    @cuda threads=num_threads blocks=numblocks gpu_dif!(yDif, atomsDifY', atomsDifY)

    ### Find Distance and Angle between Atoms
    @cuda threads=num_threads blocks=numblocks gpu_dist!(distOG, yDif, xDif)
    @cuda threads=num_threads blocks=numblocks gpu_angle!(anglesOG, yDif, xDif)
    ### Assign Resutls to Global Variables
    global dist = distOG
    global angles = anglesOG

    ### Find Distance and Angle between Atoms and their Periodic Boundary Pair Atoms
    ### Cardinals
    ### X + L
    @cuda threads=num_threads blocks=numblocks gpu_dif!(xDif, atomsDifX', (atomsDifX.+L))
    @cuda threads=num_threads blocks=numblocks gpu_dif!(yDif, atomsDifY', atomsDifY)
    @cuda threads=num_threads blocks=numblocks gpu_dist!(dist1, yDif, xDif)
    @cuda threads=num_threads blocks=numblocks gpu_angle!(angles1, yDif, xDif)

    ### X - L
    @cuda threads=num_threads blocks=numblocks gpu_dif!(xDif, atomsDifX', (atomsDifX.-L))
    @cuda threads=num_threads blocks=numblocks gpu_dif!(yDif, atomsDifY', atomsDifY)
    @cuda threads=num_threads blocks=numblocks gpu_dist!(dist2, yDif, xDif)
    @cuda threads=num_threads blocks=numblocks gpu_angle!(angles2, yDif, xDif)

    ### Y + L
    @cuda threads=num_threads blocks=numblocks gpu_dif!(xDif, atomsDifX', atomsDifX)
    @cuda threads=num_threads blocks=numblocks gpu_dif!(yDif, atomsDifY', (atomsDifY.+L))
    @cuda threads=num_threads blocks=numblocks gpu_dist!(dist3, yDif, xDif)
    @cuda threads=num_threads blocks=numblocks gpu_angle!(angles3, yDif, xDif)

    ### Y - L
    @cuda threads=num_threads blocks=numblocks gpu_dif!(xDif, atomsDifX', atomsDifX)
    @cuda threads=num_threads blocks=numblocks gpu_dif!(yDif, atomsDifY', (atomsDifY.-L))
    @cuda threads=num_threads blocks=numblocks gpu_dist!(dist4, yDif, xDif)
    @cuda threads=num_threads blocks=numblocks gpu_angle!(angles4, yDif, xDif)

    ### Corners
    ### X + L, Y + L
    @cuda threads=num_threads blocks=numblocks gpu_dif!(xDif, atomsDifX', (atomsDifX.+L))
    @cuda threads=num_threads blocks=numblocks gpu_dif!(yDif, atomsDifY', (atomsDifY.+L))
    @cuda threads=num_threads blocks=numblocks gpu_dist!(dist5, yDif, xDif)
    @cuda threads=num_threads blocks=numblocks gpu_angle!(angles5, yDif, xDif)

    ### X + L, Y - L
    @cuda threads=num_threads blocks=numblocks gpu_dif!(xDif, atomsDifX', (atomsDifX.+L))
    @cuda threads=num_threads blocks=numblocks gpu_dif!(yDif, atomsDifY', (atomsDifY.-L))
    @cuda threads=num_threads blocks=numblocks gpu_dist!(dist6, yDif, xDif)
    @cuda threads=num_threads blocks=numblocks gpu_angle!(angles6, yDif, xDif)

    ### X - L, Y + L
    @cuda threads=num_threads blocks=numblocks gpu_dif!(xDif, atomsDifX', (atomsDifX.-L))
    @cuda threads=num_threads blocks=numblocks gpu_dif!(yDif, atomsDifY', (atomsDifY.+L))
    @cuda threads=num_threads blocks=numblocks gpu_dist!(dist7, yDif, xDif)
    @cuda threads=num_threads blocks=numblocks gpu_angle!(angles7, yDif, xDif)

    ### X - L, Y - L
    @cuda threads=num_threads blocks=numblocks gpu_dif!(xDif, atomsDifX', (atomsDifX.-L))
    @cuda threads=num_threads blocks=numblocks gpu_dif!(yDif, atomsDifY', (atomsDifY.-L))
    @cuda threads=num_threads blocks=numblocks gpu_dist!(dist8, yDif, xDif)
    @cuda threads=num_threads blocks=numblocks gpu_angle!(angles8, yDif, xDif)


    ### Choose Pair that is shortest distance from atom
    d1 = findall(dist .> PBCDist)
    dist[d1] = dist1[d1]
    angles[d1] = angles1[d1]

    d2 = findall(dist .> PBCDist)
    dist[d2] = dist2[d2]
    angles[d2] = angles2[d2]

    d3 = findall(dist .> PBCDist)
    dist[d3] = dist3[d3]
    angles[d3] = angles3[d3]

    d4 = findall(dist .> PBCDist)
    dist[d4] = dist4[d4]
    angles[d4] = angles4[d4]

    d5 = findall(dist .> PBCDist)
    dist[d5] = dist5[d5]
    angles[d5] = angles5[d5]

    d6 = findall(dist .> PBCDist)
    dist[d6] = dist6[d6]
    angles[d6] = angles6[d6]

    d7 = findall(dist .> PBCDist)
    dist[d7] = dist7[d7]
    angles[d7] = angles7[d7]

    d8 = findall(dist .> PBCDist)
    dist[d8] = dist8[d8]
    angles[d8] = angles8[d8]

    dOG = findall(dist .> PBCDist)
    dist[dOG] = distOG[dOG]
    angles[dOG] = anglesOG[dOG]


    ### Calculate Force in X and Y Plane
    @cuda threads=num_threads blocks=numblocks gpu_forceXCalc!(ForceX, dist, angles)
    @cuda threads=num_threads blocks=numblocks gpu_forceYCalc!(ForceY, dist, angles)

    ### Remove Forces from greater than Cutoff Distance
    @cuda threads=num_threads blocks=numblocks gpu_rCheck!(dist, rCut, ForceX, ForceY)

    ### Remove Div 0 Error from atoms effect on themselves
    @cuda threads=num_threads blocks=numblocks gpu_divZeroCheck!(dist, rCut, ForceX, ForceY)

    ### Sum Forces to get total per atom
    LeForceX = (CUDA.sum(ForceX,dims= 2)) .* (Δt*Δt)
    LeForceY = (CUDA.sum(ForceY,dims= 2)) .* (Δt*Δt)

    ### Apply Forces
    @cuda threads=num_threads blocks=numblocks gpu_move!(atoms_newX, atomsX, atoms_prevX, LeForceX)
    @cuda threads=num_threads blocks=numblocks gpu_move!(atoms_newY, atomsY, atoms_prevY, LeForceY)

    ### Check if Atoms go Out of Boundary
    @cuda threads=num_threads blocks=numblocks gpu_OverCheck!(atoms_newX, atomsX, atoms_prevX, atomsOriginsX, L)
    @cuda threads=num_threads blocks=numblocks gpu_OverCheck!(atoms_newY, atomsY, atoms_prevY, atomsOriginsY, L)

    @cuda threads=num_threads blocks=numblocks gpu_UnderCheck!(atoms_newX, atomsX, atoms_prevX, atomsOriginsX, L)
    @cuda threads=num_threads blocks=numblocks gpu_UnderCheck!(atoms_newY, atomsY, atoms_prevY, atomsOriginsY, L)

    ## Plot Output --- Currently use CPU
    plotX = Array(atoms_newX)
    plotY = Array(atoms_newY)
    Plots.scatter(plotX, plotY, label = "Step: $iterator" ,xlim = (0,L),ylim = (0,L),  title = "2-Dimensional Molecular Dynamics on GPU", background_color = :black, seriescolor = [:dodgerblue], size=(1920,1080)) 

    ### Calculate Velocities
    @cuda threads=num_threads blocks=numblocks gpu_VelCalc!(xVelocity, atoms_newX, atoms_prevX)
    global xVelocity = xVelocity ./ (2*Δt)
    @cuda threads=num_threads blocks=numblocks gpu_VelCalc!(yVelocity, atoms_newY, atoms_prevY)
    global yVelocity = yVelocity ./ (2*Δt)


    ### Calculate Per Atom Energies
    @cuda threads=num_threads blocks=numblocks gpu_KeCalc!(KE, xVelocity, yVelocity)
    @cuda threads=num_threads blocks=numblocks gpu_PeCalc!(PE, ForceX, ForceY)
    @cuda threads=num_threads blocks=numblocks gpu_MSDCalc!(MSD, atoms_newX, atomsOriginsX, atoms_newY, atomsOriginsY)

    ### Calculate Average Energies
    AvgMSD[iterator] = (sum(MSD)/NAtoms)
    KineticEnergy[iterator] = (sum(KE)/NAtoms)
    PotentialEnergy[iterator] = (sum(PE)/NAtoms)
    TotalEnergy[iterator] = KineticEnergy[iterator] + PotentialEnergy[iterator]

    ## Set Atoms Position for Next Timestep
    @cuda threads=num_threads blocks=numblocks gpu_setAtoms!(atoms_prevX, atomsX)
    @cuda threads=num_threads blocks=numblocks gpu_setAtoms!(atoms_prevY, atomsY)

    @cuda threads=num_threads blocks=numblocks gpu_setAtoms!(atomsX, atoms_newX)
    @cuda threads=num_threads blocks=numblocks gpu_setAtoms!(atomsY, atoms_newY)

    ### Completion Indicator
    if mod(iterator,iters/10) == 0
        println("Iterations Done: ",iterator, "    Time:  ", now())
    end

end


t2 = time()

println("Seconds Taken: ", t2-t1)
println("In Minutes: ", (t2-t1)/60)
println("In Hours: ", ((t2-t1)/60)/60)

### Create Video Output of Simulation
mp4(anim, "JuliaPlot/GPU.mp4", fps = 24)
gif(anim, "JuliaPlot/GPU.gif", fps = 24)

