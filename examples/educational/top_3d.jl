using Comodo
using Comodo.LinearAlgebra
using Comodo.GLMakie
using Comodo.GeometryBasics
using Statistics
using SparseArrays
using Profile 
using Printf 

#3D FEA analysis

function element_stiffness_3d(nu)
    A = [32 -48; 6 0; -8 0; 6 -24; -6 24; 4 0; 3 0; -6 0; -10 12; 3 -12; -3 0; -3 12; -4 12; -8 12]
    k = 1/144*A*[1; nu]
  
    #Six sub-matricrices for obtaining KE Matrix
    K1 = [k[1] k[2] k[2] k[3] k[5] k[5];
          k[2] k[1] k[2] k[4] k[6] k[7];
          k[2] k[2] k[1] k[4] k[7] k[6];
          k[3] k[4] k[4] k[1] k[8] k[8];
          k[5] k[6] k[7] k[8] k[1] k[2];
          k[5] k[7] k[6] k[8] k[2] k[1]]
    K2 = [k[9]  k[8]  k[12] k[6]  k[4]  k[7];
          k[8]  k[9]  k[12] k[5]  k[3]  k[5];
          k[10] k[10] k[13] k[7]  k[4]  k[6];
          k[6]  k[5]  k[11] k[9]  k[2]  k[10];
          k[4]  k[3]  k[5]  k[2]  k[9]  k[12]
          k[11] k[4]  k[6]  k[12] k[10] k[13]]
    K3 = [k[6]  k[7]  k[4]  k[9]  k[12] k[8];
          k[7]  k[6]  k[4]  k[10] k[13] k[10];
          k[5]  k[5]  k[3]  k[8]  k[12] k[9];
          k[9]  k[10] k[2]  k[6]  k[11] k[5];
          k[12] k[13] k[10] k[11] k[6]  k[4];
          k[2]  k[12] k[9]  k[4]  k[5]  k[3]]
    K4 = [k[14] k[11] k[11] k[13] k[10] k[10];
          k[11] k[14] k[11] k[12] k[9]  k[8];
          k[11] k[11] k[14] k[12] k[8]  k[9];
          k[13] k[12] k[12] k[14] k[7]  k[7];
          k[10] k[9]  k[8]  k[7]  k[14] k[11];
          k[10] k[8]  k[9]  k[7]  k[11] k[14]]
    K5 = [k[1] k[2]  k[8]  k[3] k[5]  k[4];
          k[2] k[1]  k[8]  k[4] k[6]  k[11];
          k[8] k[8]  k[1]  k[5] k[11] k[6];
          k[3] k[4]  k[5]  k[1] k[8]  k[2];
          k[5] k[6]  k[11] k[8] k[1]  k[8];
          k[4] k[11] k[6]  k[2] k[8]  k[1]]
    K6 = [k[14] k[11] k[7]  k[13] k[10] k[12];
          k[11] k[14] k[7]  k[12] k[9]  k[2];
          k[7]  k[7]  k[14] k[10] k[2]  k[9];
          k[13] k[12] k[10] k[14] k[7]  k[11];
          k[10] k[9]  k[2]  k[7]  k[14] k[7];
          k[12] k[2]  k[9]  k[11] k[7]  k[14]]
    KE = 1/((nu+1)*(1-2*nu))*[K1  K2  K3  K4; 
                              K2' K5  K6  K3';
                              K3' K6  K5' K2'; 
                              K4  K3  K2  K1']

    return KE 
end

#user defined load DOF
function K_map_3d(nelx,nely,nelz)
    nele = nelx*nely*nelz
    edofVec = [1+3*(LinearIndices((nely+1,nelx+1))[j,i]+k) for k in 0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1) for i in 1:nelx for j in 1:nely]
    
    edofMat = repeat(edofVec,1,24) + 
              repeat([0 1 2 3*nely .+ [3 4 5 0 1 2] -3 -2 -1 3*(nely+1)*(nelx+1) .+ [0 1 2 3*nely .+ [3 4 5 0 1 2 ] -3 -2 -1]], nele)
    
    iK = reshape(kron(edofMat, ones(Int,24,1))', 24*24*nele)    
    jK = reshape(kron(edofMat, ones(Int,1,24))', 24*24*nele)    
    
    return iK,jK,edofMat
end

function Bcs(nelx,nely,nelz)
    il = nelx .* fill(1,(1,1,nelz+1))
    jl = fill(0,(1,1,nelz+1))
    kl = fill(0,(1,1,nelz+1))
    for i in 0:nelz
        kl[1,1,i+1] = i 
    end 
    loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1 .- jl)
    loaddof = 3*loadnid[:] .- 1
    
    # User defined support fixed dofs 
    iif = fill(0,(nely+1,1,nelz+1))
    jf = fill(0,(nely+1,1,nelz+1))
    kf = fill(0,(nely+1,1,nelz+1))
    for i in 1:nelz+1
        jf[:,:,i] = 0:nely
    end
    for i in 0:nelz
        kf[:,:,i+1] .= i
    end
    fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1 .- jf)
    fixeddof = [3*fixednid[:]; 
                3*fixednid[:] .- 1; 
                3*fixednid[:] .- 2]  
    ndof = 3*(nelx+1)*(nely+1)*(nelz+1)
    F = zeros(ndof,1)
    F[loaddof] .= -1.0
    U = zeros(ndof,1)
    freedofs = setdiff(1:ndof,fixeddof)
   return freedofs,F,U,loaddof
end   

#preparing filter 
function filter_Hs(rmin, nelx,nely,nelz)
    rf = ceil(Int64,rmin)-1    
    iH = Vector{Int}()
    jH = Vector{Int}()
    sH = Vector{Float64}()
    k = 0 
    for kl in  range(1,nelz)
        for il in range(1,nelx)
            for jl in  range(1,nely)
                e1 = (kl-1)*nelx*nely + (il-1)*nely+jl
                for k2 in  range(max(kl-(rf),1),min(kl+(rf),nelz))
                    for i2 in  range(max(il-(rf),1),min(il+(rf),nelx))
                        for j2 in range(max(jl-(rf),1),min(jl+(rf),nely))
                            e2 = (k2-1)*nelx*nely + (i2-1)*nely+j2
                            k += 1
                            push!(iH,e1)
                            push!(jH,e2)
                            push!(sH,max(0,rmin-(sqrt((il-i2)^2+(jl-j2)^2+(kl-k2)^2))))
                            #   iH[k] = e1
                            #   jH[k] = e2
                            #   sH[k]  = max(0,rmin-(sqrt((il-i2)^2+(jl-j2)^2+(kl-k2)^2)))                        
                        end
                    end
                end
            end
        end
    end
    H = sparse(iH,jH,sH)
    Hs = sum(H,dims=2)
    return Hs,H
end

function fea_3d(xPhs,Emin,E0,nu,iK,jK,freedofs,penal,nelx,nelz,nely,F,U,loaddof)
    nele = nelx*nely*nelz
    KE = element_stiffness_3d(nu)
    sK = reshape(reshape(KE,576).*(Emin .+ reshape(xPhs,1,nele).^penal.*(E0-Emin)), 576*nele)     
    K = sparse(iK,jK,sK)
    K = (K+K') ./ 2.0
    U[freedofs] = K[freedofs,freedofs]\F[freedofs]       
    return U
end

function OC_3d(nelx,nely,nelz,dc,dv,x,volfrac,H,Hs)
    l1 = 0.0 
    l2 = 1e9 
    move = 0.2
    condition = (l2-l1)/(l1+l2)
    xnew = zeros(Float64,nely,nelx,nelz)
    xPhs = zeros(Float64,nely,nelx,nelz)
    # nele = nelx*nely*nelz
    while condition> 1e-3
        lmid = 0.5*(l2+l1)
        h = (-dc./dv)./lmid
        g = x.*(h.^0.5)
        xnew = max.(0.0, max.(x .- move, min.(1.0, min.(x .+ move, g))))   
        xPhs[:] = (H*xnew[:])./Hs
        if sum(xPhs[:]) > volfrac*nelx*nely*nelz
            l1 =lmid
        else 
            l2 =lmid  
        end
        condition = (l2-l1)/(l1+l2)
    end    
    return xnew,xPhs
end

function top_opt_3d(nelx, nely, nelz, volfrac, penal, rmin, changeTol, maxIter; E0=1.0, nu= 0.3, Emin=1e-9)
    
    # nele = nelx*nely*nelz
    x = fill(0.5,(nely,nelx,nelz))
    xPhs = x
    
    KE = element_stiffness_3d(nu)
    iK,jK,edofMat = K_map_3d(nelx,nely,nelz)
    freedofs,F,U,loaddof = Bcs(nelx,nely,nelz)
    Hs, H = filter_Hs(rmin, nelx,nely,nelz)
    
    # Start iterations
    loop = 0 # Iteration loop counter 
    change = 10.0 * changeTol # Change size initialised as larger than tolerance 
    while change > changeTol && loop < maxIter
        loop += 1 
        
        U = fea_3d(xPhs,Emin,E0,nu,iK,jK,freedofs,penal,nelx,nelz,nely,F,U,loaddof)
        
        # Objective function
        ce = reshape(sum((U[edofMat]*KE).*U[edofMat],dims=2),(nely,nelx,nelz))
        c = sum((Emin .+ xPhs.^penal*(E0-Emin)).*ce)
        
        dc = -penal*(E0 - Emin) .* xPhs.^(penal - 1).*ce
        dv = ones(nely,nelx,nelz)

        # Filtering and modification of sensitivities        
        dc[:] = H * (dc[:]./ Hs)
        dv[:] = H * (dv[:]./ Hs) # CHECK IF NEEDED

        # Optimality criteria 
        xnew,xPhs = OC_3d(nelx, nely, nelz, dc, dv, x, volfrac, H, Hs)
        
        change = maximum(abs.(xnew.-x))
        x = xnew
        vol = mean(xPhs)

        println(@sprintf("Iter. = %i, Change = %.6f, Obj.: = %.6f, Vol. = %.3f", loop, change, c, vol))         
    end 
    return xPhs,loop
end

nelx = 60
nely = 20
nelz = 4
rmin = 2.4
volfrac = 0.5
penal = 3.0
changeTol = 1e-2
maxIter = 16

boxEl = [nely,nelx,nelz]
E,V,_,_,_ = hexbox(boxEl,boxEl)

@time xPhs,loop =  top_opt_3d(nelx,nely,nelz, volfrac,penal, rmin, changeTol, maxIter)

xPhs_E = reshape(xPhs,prod(boxEl)) # a col corresponding to each element 
# xPhs_F = repeat(xPhs_E,inner=6)

bool3D_E = xPhs_E .>= 0.5 #Boolean  
bool3D_F = repeat(bool3D_E, inner=6) #repeat for all faces (1 1 1 1 1 1 2 2 2 2 2 2 ...) hex = 6

E_keep = E[bool3D_E] #elements to keep 
xPhs_E_keep = xPhs_E[bool3D_E] 

F_keep = element2faces(E_keep)# get faces for desired elements
xPhs_F_keep = repeat(xPhs_E_keep, inner=6) # get Xphs for boundary faces 

boolBoundaryFaces_F_keep = occursonce(F_keep; sort_entries=true)

Fs = F_keep[boolBoundaryFaces_F_keep]
xPhs_Fs = xPhs_F_keep[boolBoundaryFaces_F_keep]

# Visualisation

Fss,Vss = separate_vertices(Fs, V)
Css_V = simplex2vertexdata(Fss, xPhs_Fs)

fig = Figure(size=(1600,800))
ax1 = AxisGeom(fig[1, 1], title = "Boundary faces with boundary markers for the hexahedral mesh")
hp2 = meshplot!(ax1, Fss, Vss, color=Css_V, colormap=:tableau_red_blue, colorrange=(0,1), strokewidth=0.0)

#center!(scene.scene)
Colorbar(fig[1, 2], hp2)
fig
