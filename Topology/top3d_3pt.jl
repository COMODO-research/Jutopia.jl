using LinearAlgebra
using Statistics
using SparseArrays
using GLMakie
using Comodo
using Comodo.GeometryBasics

#3D FEA analysis

function elment_stiff_3d(nu)

    A = [32 6 -8 6 -6 4 3 -6 -10 3 -3 -3 -4 -8; -48 0 0 -24 24 0 0 0 12 -12 0 12 12 12];
    k = 1/144*A'*[1; nu];
  
#Six sub-matricrices for obtaining KE Matrix
    K1 = [k[1] k[2] k[2] k[3] k[5] k[5];
          k[2] k[1] k[2] k[4] k[6] k[7];
          k[2] k[2] k[1] k[4] k[7] k[6];
          k[3] k[4] k[4] k[1] k[8] k[8];
          k[5] k[6] k[7] k[8] k[1] k[2];
          k[5] k[7] k[6] k[8] k[2] k[1]];
    K2 = [k[9]  k[8]  k[12] k[6]  k[4]  k[7];
          k[8]  k[9]  k[12] k[5]  k[3]  k[5];
          k[10] k[10] k[13] k[7]  k[4]  k[6];
          k[6]  k[5]  k[11] k[9]  k[2]  k[10];
          k[4]  k[3]  k[5]  k[2]  k[9]  k[12]
          k[11] k[4]  k[6]  k[12] k[10] k[13]];
    K3 = [k[6]  k[7]  k[4]  k[9]  k[12] k[8];
          k[7]  k[6]  k[4]  k[10] k[13] k[10];
          k[5]  k[5]  k[3]  k[8]  k[12] k[9];
          k[9]  k[10] k[2]  k[6]  k[11] k[5];
          k[12] k[13] k[10] k[11] k[6]  k[4];
          k[2]  k[12] k[9]  k[4]  k[5]  k[3]];
    K4 = [k[14] k[11] k[11] k[13] k[10] k[10];
          k[11] k[14] k[11] k[12] k[9]  k[8];
          k[11] k[11] k[14] k[12] k[8]  k[9];
          k[13] k[12] k[12] k[14] k[7]  k[7];
          k[10] k[9]  k[8]  k[7]  k[14] k[11];
          k[10] k[8]  k[9]  k[7]  k[11] k[14]];
    K5 = [k[1] k[2]  k[8]  k[3] k[5]  k[4];
          k[2] k[1]  k[8]  k[4] k[6]  k[11];
          k[8] k[8]  k[1]  k[5] k[11] k[6];
          k[3] k[4]  k[5]  k[1] k[8]  k[2];
          k[5] k[6]  k[11] k[8] k[1]  k[8];
          k[4] k[11] k[6]  k[2] k[8]  k[1]];
    K6 = [k[14] k[11] k[7]  k[13] k[10] k[12];
          k[11] k[14] k[7]  k[12] k[9]  k[2];
          k[7]  k[7]  k[14] k[10] k[2]  k[9];
          k[13] k[12] k[10] k[14] k[7]  k[11];
          k[10] k[9]  k[2]  k[7]  k[14] k[7];
          k[12] k[2]  k[9]  k[11] k[7]  k[14]];
    KE = 1/((nu+1)*(1-2*nu))*[K1 K2 K3 K4; K2' K5 K6 K3';K3' K6 K5' K2'; K4 K3 K2 K1']

    return KE 
end
#user defined load DOF
function K_map_3d(nelx,nely,nelz)
    nele = nelx*nely*nelz
    nodegrd = reshape(1:(nely+1)*(nelx+1),nely+1,nelx+1)
    nodeids = reshape(nodegrd[1:end-1,1:end-1],nely*nelx,1)
    nodeidz = collect(0:(nely+1)*(nelx+1):(nelz-1)*(nely+1)*(nelx+1))
    nodeids = kron(nodeids,ones(size(nodeidz))') + kron(nodeidz',ones(size(nodeids)[1]))   
    edofVec = 3*nodeids[:] .+ 1
    edofMat = Matrix{Int}
    edofMat = kron(edofVec,ones(1,24)) .+ kron([0 1 2 3*nely .+ [3 4 5 0 1 2] -3 -2 -1 3*(nely+1)*(nelx+1) .+ [0 1 2 3*nely .+ [3 4 5 0 1 2 ] -3 -2 -1]],ones(nele,1))
    iK = reshape(kron(edofMat,ones(24,1))',24*24*nele)
    iK = Int.(iK)
    jK = reshape(kron(edofMat,ones(1,24))',24*24*nele)
    jK = Int.(jK)
   
    return iK,jK,edofMat
end
function Bcs(nelx,nely,nelz)
    il = nelx/2
    jl = nelz
    kl = nelz/2
    loadnid = kl*(nelx+1)*(nely+1)+il*(nely+1)+(nely+1 - jl)
    loaddof = Int.(3*loadnid - 1)
#user defined support fixed dofs 
    
    iif = [0 0 nelx nelx]
    jf = [0 0 0 0]
    kf = [0 nelz 0 nelz]
    fixednid = kf*(nelx+1)*(nely+1)+iif*(nely+1)+(nely+1 .- jf)
    fixeddof = [3*fixednid[:]; 3*fixednid[:] .- 1; 3*fixednid[:] .- 2]  
    ndof = 3*(nelx+1)*(nely+1)*(nelz+1)
    F = zeros(ndof,1)
    F[loaddof] = -1
    U = zeros(ndof,1)
    freedofs = setdiff(collect(1:ndof),fixeddof)
   return freedofs,F,U,loaddof
end   
#preparing filter 
function filter_Hs(rmin, nelx,nely,nelz)


    rf = ceil(Int64,rmin)-1
    # iH = ones(Int64,nele*(2*rf+1)^2)    
    # iH = ones(Int(nele*(2*rf+1)^2),1)
    # jH = ones(Int64,size(iH))
    sH=Float64[]
    iH=Int[]
    jH=Int[]
    # sH = zeros(Float64,size(iH))
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
    KE = elment_stiff_3d(nu)
    sK = reshape(reshape(KE,576).*(Emin .+ reshape(xPhs,1,nele).^penal.*(E0-Emin)),576*nele)     
    K = sparse(iK,jK,sK)
    K = (K+K')/2
    U[freedofs] = K[freedofs,freedofs]\F[freedofs]    
   
    return U
end
function OC_3d(nelx,nely,nelz,dc,dv,x,volfrac,H,Hs)
    l1 = 0 
    l2 = 1e9 
    move = 0.2
    condition = (l2-l1)/(l1+l2)
    xnew = zeros(Float64,nely,nelx,nelz)
    xPhs = zeros(Float64,nely,nelx,nelz)
    #nele = nelx*nely*nelz
    while condition> 1e-3
        lmid = 0.5*(l2+l1)
        h = (-dc./dv)./lmid
        g = x.*(h.^0.5)
        xnew = max.(0, max.(x .- move, min.(1.0, min.(x .+ move, g))))   
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

function top_opt_3d(nelx,nely,nelz, volfrac,penal, rmin)
    E0 = 1
    Emin = 1e-9 
    nu = 0.3
    nele = nelx*nely*nelz
    x = fill(volfrac,(nely,nelx,nelz))
    xPhs = x
    loop =0
    change =1
    KE = elment_stiff_3d(nu)
    iK,jK,edofMat = K_map_3d(nelx,nely,nelz)
    freedofs,F,U,loaddof = Bcs(nelx,nely,nelz)
    Hs, H = filter_Hs(rmin, nelx,nely,nelz)
    
    #Iteration 
    maxloop = 200
    tolx = 0.01
    displayflag = 0  
    while change > tolx && loop < maxloop 
        loop +=1 
        U = fea_3d(xPhs,Emin,E0,nu,iK,jK,freedofs,penal,nelx,nelz,nely,F,U,loaddof)
        #Objective function
        ce = reshape(sum((U[Int.(edofMat)]*KE).*U[Int.(edofMat)],dims=2),(nely,nelx,nelz))
        c = sum((Emin .+ xPhs.^penal*(E0-Emin)).*ce)
        
        dc = -penal*(E0 - Emin)*xPhs.^(penal - 1).*ce
        dv = ones(nely,nelx,nelz)
        # Filtering and modifications of sensitivites
        dc[:] = H * (dc[:]./ Hs)
        dv[:] = H * (dv[:]./ Hs)
        #optimality criteria 
        xnew,xPhs = OC_3d(nelx,nely,nelz,dc,dv,x,volfrac,H,Hs)
        
        change = maximum(abs.(xnew.-x))
        x=xnew
        vol = mean(xPhs)
        
        println("Inter = $loop, ","Change = $change, ", "Obj: = $c, ", "Vol = $vol, ")   
    end 

    return xPhs,loop
end

nelx = 60
nely = 10
nelz = 10
volfrac = 0.2
penal = 3
rmin = 1.5

boxEl = (nely,nelx,nelz)
E,V,_,_,_ = hexbox(boxEl,boxEl)

xPhs,loop = top_opt_3d(nelx,nely,nelz, volfrac,penal, rmin)

xPhs_E = reshape(xPhs,prod(boxEl)) # a col corresponding to each element 
# xPhs_F = repeat(xPhs_E,inner=6)

bool3D_E = xPhs_E .>= 0.4 #Boolean  
bool3D_F = repeat(bool3D_E,inner=6) #repeat for all faces (1 1 1 1 1 1 2 2 2 2 2 2 ...) hex = 6

E_keep = E[bool3D_E] #elements to keep 
xPhs_E_keep = xPhs_E[bool3D_E] 

F_keep = element2faces(E_keep)# get faces for desired elements
xPhs_F_keep = repeat(xPhs_E_keep,inner=6) # get Xphs for boundary faces 

boolBoundaryFaces_F_keep = occursonce(F_keep; sort_entries=true)

Fs = F_keep[boolBoundaryFaces_F_keep]
xPhs_Fs = xPhs_F_keep[boolBoundaryFaces_F_keep]

# Visualisation

Fss,Vss = separate_vertices(Fs,V)
Css_V = simplex2vertexdata(Fss,xPhs_Fs)

M = Comodo.GeometryBasics.Mesh(Vss,Fss)

fig = Figure(size=(1600,800))
#scene =LScene(fig[1,1])
#cam = Makie.Camera3D(scene.scene, projectiontype = Makie.Perspective)
ax1 = Axis3(fig[1, 1], aspect = :data, xlabel = "X", ylabel = "Y", zlabel = "Z", title = "Boundary faces with boundary markers for the hexahedral mesh")
hp2 = poly!(ax1,M,shading=FastShading, color=Css_V, transparency=false, overdraw=false,colormap=:tableau_red_blue,colorrange=(0,1))
#hp2 = poly!(scene,M, strokewidth=3,shading=FastShading,strokecolor=:black, color=Css_V, transparency=false, overdraw=false,colormap=:grays,colorrange=(0,1))

#center!(scene.scene)
Colorbar(fig[1, 2], hp2)
fig
