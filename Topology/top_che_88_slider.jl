using Comodo
using Comodo.LinearAlgebra
using Comodo.GLMakie
using SparseArrays
using Statistics
using Profile
using Printf 

# Finite element analysis 
function elem_stiff(nu)
  A11 = [12 3 -6 -3;3 12 3 0; -6 3 12 -3; -3 0 -3 12]
  A12 = [-6 -3 0 3;-3 -6 -3 -6; 0 -3 -6 3; 3 -6 3 -6]
  B11 = [-4 3 -2 9; 3 -4 -9 4; -2 -9 -4 -3; 9 4 -3 -4]
  B12 = [2 -3 4 -9; -3 2 9 -2; 4 9 2 3; -9 -2 3 2]
  KE = 1/(1-nu^2)/24*([A11 A12; A12' A11]+ nu*[B11 B12; B12' B11])
  return KE
end 

function K_map(nelx,nely) 
    grid_view=reshape(collect(1:nelx*nely),nely,nelx)        
    global_node_numbering =reshape(collect(1:(nelx+1)*(nely+1)),nely+1,nelx+1)
    edofmat = zeros(Int64, nelx*nely,8)
    for i in axes(grid_view,2), j in axes(grid_view,1)    
        elm = grid_view[j,i]
        n1= global_node_numbering[j,i]
        n2= global_node_numbering[j,i+1]
        edofmat[elm:elm,:] = [2*n1+1, 2*n1+2, 2*n2+1, 2*n2+2, 2*n2-1, 2*n2, 2*n1-1, 2*n1]
    end
    iK = reshape(kron(edofmat,ones(8,1))',64*nelx*nely)
    jK = reshape(kron(edofmat,ones(1,8))',64*nelx*nely)    
    return edofmat,iK,jK
end 

function BC(nelx,nely)
  F=zeros(Float64,2*(nely+1)*(nelx+1),1)
  U=zeros(Float64,2*(nely+1)*(nelx+1))
  F[2,1] = -1
  fixeddofs = push!(collect(1:2:2*(nely+1)),  2*(nely+1)*(nelx+1))
  alldofs= collect(1:2*(nelx+1)*(nely+1))
  freedofs = setdiff(alldofs,fixeddofs) 
return F,U,freedofs
end 

## FEA solve
##
function FEA_solve(iK,jK,nelx,nely,penal,xphs,E0,Emin,KE,F,U,freedofs)     
    sK = reshape((KE[:]).*(Emin.+xphs[:]'.^penal.*(E0-Emin)),64*nelx*nely)   
    K = sparse(iK,jK,sK)  
    K = (K+K')/2   
#    K2= cholesky(K)         
    #U = cholesky(K)\F    
#    U[freedofs] = K[freedofs,freedofs]\F[freedofs]  
    U[freedofs] = cholesky(K[freedofs,freedofs])\F[freedofs]  
    return U 
end

function filter_weights(rmin,nelx,nely)
    iH = ones(nelx*nely*(2*(Int64(ceil(rmin))-1)+1)^2)
    jH = ones(size(iH))
    sH = zeros(size(iH))
    k = 0
    for i1 = 1:nelx
      for j1 = 1:nely
        e1 = (i1-1)*nely+j1
        for i2 = max(i1-(Int64(ceil(rmin))-1),1):min(i1+(Int64(ceil(rmin))-1),nelx)
          for j2 = max(j1-(Int64(ceil(rmin))-1),1):min(j1+(Int64(ceil(rmin))-1),nely)
            e2 = (i2-1)*nely+j2
            k = k+1
            iH[k] = e1
            jH[k] = e2
            sH[k] = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2)) #Hei = max(0,r_min - \Delta(e,i) (centre to centre dist between ))
          end
        end
      end
    end
    H = sparse(iH,jH,sH) # He values at every filtered region for the elemnt 
    Hs = sum(H,dims=2) #normalised constatnts 
    return H ,Hs
end

function OC_88(dc,dv,x,volfrac,nelx,nely,ft,H,Hs)
    l1= 0
    l2= 1e9
    move = 0.2
    condition = (l2-l1)/(l1+l2)
    xnew = zeros(Float64,nely,nelx)
    xphs = zeros(Float64,nely,nelx)
    
    while condition >  1e-3
      lmid = 0.5*(l1+l2)
      #h= -dc./dv
      h = (-dc./dv)./lmid
      g = x.*(h.^0.5)
      xnew = max.(0, max.(x .- move, min.(1.0, min.(x .+ move, g))))
      if ft == 1 
         xphs = xnew
        elseif ft == 2 
         xphs[:] = (H*xnew[:])./Hs    
      end    
      if sum(xphs[:]) > volfrac*nelx*nely
        l1 =lmid
        else 
        l2 =lmid 
      end 
      condition = (l2-l1)/(l1+l2)      
     
    end
    return xnew,xphs
end

function top_opt_88_che(nelx, nely, volfrac, penal, rmin, ft, changeTol)
    E0 = 1.0
    Emin = 1e-9
    nu = 0.3
    x=fill(volfrac,nely,nelx)
    xphs = x   
    loop = 1
    change = 1.0
    KE = elem_stiff(nu)
    edofmat,iK,jK = K_map(nelx,nely)
    F,U,freedofs = BC(nelx,nely)
    H, Hs = filter_weights(rmin, nelx, nely)
    data = Vector{Matrix{Float64}}()
    push!(data,x) # Add initial    
    while change > changeTol    
      #fe analysis
      U = FEA_solve(iK, jK, nelx, nely, penal, xphs, E0, Emin, KE, F, U, freedofs)    
      ce = reshape(sum(U[edofmat]*KE.*U[edofmat],dims=2),nely,nelx) #equation 2 
      c = sum((Emin.+xphs).^penal*(E0-Emin).*ce) #compliance Objective funciton debug
      dc = -penal*(E0-Emin).*xphs.^(penal-1).*ce #objective function 
      dv = ones(nely,nelx)

      ##### sensitivity       
      if ft == 1
        dc[:] = H*(x[:].*dc[:])./Hs./max.(1e-3,x[:]) #sensitivity
      elseif ft == 2
        dc[:] =  H*(dc[:]./Hs)
        dv[:] = H*(dv[:]./Hs)  
      end

      ## optimally criteria update 
      xnew,xphs = OC_88(dc,dv,x,volfrac,nelx,nely,ft,H,Hs)
      change = maximum(abs.(xnew[:].-x[:]))
      x = xnew
      vol=mean(xphs[:])
      println(@sprintf("Iter. = %i, Change = %.6f, Obj.: = %.6f, Vol. = %.3f",loop, change, c, vol))  
      # display(heatmap(1 .- xphs, color=:greys, yflip=true))     
      loop += 1
      push!(data,xphs)     
    end  
    return xphs, data, loop
end

function final_figure(data,loop)
    f = Figure(size =(600,200))
    ax = Axis(f[1,1], aspect = DataAspect(), yreversed=true)
    hm = heatmap!(ax, data[1]', colormap=:Spectral, colorrange=(0.0,1.0))
    Colorbar(f[:,end+1],hm)
    hSlider = Slider(f[2,1], range = 1:loop, startvalue = 1, linewidth=30)
   
    on(hSlider.value) do i 
        hm[1] = data[i]'        
    end 
    display(f)
    slidercontrol(hSlider,ax)

    return f 
end 

nelx = 60
nely = 20
volfrac = 0.5
penal = 3
ft = 1
rmin = 2.4
changeTol = 1e-2
@time xphs,data,loop = top_opt_88_che(nelx, nely, volfrac, penal, rmin, ft, changeTol)

f = final_figure(data,loop)
display(f) 