N, T = 4, 10;
function buildLowerTri(n::Int64)
    X = zeros(n,n)
    for i=1:n
        for j=1:i
            X[i,j] = rand()
        end
    end
    return X
end
##

function mean_u(t::Int,u::Vector{Float64})
    u[t]
end

function var_u(t::Int,U::Array{Float64,3})
    T,~,N=size(U)
    sum( U[t,k,j]^2 for j=1:N for k=1:t )
end

function var_u_glob(t::Int,U::Array{Float64,2},N::Int64)
    # sum( N*U[t,k]^2 for k=1:t )
    v = [ sqrt(N)*U[t,k] for k=1:t ]
    return norm(v)^2
end

function var_u_(t::Int,U::Array{Float64,3})
    T,~,N=size(U)
    v = [ U[t,k,j] for j=1:N for k=1:t ]
    return norm(v)^2
end

function mean_Δu(t::Int,u::Vector{Float64})
    @assert t>=2
    mean_u(t,u)-mean_u(t-1,u)
end
function var_Δu(t::Int,U::Array{Float64,3})
    @assert t>=2
    T,~,N=size(U)
    sum( U[t,t,j]^2 + sum( (U[t,k,j]-U[t-1,k,j])^2 for k=1:t-1 ) for j=1:N)
end

function var_Δu_glob(t::Int,U::Array{Float64,2},N::Int64)
    @assert t>=2
    # N*(U[t,t]^2 + sum( (U[t,k]-Ui[t-1,k])^2 for k=1:t-1 ))
    v = [ sqrt(N)*U[t,t];
          sqrt(N)*[ U[t,k]-Ui[t-1,k] for k=1:t-1] ]
    return norm(v)^2
end

function var_Δu_(t::Int,U::Array{Float64,3})
    @assert t>=2
    T,~,N=size(U)
    v = [ [ U[t,k,j]-Ui[t-1,k,j] for k=1:t-1 for j=1:N ]; [U[t,t,j] for j=1:N] ]
    return norm(v)^2
end

function mean_e(t::Int,h::Float64,ic::Float64,s::Vector{Float64})
   if t==1 return ic
   else
       return ic + h*sum(s[k] for k=1:t-1)
   end
end
function var_e(t::Int,h::Float64,S::Array{Float64,3})
    if t==1 return 0.
    else
        T,~,N=size(S)
        return h^2*sum( sum(S[l,k,j] for l=k:t-1)^2  for k=1:t-1, j=1:N )
    end
end

function var_e_glob(t::Int,h::Float64,S::Array{Float64,2},N::Int64)
    if t==1 return 0.
    else
        # return h^2*N*sum( sum(S[l,k] for l=k:t-1)^2  for k=1:t-1)
        v = [ h*sqrt(N)*sum(S[l,k] for l=k:t-1)  for k=1:t-1 ]
        return norm(v)^2
    end
end

function var_e_(t::Int,h::Float64,S::Array{Float64,3})
    if t==1 return 0.
    else
        T,~,N=size(S)
        v = [ h*sum(S[l,k,j] for l=k:t-1)  for k=1:t-1, j=1:N  ]
        return dot(v,v)
    end
end



ui = collect(1.:T)
Ui = zeros(T,T,N)

Lu = buildLowerTri(T) # global balancing
for n=1:N
    # local balancing
    # Ui[:,:,n] = buildLowerTri(T)
    # global balancing
    Ui[:,:,n] = Lu
end

##########################################################################
##########################################################################
##########################################################################
u(ξ) = ui + sum( Ui[:,:,n]*ξ[:,n] for n=1:N )
Δu(t,ξ) = u(ξ)[t]-u(ξ)[t-1]

# moments for u
mean_u(t) = mean_u(t,ui)
var_u(t) = var_u(t,Ui)
var_u_glob(t) = var_u_glob(t,Lu,N)
# moments for Δu
mean_Δu(t) = mean_Δu(t,ui)
var_Δu(t) = var_Δu(t,Ui)
var_Δu_wrong(t) = var_u(t)+var_u(t-1)
var_Δu_glob(t) = var_Δu_glob(t,Lu,N)
# samples
Nsmpls = 100000;
U, ΔU = zeros(T,Nsmpls), zeros(T-1,Nsmpls)
for i=1:Nsmpls
    xi = randn(T,N)
    U[:,i] = u(xi)
    ΔU[:,i] = [ Δu(t,xi) for t=2:T ]
end
# compare
display("Mean of u:")
display([ [mean_u(t) for t=1:T] mean(U,2) ])
display("Variance of u:")
display([ [var_u(t) for t=1:T] var(U,2) ])

display("Mean of Δu:")
display([ [mean_Δu(t) for t=2:T] mean(ΔU,2) ])
display("Variance of Δu:")
display([ [var_Δu(t) for t=2:T] var(ΔU,2) ])
display("Variance of Δu (global):")
display([ [var_Δu_glob(t) for t=2:T] var(ΔU,2) ])

##########################################################################
##########################################################################
##########################################################################
si, Si = ones(T), zeros(T,T,N)
h = 0.1
Ls = buildLowerTri(T)
for n=1:N
    # Si[:,:,n] = buildLowerTri(T)
    Si[:,:,n] = Ls
end
s(ξ) = si + sum( Si[:,:,n]*ξ[:,n] for n=1:N )
eiIC = rand()
function ei_rec(t::Int,ξ)
    x = eiIC
    if t>1
        for t_=1:t-1
            x += h*s(ξ)[t_]
        end
    end
    return x
end
function ei(t::Int,ξ)
    if t==1
    return eiIC
    elseif t>1
        x = eiIC + h*sum(si[t_] for t_=1:t-1) + h*sum( sum(Si[l,k,j] for l=k:t-1)*ξ[k,j] for j=1:N, k=1:t-1 )
        return x
    end
end
E = zeros(T,Nsmpls)
for i=1:Nsmpls
    xi = randn(T,N)
    E[:,i] = [ ei(t,xi) for t=1:T ]
end
mean_ei(t) = mean_e(t,h,eiIC,si)
var_ei(t) = var_e(t,h,Si)
var_ei_glob(t) = var_e_glob(t,h,Ls,N)
## compare
display("Mean of e:")
display([ [ eiIC; [mean_ei(t) for t=2:T]] mean(E,2) ])
display("Variance of e:")
display([ [ 0; [var_ei(t) for t=2:T]] var(E,2) ])
display("Variance of e (global):")
display([ [ 0; [var_ei_glob(t) for t=2:T]] var(E,2) ])
