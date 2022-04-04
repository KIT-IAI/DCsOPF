function buildLowerTri(n::Int64)
    X = zeros(n,n)
    for i=1:n
        for j=1:i
            X[i,j] = randn()
        end
    end
    return X
end
function ui(t::Int,ξ::Matrix{Float64},uhat::Vector{Float64},U::Array{Float64,3})
    T,N = size(ξ)
    ~,~,N_ = size(U)
    @assert N==N_ "wrong inputs"
    return uhat[t] + sum( sum( U[t,k,j]*ξ[k,j] for k=1:t ) for j=1:N )
end

function di(t::Int,ξ::Vector{Float64},d::Vector{Float64},D::Matrix{Float64})
    a, b = size(D)
    D_, ξ_ = zeros(a,b,1), zeros(length(ξ),1)
    D_[:,:,1], ξ_[:,1] = D, ξ
    ui(t,ξ_,d,D_)
end

function cl(ϕ,t,ξ::Matrix{Float64},d,D,u,U)
    N = length(ϕ)
    T,N_ = size(ξ)
    @assert N==N_
    sum( ϕ[i]*( di(t,ξ[:,i],d[:,i],D[:,:,i]) + ui(t,ξ,u[:,i],U[:,:,i,:]) ) for i=1:N )
end

function cl_samples(l::Int,Φ,d,D,u,U,Nsamples::Int=1000)
    Nline, N = size(Φ)
    T, N_ = size(d)
    @assert N==N_
    Cl = zeros(T,Nsamples)
    for i=1:Nsamples
        ξ = randn(T,N)
        Cl[:,i] = [ cl(Φ[l,:],t,ξ,d,D,u,U) for t=1:T ]
    end
    return Cl
end

function mean_cl(ϕ,t::Int,d,u)
    N = length(ϕ)
    T,N_ = size(d)
    @assert N==N_
    sum( ϕ[i]*( d[t,i]+u[t,i] ) for i=1:N )
end

function var_cl(ϕ,t::Int,D,U)
    N = length(ϕ)
    T,~,N_ = size(D)
    @assert N==N_
    sum( ( ϕ[i]*D[t,k,i] + sum(ϕ[j]*U[t,k,j,i] for j=1:N) )^2 for i=1:N, k=1:t)
end

function var_cl_glob(ϕ,t::Int,D,U)
    N = length(ϕ)
    T,~,N_ = size(D)
    @assert N==N_
    sum( ( ϕ[i]*D[t,k,i] + sum(ϕ[j]*U[t,k,j] for j=1:N) )^2 for i=1:N, k=1:t)
end

function var_cl_(ϕ,t::Int,D,U)
    N = length(ϕ)
    T,~,N_ = size(D)
    @assert N==N_
    VAR = zeros(N,T)
    for i=1:N
        for k=1:T
            VAR[i,k] = ϕ[i]*D[t,k,i]
            for j=1:N
                VAR[i,k] += ϕ[j]*U[t,k,j,i]
            end
        end
    end
    VAR_ = reshape(VAR,N*T,1)
    return norm(VAR_)^2
end



##
N, Nline, T = 4, 5, 10
Φ = ones(Nline,N)
dhat, D = rand(T,N), zeros(T,T,N)
[ D[:,:,i] = buildLowerTri(T) for i=1:N ]
uhat, U = rand(T,N), zeros(T,T,N,N)
[ U[:,:,i,j] = buildLowerTri(T) for i=1:N, j=1:N ]
# solution via sampling
Nsmpls, l = 100000, 3
Cl = cl_samples(l,Φ,dhat,D,uhat,U,Nsmpls)

display([ [ mean_cl(Φ[l,:],t,dhat,uhat) for t=1:T ] mean(Cl,2) ])
display([ [ var_cl(Φ[l,:],t,D,U) for t=1:T ] var(Cl,2) ])

## global balancing
N, Nline, T = 4, 5, 10
Φ = ones(Nline,N)
dhat, D = rand(T,N), zeros(T,T,N)
[ D[:,:,i] = buildLowerTri(T) for i=1:N ]
uhat, U, UU = rand(T,N), zeros(T,T,N), zeros(T,T,N,N)
[ U[:,:,i] = buildLowerTri(T) for i=1:N ]
[ UU[:,:,i,j] = U[:,:,i] for i=1:N for j=1:N ]
# solution via sampling
Nsmpls, l = 100000, 3
Cl = cl_samples(l,Φ,dhat,D,uhat,UU,Nsmpls)

display([ [ mean_cl(Φ[l,:],t,dhat,uhat) for t=1:T ] mean(Cl,2) ])
display([ [ var_cl_glob(Φ[l,:],t,D,U) for t=1:T ] var(Cl,2) ])
display([ [ var_cl(Φ[l,:],t,D,UU) for t=1:T ] var(Cl,2) ])
