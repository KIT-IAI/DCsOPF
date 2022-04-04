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
## Global balancing
Ns = 300000
u, U = rand(T,N), zeros(T,T,N)
s, S = rand(T,N), zeros(T,T,N)
d, D = rand(T,N), zeros(T,T,N)
for n=1:N
    D[:,:,n] = buildLowerTri(T)
    U[:,:,n] = buildLowerTri(T)
    S[:,:,n] = buildLowerTri(T)
end

function a(t,n::Int64=1)
    term1 = sum( d[t,i]+u[t,i]+s[t,i] for i=1:N )
    term2 = sum( sum( D[t,k,i]*randn(n)+(U[t,k,i]+S[t,k,i])*sum(randn(n) for l=1:N) for k=1:t ) for i=1:N )
    term1*ones(n) + term2
end

function mean_a(t)
    sum( d[t,i]+u[t,i]+s[t,i] for i=1:N )
end

function var_a(t)
    sum( sum( D[t,k,i]^2+(S[t,k,i]+U[t,k,i])^2*N for k=1:t ) for i=1:N )
    # d = [ D[t,k,i] for i=1:N for k=1:t ]
    # su = [ sqrt(N)*(U[t,k,i]+S[t,k,i]) for i=1:N for k=1:t ]
end

function var_a_(t)
    # sum( sum( D[t,k,i]^2+(S[t,k,i]+U[t,k,i])^2*N for k=1:t ) for i=1:N )
    d = [ D[t,k,i] for i=1:N for k=1:t ]
    su = [ sqrt(N)*(U[t,k,i]+S[t,k,i]) for i=1:N for k=1:t ]
    return norm([d;su])^2
end

display("Mean of a:")
display(([mean(a(t,Ns)) for t=1:T],[mean_a(t) for t=1:T]))

display("Variance of a:")
display(([var(a(t,Ns)) for t=1:T],[var_a(t) for t=1:T]))
