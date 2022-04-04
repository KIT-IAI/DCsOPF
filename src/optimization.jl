export checkGenerationCapacity,
        buildDCsOPF_core,
        addAllChanceConstraints!,
        balancing

function buildDCsOPF_core(ps::PowerSystem, solver; storage::Bool = true, glob::Bool = false)
    # extract variables
    T, Ψ = ps.horizon, ps.ptdf
    nu, nline = ps.inds[:N][:gen], ps.inds[:N][:line]
    ndunc, ndcert = ps.inds[:N][:dist][:unc], ps.inds[:N][:dist][:cert]
    umin, umax = ps.con[:val][:u][:min], ps.con[:val][:u][:max]
    Δumin, Δumax = ps.con[:val][:Δu][:min], ps.con[:val][:Δu][:max]
    plmin, plmax = ps.con[:val][:pl][:min], ps.con[:val][:pl][:max]
    h, H = ps.cost[:h], ps.cost[:H]
    certDist, uncDist = ps.dist[:cert], ps.dist[:unc]
    ns = storage ? ps.inds[:N][:stor] : 0

    ###################################################################################
    # storage variables: s, S
    mod = Model(solver)
    @variable(mod, u0[1:nu, 1:T])
    s0 = storage ? (@variable(mod, s0[1:ns, 1:T])) : similar(u0)
    if glob
        # global balancing
        @variable(mod, U[1:nu,i = 1:T,j = 1:i])
        S = storage ? @variable(mod, S[1:ns,i = 1:T,j = 1:i]) : similar(U)
    else
        # local balancing
        @variable(mod, U[1:nu,i = 1:T,j = 1:i,1:ndunc])
        S = storage ? @variable(mod, S[1:ns,i = 1:T,j = 1:i,1:ndunc]) : similar(U)
    end
    @variable(mod, γ >= 0)

    ###################################################################################
    # constraints: energy balance
    if storage
        [ @constraint(mod, sum(u0[:,k]) + sum(s0[:,k]) + sum(certDist[i][k] for i in 1:ndcert) + sum(uncDist[i][:μ][k] for i in 1:ndunc) == 0) for k in 1:T ]
    else
        [ @constraint(mod, sum(u0[:,k]) + sum(certDist[i][k] for i in 1:ndcert) + sum(uncDist[i][:μ][k] for i in 1:ndunc) == 0) for k in 1:T ]
    end
    if glob
        # global balancing
        if storage
            [ @constraint(mod, sum(U[i,k,l] for i in 1:nu) + sum(S[i,k,l] for i in 1:ns) + uncDist[n][:Σ][k,l] == 0) for k in 1:T, l in 1:T if k >= l for n in 1:ndunc]
        else
            [ @constraint(mod, sum(U[i,k,l] for i in 1:nu) + uncDist[n][:Σ][k,l] == 0) for k in 1:T, l in 1:T if k >= l for n in 1:ndunc]
        end
    else
        # local balancing
        if storage
            # with storage
            [ @constraint(mod, sum(U[i,k,l,n] for i in 1:nu) + sum(S[i,k,l,n] for i in 1:ns) + uncDist[n][:Σ][k,l] == 0) for k in 1:T, l in 1:T if k >= l for n in 1:ndunc]
        else
            # without storage
            [ @constraint(mod, sum(U[i,k,l,n] for i = 1:nu) + uncDist[n][:Σ][k,l] == 0) for k in 1:T, l in 1:T if k >= l for n in 1:ndunc]
        end
    end

    ###################################################################################
    # cost function constraint
    U_, h_, H_, Hsq, Hsqi = buildCostSOC(u0, U, diag(H), h, ndunc;glob = glob)
    @constraint(mod, [γ; Hsq .* U_ + Hsqi .* h_] in SecondOrderCone())
    # # @constraint(mod, norm(Hsq.*U_+Hsqi.*h_,2) <= γ)
    @objective(mod, Min, γ)
    fval(m) = objective_value(m)^2 - dot(h_, Hsqi .* Hsqi .* h_) + ps.cost[:c0]
    # fval(m) = m^2-dot(h_,Hsqi.*Hsqi.*h_)+ps.cost[:c0]

    ###################################################################################
    conrefs = Dict(:u => generate_undef_dict(nu, T), :Δu => generate_undef_dict(nu, T), :pl => generate_undef_dict(nline, T))

    if storage
        conrefs[:e] = generate_undef_dict(ns, T)
        conrefs[:s] = generate_undef_dict(ns, T)
    end

    OPF(mod, u0, U, s0, S, conrefs, fval, storage)
end

generate_undef_dict(n::Int, T::Int) = Dict(:max => build_undef_array(n, T), :min => build_undef_array(n, T))

build_undef_array(n::Int, T::Int) = Array{JuMP.ConstraintRef,2}(undef, n, T)

################################################
################################################
################################################
# needed only for comparisons
################################################
################################################
################################################
function buildDCOPF(ps::PowerSystem, solver)
    T, Ψ = ps.horizon, ps.ptdf
    nu, ndcert, ndunc, nline = ps.inds[:N][:gen], ps.inds[:N][:dist][:cert], ps.inds[:N][:dist][:unc], ps.inds[:N][:line]
    umin, umax = ps.con[:val][:u][:min], ps.con[:val][:u][:max]
    Δumin, Δumax = ps.con[:val][:Δu][:min], ps.con[:val][:Δu][:max]
    plmin, plmax = ps.con[:val][:pl][:min], ps.con[:val][:pl][:max]
    h, H = ps.cost[:h], ps.cost[:H]
    certDist, uncDist = ps.dist[:cert], ps.dist[:unc]
    ###################################################################
    m = Model(solver)
    @variable(m, u[1:nu,1:T])
    # energy balance
    [ @constraint(m, sum(u[:,k]) + sum(certDist[i][k] for i = 1:ndcert) + sum(uncDist[i][:μ][k] for i = 1:ndunc) == 0) for k = 1:T ]
    # inequalities
    @constraint(m, Umax[i = 1:nu,k = 1:T], u[i,k] <= umax[i,k])
    @constraint(m, Umin[i = 1:nu,k = 1:T], u[i,k] >= umin[i,k])
    @constraint(m, ΔUmax[i = 1:nu,k = 2:T], u[i,k] - u[i,k - 1] <= Δumax[i,k])
    @constraint(m, ΔUmin[i = 1:nu,k = 2:T], u[i,k] - u[i,k - 1] >= Δumin[i,k])
    @constraint(m, Plmin[l = 1:nline,k = 1:T], flow_pl(l, k, u, ps) >= plmin[l,k])
    @constraint(m, Plmax[l = 1:nline,k = 1:T], flow_pl(l, k, u, ps) <= plmax[l,k])
    # [ @constraint(m, umin[i,k]<=u[i,k]<=umax[i,k]) for i=1:nu, k=1:T ]
    # [ @constraint(m, Δumin[i,k]<=u[i,k]-u[i,k-1]<=Δumax[i,k]) for i=1:nu, k=2:T ]
    # [ @constraint(m, plmin[l,k]<=flow_pl(l,k,u,ps)<=plmax[l,k]) for l=1:nline, k=1:T ]
    @objective(m, Min, sum(dot(u[:,k], H * u[:,k]) + dot(h, u[:,k]) for k = 1:T) + ps.cost[:c0])
    return m
end

# function balancing(mysys::PowerSystem, optimizer; glob::Bool=false, storage::Bool=false, variance::Bool=false, text::String="Start balancing")
#     # optimizer_with_attributes(Mosek.Optimizer,  "Threads" => 2)
#     # options: Mosek.MSK_IINF_INTPNT_NUM_THREADS (interior-point algorithm), Mosek.MSK_IPAR_NUM_THREADS (optimizer)
#     opf = buildDCsOPF_core(mysys, optimizer; storage = storage, glob = glob)
#     if storage && variance
#         addSpecificConstraints_Storage_WithVariance!(mysys, opf; glob = glob)
#     elseif storage
#         addSpecificConstraints_Storage!(mysys, opf; glob = glob)
#     else
#         addSpecificConstraints_NoStorage!(mysys, opf; glob = glob)
#     end
#     addAllChanceConstraints!(mysys, opf; glob = glob)
#     optimize!(opf.m)
#     return opf
# end

function addAllChanceConstraints!(ps::PowerSystem, opf::OPF;glob::Bool = false)
    T = ps.horizon
    inds = ps.inds
    nu, nl, ns = ps.inds[:N][:gen], ps.inds[:N][:line], ps.inds[:N][:stor]
    bnds = [:lb :ub]
    storage = opf.stor
    for t in 1:T
        # generation constraints u
        for i in 1:nu
            addCC!(opf, ps, i, t, :u, :lb, glob = glob)
            addCC!(opf, ps, i, t, :u, :ub, glob = glob)
        end
        # Δ-generation constraints Δu
        if t >= 2
            for i in 1:nu
                addCC!(opf, ps, i, t, :Δu, :lb, glob = glob)
                addCC!(opf, ps, i, t, :Δu, :ub, glob = glob)
            end
        end
        # if storages are present, then add storage and energy constraints
        if storage
            # storage constraints s
            for i in 1:ns
                addCC!(opf, ps, i, t, :s, :lb, glob = glob)
                addCC!(opf, ps, i, t, :s, :ub, glob = glob)
            end
            # energy constraints e
            for i in 1:ns
                addCC!(opf, ps, i, t, :e, :lb, glob = glob)
                addCC!(opf, ps, i, t, :e, :ub, glob = glob)
            end
        end
        # line flows pl
        for i in 1:nl
            addCC!(opf, ps, i, t, :pl, :lb, storage = storage, glob = glob)
            addCC!(opf, ps, i, t, :pl, :ub, storage = storage, glob = glob)
        end

    end
end

function addCC!(opf::OPF, s::PowerSystem, i::Int64, t::Int64, x::Symbol, kind::Symbol;eval::Bool = false,storage::Bool = true,glob::Bool = false)
    @assert kind in [:lb, :ub]
    @assert x in [:u, :Δu, :pl, :e, :s] "variable $x not known"
    bnd = kind == :ub ? :max : :min
    ref, val, λ = opf.con[x][bnd], s.con[:val][x][bnd], s.con[:λ][x][bnd]

    # compute expected value and standard deviation
    ev, sig = 
    if x == :u
        moments_u(i, t, opf, s;eval = eval,glob = glob)
    elseif x == :Δu
        moments_Δu(i, t, opf, s;eval = eval,glob = glob)
    elseif x == :pl
        # implementation depends on presence of storages
        moments_pl(i, t, opf, s;eval = eval,storage = storage,glob = glob)
    elseif x == :e
        moments_e(i, t + 1, opf, s;eval = eval,glob = glob)
    elseif x == :s
        moments_s(i, t, opf, s;eval = eval,glob = glob)
    else
        throw(error("variable $x not known"))
    end

    # add constraint upper/lower bound
    if kind == :ub
        if eval
            # return ev+λ[i,t]*sqrt(v)-val[i,t]
            throw(error("not yet implemented"))
            return ev + λ[i,t] * sig - val[i,t]
        else
            ref[i,t] = @constraint(opf.m, [val[i,t] - ev; λ[i,t] * sig] in SecondOrderCone())
            # ref[i,t] = @constraint(opf.m, ev+λ[i,t]*sig<=val[i,t])
        end
    else
        if eval
            throw(error("not yet implemented"))
            return -(ev - λ[i,t] * sig - val[i,t])
        else
            ref[i,t] = @constraint(opf.m, [ev - val[i,t]; λ[i,t] * sig ] in SecondOrderCone())
            # ref[i,t] = @constraint(opf.m, ev-λ[i,t]*sig>=val[i,t])
        end
    end
end