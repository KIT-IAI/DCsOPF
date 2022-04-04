export moments,
        moments_e,
        moments_pl,
        moments_pl_WithoutStorage,
        moments_pl_WithStorage,
        moments_s,
        moments_u,
        moments_Δu

function moments(i::Int64, t::Int64, opf::OPF, ps::PowerSystem;eval::Bool = false,kind::Symbol = :undef,glob::Bool = false)
    if kind == :u
        return moments_u(i, t, opf, ps;eval = eval,glob = glob)
    elseif kind == :Δu
        return moments_Δu(i, t, opf, ps;eval = eval,glob = glob)
    elseif kind == :pl
        return moments_pl(i, t, opf, ps;eval = eval,storage = opf.stor,glob = glob)
    elseif kind == :e
        return moments_e(i, t, opf, ps;eval = eval,glob = glob)
    elseif kind == :s
        return moments_s(i, t, opf, ps;eval = eval,glob = glob)
    end
end
#################################################################
#################################################################
# Generation
#################################################################
#################################################################
function moments_u(i::Int64, t::Int64, opf::OPF, ps::PowerSystem;eval::Bool = false,glob::Bool = false)
    u, U =
    if eval
        value.(opf.u), value.(opf.U)
    else
        opf.u, opf.U
    end
    moments_u(i, t, u, U, ps; glob = glob)
end
# core function
function moments_u(i::Int64, t::Int64, u, U, ps::PowerSystem;glob::Bool = false)
    N = ps.inds[:N][:dist][:unc]
    EV = u[i,t]
    VAR =
    if glob
        # global balancing
        [ sqrt(N) * U[i,t,k] for k = 1:t ]
    else
        # local balancing
        [ U[i,t,k,j] for k = 1:t for j = 1:N ]
    end
    EV, VAR
end

#################################################################
#################################################################
## Storage
#################################################################
#################################################################
function moments_s(i::Int64, t::Int64, opf::OPF, ps::PowerSystem;eval::Bool = false,glob::Bool = false)
    s, S = 
    if eval
        value.(opf.s), value.(opf.S)
    else
        opf.s, opf.S
    end
    moments_u(i, t, s, S, ps; glob = glob)
end
#################################################################
#################################################################
## Δ-Generation
#################################################################
#################################################################
function moments_Δu(i::Int64, τ::Int64, opf::OPF, ps::PowerSystem;eval::Bool = false,glob::Bool = false)
    @assert τ >= 2 "Constraints for Δu only valid for times>=2"
    u, U =
    if eval
        value.(opf.u), value.(opf.U)
    else
        opf.u, opf.U
    end
    moments_Δu(i, τ, u, U, ps; glob = glob)
end
# core function
function moments_Δu(i::Int64, τ::Int64, u, U, ps::PowerSystem;glob::Bool = false)
    N = ps.inds[:N][:dist][:unc]
    EV  = u[i,τ] - u[i,τ - 1]
    VAR =
    if glob
        # global balancing
        [ sqrt(N) * U[i,τ,τ];
                sqrt(N) * [ U[i,τ,k] - U[i,τ - 1,k] for k = 1:τ - 1 ]]
    else
        # local balancing
        [ [ U[i,τ,τ,j] for j = 1:N ];
                [ U[i,τ,k,j] - U[i,τ - 1,k,j] for k = 1:τ - 1 for j = 1:N ]]
    end
    EV, VAR
end
#################################################################
#################################################################
## Energy
#################################################################
#################################################################
function moments_e(i::Int64, t::Int64, pf::OPF, ps::PowerSystem;eval::Bool = false,glob::Bool = false)
    s, S =
    if eval
        value.(pf.s), value.(pf.S)
    else
        pf.s, pf.S
    end
    moments_e(i, t, s, S, ps;glob = glob)
end
# core function
function moments_e(i::Int64, t::Int64, s, S, ps::PowerSystem; glob::Bool = false)
    @assert t >= 2
    EVic, VARic, h = ps.stor[i][:ic][:μ], ps.stor[i][:ic][:σ²], ps.stor[i][:h]
    @assert abs(VARic) <= 1e-10 "Storages with initial uncertainty are not supported--double-check the formula."
    N = ps.inds[:N][:dist][:unc]
    # expected value
    EV = EVic - h * sum(s[i,k] for k in 1:t - 1)
    # variance (or rather: standard deviation (=sqrt(VAR)))
    VAR =
    if glob
        # global balancing (all j=1:N are the same, hence take 1 and factor sqrt(N))
        -[ h * sqrt(N) * sum(S[i,l,k] for l in k:t - 1) for k in 1:t - 1 ]
    else
        # local balancing (take all j=1:N for themselves)
        -[ h * sum(S[i,l,k,j] for l in k:t - 1) for j in 1:N for k in 1:t - 1 ]
    end
    EV, [VARic; VAR]
end
#################################################################
#################################################################
# Line flows
#################################################################
#################################################################
function moments_pl(l::Int64, t::Int64, opf::OPF, ps::PowerSystem;eval::Bool = false,storage::Bool = true,glob::Bool = false)
    if eval
        u, U = value.(opf.u), value.(opf.U)
        if storage
            s, S = value.(opf.s), value.(opf.S)
        end
    else
        u, U = opf.u, opf.U
        if storage
            s, S = opf.s, opf.S
        end
    end
    storage ? moments_pl_WithStorage(l, t, u, U, s, S, ps;eval = eval,glob = glob) : moments_pl_WithoutStorage(l, t, u, U, ps;eval = eval,glob = glob)
end
# core function
function moments_pl_WithoutStorage(l::Int64, t::Int64, u, U, ps::PowerSystem;eval::Bool = false,glob::Bool = false)
    Nbus = ps.inds[:N][:bus]
    gen, certDist, uncDist = ps.inds[:gen], ps.inds[:dist][:cert], ps.inds[:dist][:unc]
    Nunc = ps.inds[:N][:dist][:unc]
    EV, VAR = 
    if eval
        0., zeros(Nunc, t)
    else
        JuMP.GenericAffExpr{Float64,JuMP.VariableRef}(0), Array{JuMP.GenericAffExpr{Float64,JuMP.VariableRef},2}(undef, Nunc, t)
    end
    #################################################################
    # EV
    # add certain disturbance
    for i in certDist[:buses]
        @assert length(certDist[:at_bus][i]) == 1 "Multiple sources not supported"
        no = certDist[:at_bus][i][1]
        EV += ps.ptdf[l,i] * ps.dist[:cert][no][t]
    end
    # add uncertain disturbance
    for i in uncDist[:buses]
        @assert length(uncDist[:at_bus][i]) == 1 "Multiple sources not supported"
        no = uncDist[:at_bus][i][1]
        EV += ps.ptdf[l,i] * ps.dist[:unc][no][:μ][t]
    end
    # add generation
    for i in gen[:buses]
        @assert length(gen[:at_bus][i]) == 1 "Multiple sources not supported"
        no = gen[:at_bus][i][1]
        EV += ps.ptdf[l,i] * u[no,t]
    end
    #################################################################
    # VAR
    for k in 1:t
        # sum only over those buses that have an uncertain disturbance
        for (nξ, i) in enumerate(uncDist[:buses])
            
            # add uncertain disturbance
            no = uncDist[:at_bus][i][1]
            # @assert no==nξ "something is fishy with the numbering of the uncertain disturbances"
            VAR[nξ,k] = ps.ptdf[l,i] * ps.dist[:unc][no][:Σ][t,k]
            
            # add generation
            for j in gen[:buses]
                no_g = gen[:at_bus][j][1]
                if glob
                    # global balancing
                    VAR[nξ,k] += ps.ptdf[l,j] * U[no_g,t,k]
                else
                    # local balancing
                    VAR[nξ,k] += ps.ptdf[l,j] * U[no_g,t,k,no]
                end
            end
        end
    end

    EV, reshape(VAR, Nunc * t)
end

function moments_pl_WithStorage(l::Int64, t::Int64, u, U, s, S, ps::PowerSystem;eval::Bool = false,glob::Bool = false)
    Nunc = ps.inds[:N][:dist][:unc]
    uncDist = ps.inds[:dist][:unc]
    stor = ps.inds[:stor]
    # get EV, VAR for line flow without storage
    EV, VAR_ = moments_pl_WithoutStorage(l, t, u, U, ps, eval = eval, glob = glob)
    # re-arrange variance from vector to matrix
    VAR = reshape(VAR_, Nunc, t)
    @assert size(VAR) == (Nunc, t) "dimensions do not agree"

    # add storage
    for i in stor[:buses]
        @assert length(stor[:at_bus][i]) == 1 "Multiple sources not supported"
        no = stor[:at_bus][i][1]
        EV += ps.ptdf[l,i] * s[no,t]
    end

    # add storage
    for k in 1:t
        # sum only over those buses that have an uncertain disturbance
        for (nξ, i) in enumerate(uncDist[:buses])
            no = uncDist[:at_bus][i][1]
            # add storage
            for j in stor[:buses]
                no_s = stor[:at_bus][j][1]
                if glob
                    # global balancing
                    VAR[nξ,k] += ps.ptdf[l,j] * S[no_s,t,k]
                else
                    # local balancing
                    VAR[nξ,k] += ps.ptdf[l,j] * S[no_s,t,k,no]
                end
            end
        end
    end
    EV, reshape(VAR, Nunc * t)
end
