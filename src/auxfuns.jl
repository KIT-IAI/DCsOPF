export buildInd,
       initializeDisturbance,
       plot_costs

function _initializeCertainDistBuses(d::Dict)
    [d["load"][string(i)]["load_bus"] for i = 1:length(d["load"]) ]
end

function _initializeGenerators(d::Dict)
    [d["gen"][string(i)]["gen_bus"] for i = 1:length(d["gen"]) ] # gen_bus --> not bus?? (see casefile)
end

function _buildMap(buses::Vector{Int}, N::Int)
    """
        Save subset of buses in empty dict of length N (=all buses), 
        rest is empty, and assign new indices.
    """
    # @assert length(unique(buses))==length(buses) "Multiple sources not yet supported"
    
    #numbering = collect(1:length(buses)) # needed what for??
    # k = 1 # needed what for??
    d = Dict()
    [ d[i] = [] for i = 1:N ]
    for (a, b) in enumerate(buses)
        # b is original index in set of all buses (length N)
        # a assigns new indices to subset 'buses::Vector{Int}'
        push!(d[b], a) 
    end
    d
end

function buildInd(data::Dict, dists_unc::Union{Vector{Int},Vector{Any}}, stors::Union{Vector{Int},Vector{Any}})
    
    N = length(data["bus"])
    gens = _initializeGenerators(data)
    dict_gen = Dict(:buses => gens, :at_bus => _buildMap(gens, N))
    
    # ensure that a bus does not have a certain + uncertain part
    length(dists_unc) == 0 ? dists_cert = _initializeCertainDistBuses(data) : dists_cert = setdiff(_initializeCertainDistBuses(data), dists_unc)
    # indices for certain disturbances
    dict_distcert = Dict(:buses => dists_cert, :at_bus => _buildMap(dists_cert, N))
    # indices for uncertain disturbances
    dict_distunc = 
    if length(dists_unc) == 0
        Dict(:buses => [], :at_bus => [])
    else
        Dict(:buses => dists_unc, :at_bus => _buildMap(dists_unc, N))
    end
    
    # indices for storages
    dict_stor = 
    if length(stors) == 0
        Dict(:buses => [], :at_bus => [])
    else
        Dict(:buses => stors, :at_bus => _buildMap(stors, N))
    end
    
    # all in dict
    dict_N = Dict(:bus => N,
                  :gen => length(gens),
                  :dist => Dict(:cert => length(dists_cert),
                              :unc => length(dists_unc)),
                  :stor => length(stors),
                  :line => length(data["branch"]))
    
    return Dict(:gen => dict_gen,
                :dist => Dict(:cert => dict_distcert, :unc => dict_distunc),
                :stor => dict_stor,
                :N => dict_N)
end

buildInd(d::Dict) = buildInd(d, [], [])

# take nominal value and stretch over horizon T
function initializeDisturbance(d::Dict, T::Int, ind::Dict, kind;α::Real = 0.1)
    @assert kind in [:cert, :unc] "kind $kind not valid"
    buses = ind[:dist][kind][:buses]
    p, q = _getDemand(d)
    load = Dict()
    for (i, bus) in enumerate(buses)
        load[i] = 
        if kind == :cert
            p[bus] * ones(T)
        else
            Dict(:μ => p[bus] * ones(T), :Σ => α * p[bus] / 3 * diagm(ones(T)))
        end
    end
    load
end

function buildL(X, m::Int, nd::Int, N::Int)
    Y = zeros(N, N)
    for i = 1:N
        for j = 1:i
            Y[i,j] = getvalue(X[m,i,j,nd])
        end
    end
    Y
end

function checkGenerationCapacity(ps::PowerSystem)
    width = 0.35
    ind = ps.inds
    T = ps.horizon
    tspan = 1:T
    p_g = [ sum(ps.con[:val][:u][:max][:,k]) for k = 1:T ]
    lst_distunc = zeros(T)
    lst_distcert = zeros(T)
    for k = 1:T
        # for testing
        # for i = 1:ind[:N][:dist][:unc]
        #     display(ps.dist[:unc][i][:μ][k])
        # end
        # for i = 1:ind[:N][:dist][:cert]
        #     display(ps.dist[:cert][i][k])
        # end
        l_unc = [ps.dist[:unc][i][:μ][k] for i = 1:ind[:N][:dist][:unc]]
        l_cert = [ps.dist[:cert][i][k] for i = 1:ind[:N][:dist][:cert]]
        if !isempty(l_unc)
            lst_distunc[k] = sum(l_unc)
        end
        if !isempty(l_cert)
            lst_distcert[k] = sum(l_cert)
        end
    end
    # p_d_mean = abs.([(isempty(lst_distunc[k]) ? 0 : sum(lst_distunc[k])) + (isempty(lst_distcert[k]) ? 0 : sum(lst_distcert[k])) for k = 1:T])
    # OLD
    # display(lst_distcert)
    # display(lst_distunc)
    if ! (lst_distcert[1] == 0.0)
        p_d_mean = [ -sum(ps.dist[:unc][i][:μ][k] for i = 1:ind[:N][:dist][:unc]) - sum(ps.dist[:cert][i][k] for i = 1:ind[:N][:dist][:cert]) for k = 1:T ]
    else
        p_d_mean = [ -sum(ps.dist[:unc][i][:μ][k] for i = 1:ind[:N][:dist][:unc]) for k = 1:T ]
    end
    # display("pg")
    # display(p_g)
    # display("pd mean")
    # display(p_d_mean)
    enough_generation = (p_g - p_d_mean .>= 0)
    # display("enough generation")
    # display(enough_generation)
    res = (p_g - p_d_mean) ./ p_g
    min_res, ind_min_res = findmin(res)
    res_color = 
    if min_res < 0.1
        "red"
    elseif 0.1 <= min_res <= 0.20
        "yellow"
    else
        "green"
    end

    figure(444)
    bar(tspan[ind_min_res], 1.05 * p_g[ind_min_res], width = 2.5 * width, color = res_color, alpha = 0.45)
    grid(true)
    bar(tspan .- width / 2, p_g, width = width, label = "Total generation")

    bar(tspan .+ width / 2, p_d_mean, width = width, label = "Total demand")
    legend()
    xlabel("Time"); ylabel("Active power")

    length(findall(x->x == false, enough_generation)) > 0 && throw(error("Not enough generation capacity"))

    return "The minimum capacity reserve is $(min_res * 100) % at time t=$(tspan[ind_min_res])."
end

function plot_costs(d::Dict;fignum::Int64 = 999,figsize::Tuple = (3, 6), layout_pad = 0)
    width = 0.35
    tspan = 1:3
    cmax = maximum(d[:global])

    figure(num = fignum, figsize = figsize)
    
    bar(tspan - width / 2, d[:local], width = width, label = "Local")
    bar(tspan + width / 2, d[:global], width = width, label = "Global")
    legend()
    ylabel(L"Absolute cost")
    ylim(0.95 * maximum(d[:global]), 584000)
    tight_layout(pad = layout_pad)
    xticks([1,2,3])

    figure(num = fignum + 1, figsize = figsize)
    
    bar(tspan - width / 2, 100 * d[:local] / cmax, width = width, label = "Local")
    bar(tspan + width / 2, 100 * d[:global] / cmax, width = width, label = "Global")
    # legend()
    ylabel(L"Relative cost")
    ylim(95, 100.5)
    ax = gca()
    # [ text(tspan[i]-width,100*d[:local][i]/cmax,"$(100*round(1000*d[:local][i]/cmax)/1000)") for i=1:3 ]
    # [ text(tspan[i],100*d[:global][i]/cmax,"$(100*round(1000*d[:global][i]/cmax)/1000)") for i=2:3 ]
    xticks([1,2,3])
    tight_layout(pad = layout_pad)

end

function flow_pl(l::Int64, t::Int64, u, s::PowerSystem)
    Nbus = s.inds[:N][:bus]
    gen, certDist, uncDist = s.inds[:gen], s.inds[:dist][:cert], s.inds[:dist][:unc]
    Nunc = s.inds[:N][:dist][:unc]
    EV = JuMP.GenericAffExpr{Float64,JuMP.VariableRef}(0)
    #################################################################
    # add certain disturbance
    for i in certDist[:buses]
        @assert length(certDist[:at_bus][i]) == 1 "Multiple sources not supported"
        no = certDist[:at_bus][i][1]
        EV += s.ptdf[l,i] * s.dist[:cert][no][t]
    end
    # add uncertain disturbance
    for i in uncDist[:buses]
        @assert length(uncDist[:at_bus][i]) == 1 "Multiple sources not supported"
        no = uncDist[:at_bus][i][1]
        EV += s.ptdf[l,i] * s.dist[:unc][no][:μ][t]
    end
    # add generation
    for i in gen[:buses]
        @assert length(gen[:at_bus][i]) == 1 "Multiple sources not supported"
        no = gen[:at_bus][i][1]
        EV += s.ptdf[l,i] * u[no,t]
    end
    #################################################################
    return EV
end