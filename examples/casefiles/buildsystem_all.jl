export buildSystem,
        addSpecificConstraints_NoStorage!,
        addSpecificConstraints_Storage!,
        addSpecificConstraints_Storage_WithVariance!

function buildSystem(name::String, uncertainties::Vector{Int64}, storages::Vector{Int64}, params::CaseParams; ε::Float64=0.025)
    data = PowerModels.parse_file("casefiles/case" * name * ".m")
    T = 12

    PTDF = matread("casefiles/case" * name * "_PTDF.mat")["ptdf"]

    Nbus = length(data["bus"]);
    # specify buses with uncertain disturbance
    uncertainDisturbanceBuses = sort(uncertainties)
    # specify buses with storage
    storageBuses = sort(storages)
    # generate indices
    inds = buildInd(data, uncertainDisturbanceBuses, storageBuses)
    # specify uncertainty
    certDist, uncDist = initializeUncertainty(data, inds, T, params)
    ## Define cost & constraints
    nu, ns, nline = inds[:N][:gen], inds[:N][:stor], inds[:N][:line]
    #########################################
    # Cost Coefficients
    # The linear cost has to be multiplied by 2
    # Why? For the SOC formulation we assume a formulation
    # x'Hx + 0.5*h'x + c0
    H, ~, h, c0 = PowerModelsParsing._getCost(data)
    h *= 2
    #########################################
    # Constraints
    umin_, umax_, plmin_, plmax_, ~, ~ = PowerModelsParsing._getConstraints(data)
    umin, umax = repeat(umin_, 1, T), 1.1 * repeat(umax_, 1, T)
    γ = 0.15 # we allow a maximum change of ±γ % for the generators
    Δumax_, Δumin_ = γ * umax_, -γ * umax_
    Δumin, Δumax = repeat(Δumin_, 1, T), repeat(Δumax_, 1, T)
    emin, emax = zeros(ns, T), params.stor_ub * ones(ns, T)
    # emax[2,:] .= 1.5
    smin, smax = -10 * ones(ns, T), 10 * ones(ns, T)
    # modify line flow limits
    # plmin_[[17,24]], plmax_[[17,24]] = -[1.5,1.5], [1.5,1.5]

    plmin, plmax = repeat(plmin_, 1, T), repeat(plmax_, 1, T)
    λ = norminvcdf(1 - ε)

    # storage specification
    stor = Dict{Int64,Dict}()
    for i = 1:length(storageBuses)
        stor[i] = Dict(:h => 1, :ic => Dict(:μ => 2, :σ² => 0.0))
    end
    # stor[2][:ic][:μ] = 1.5

    ## build power system
    constraints = Dict(:val => Dict(:u => _generate_dict(umin, umax),
                                    :Δu => _generate_dict(-Δumax, Δumax),
                                    :pl => _generate_dict(plmin, plmax),
                                    :e => _generate_dict(emin, emax),
                                    :s => _generate_dict(smin, smax)),
                       :λ => Dict(:u => _generate_dict(λ, nu, T),
                                  :Δu => _generate_dict(λ, nu, T),
                                  :pl => _generate_dict(λ, nline, T),
                                  :e => _generate_dict(λ, ns, T),
                                  :s => _generate_dict(λ, ns, T))
                  )
    ps = PowerSystem(T, inds, Dict(:cert => certDist, :unc => uncDist), stor, PTDF, constraints, Dict(:h => h, :H => H, :c0 => c0))

    return ps, data
end

_generate_dict(lb, ub) = Dict(:min => lb, :max => ub)
_generate_dict(bnd) = _generate_dict(bnd, bnd)
_generate_dict(λ, n, T) = _generate_dict(λ * ones(n, T))

function my_mean(t, p, γ, h, a = 0)
    p * (1 + γ * sin(2 * pi * (t - 1 + a) / h))
end

function extract_mean(t, file)
    read(file, "mu_post")[t]
end

function initializeUncertainty(case::Dict, ind::Dict, T, params::CaseParams; figsize::Tuple = figsize, linewidth = linewidth, fontdict = fontdict)
    γ = params.γ # change in amplitude for sinusoidal mean
    tspan = 1:T
    ref = PowerModels.build_ref(case)[:it][:pm][:nw][0]
    
    certain = initializeDisturbance(case, T, ind, :cert)
    inds_certain = ind[:dist][:cert][:buses]
    N_certain = length(inds_certain)

    @show shift_certain = zeros(N_certain)
    # bus_shift_certain = sort([7,9])
    # @assert issubset(bus_shift_certain,ind[:dist][:cert][:buses]) "Some buses from $bus_shift_certain are NOT certain disturbances."
    # @show ind_shift = [ ind[:dist][:cert][:at_bus][i][1] for i in bus_shift_certain ]
    # shift_certain[ind_shift] = -T/2

    uncertain = initializeDisturbance(case, T, ind, :unc)
    inds_unc = ind[:dist][:unc][:buses]
    @show N_unc = ind[:N][:dist][:unc]
    # v = unc_factor * ones(N_unc)

    shift_uncertain = zeros(N_unc) # N_certain?
    # bus_shift_uncertain = sort([4,8])
    # @assert issubset(bus_shift_uncertain,ind[:dist][:unc][:buses]) "Some buses of $bus_shift_uncertain are NOT uncertain disturbances."
    # @show ind_shift = [ ind[:dist][:unc][:at_bus][i][1] for i in bus_shift_uncertain ]
    # shift_uncertain[ind_shift] = -T/2
    # shift_uncertain[setdiff(1:N_unc,ind_shift)] = T/2

    # assign and plot certain loads
    for (i, bus) in enumerate(inds_certain)
        file = matopen("data/CovMatrix_artificial.mat")
        pnom = case["load"][string(ref[:bus_loads][bus][1])]["pd"]
        μ1(t) = my_mean(t, pnom, params.γ, 24, shift_certain[i])
        mu = μ1.(tspan) * params.load_factor

        (mu_, sign) = params.artificial_unc ? (mu, "-") : (-mu, "")
        certain[i] = -mu

        plot_certain(mu_, T, bus, sign, figsize = figsize, linewidth = linewidth, fontdict = fontdict)
    end
    
    # assign and plot uncertain loads
    for (i, bus) in enumerate(inds_unc)

        mu, L, V = [], [], []
        if params.artificial_unc
            file = matopen("data/CovMatrix_artificial.mat")
            pnom = case["load"][string(ref[:bus_loads][bus][1])]["pd"] # TODO bus or i??
            μ1(t) = my_mean(t, pnom, 0.4, 24, shift_uncertain[i])
            L1 = params.deterministic ? zeros(T,T) : read(file, "Lpost")[1:T,1:T] *  params.unc_factor
            V1 = L1 * L1'
            mu, L, V = μ1.(tspan), L1, V1
        else
            file = matopen("data/CovMatrix_"*string(i)*".mat")
            μ2(t) = extract_mean(t, file)
            L2 = params.deterministic ? zeros(T,T) : params.unc_factor * read(file, "Lpost")[1:T,1:T]
            V2 = params.deterministic ? zeros(T,T) : params.unc_factor * params.unc_factor * read(file, "cov")[1:T,1:T]
            mu, L, V = μ2.(tspan), L2, V2
        end

        mu = mu * params.dist_factor

        if params.artificial_unc
            mu_, L_, V_, sign = -mu, -L, V, "-"
        else # uncertainties are generation
            mu_, L_, V_, sign  = mu, L, V, ""
        end

        uncertain[i][:μ], uncertain[i][:Σ], uncertain[i][:Σ_full] =  mu_, L_, V_ 

        plot_uncertain( mu_, L_, V_ , sign, T, bus, figsize = figsize, linewidth = linewidth, fontdict = fontdict)
    end

    return certain, uncertain
end

function addSpecificConstraints_NoStorage!(ps::PowerSystem, opf::OPF; glob::Bool = false)
    return opf
end

function addSpecificConstraints_Storage!(ps::PowerSystem, opf::OPF; glob::Bool = false)
    T = ps.horizon
    nu, ns = ps.inds[:N][:gen], ps.inds[:N][:stor]
    emax = ps.con[:val][:e][:max]
    umax = ps.con[:val][:u][:max]
    # terminal constraint for storages
    z, risklevel = 0.05, 0.05

    @constraint(opf.m, Eendmax[i = 1:ns], [(1 + z) * ps.stor[i][:ic][:μ] - moments_e(i, T, opf.s, opf.S, ps;glob = glob)[1]; norminvcdf(1 - risklevel) * moments_e(i, T, opf.s, opf.S, ps;glob = glob)[2]] in SecondOrderCone())
    @constraint(opf.m, Eendmin[i = 1:ns], [moments_e(i, T, opf.s, opf.S, ps;glob = glob)[1] - (1 - z) * ps.stor[i][:ic][:μ]; norminvcdf(1 - risklevel) * moments_e(i, T, opf.s, opf.S, ps;glob = glob)[2]] in SecondOrderCone())

    # @constraint(opf.m, Eendmax[i=1:ns],                              moments_e(i,T,opf.s,opf.S,ps;glob=glob)[1] + norminvcdf(1-risklevel)*moments_e(i,T,opf.s,opf.S,ps;glob=glob)[2]<= (1+z)*ps.stor[i][:ic][:μ])
    # @constraint(opf.m, Eendmin[i=1:ns], (1-z)*ps.stor[i][:ic][:μ] <= moments_e(i,T,opf.s,opf.S,ps;glob=glob)[1] - norminvcdf(1-risklevel)*moments_e(i,T,opf.s,opf.S,ps;glob=glob)[2]                            )

    # variance constraint for generator
    # [ @constraint(opf.m, moments_u(i,k,opf.u,opf.U,ps;glob=glob)[2] <= 0.01) for i in [2,4] for k in union(1:5,8:12) ]
    return opf
end

function addSpecificConstraints_Storage_WithVariance!(ps::PowerSystem, opf::OPF; glob::Bool = false)
    T = ps.horizon
    nu, ns = ps.inds[:N][:gen], ps.inds[:N][:stor]
    emax = ps.con[:val][:e][:max]
    umax = ps.con[:val][:u][:max]
    # terminal constraint for storages
    z, risklevel = 0.05, 0.05
    @constraint(opf.m, Eendmax[i = 1:ns], [ (1 + z) * ps.stor[i][:ic][:μ] - moments_e(i, T, opf.s, opf.S, ps;glob = glob)[1]; norminvcdf(1 - risklevel) * moments_e(i, T, opf.s, opf.S, ps;glob = glob)[2]] in SecondOrderCone())
    @constraint(opf.m, Eendmin[i = 1:ns], [ moments_e(i, T, opf.s, opf.S, ps;glob = glob)[1] - (1 - z) * ps.stor[i][:ic][:μ]; norminvcdf(1 - risklevel) * moments_e(i, T, opf.s, opf.S, ps;glob = glob)[2]] in SecondOrderCone())

    # @constraint(opf.m, Eendmax[i=1:ns],                              moments_e(i,T,opf.s,opf.S,ps;glob=glob)[1] + norminvcdf(1-risklevel)*moments_e(i,T,opf.s,opf.S,ps;glob=glob)[2]<= (1+z)*ps.stor[i][:ic][:μ])
    # @constraint(opf.m, Eendmin[i=1:ns], (1-z)*ps.stor[i][:ic][:μ] <= moments_e(i,T,opf.s,opf.S,ps;glob=glob)[1] - norminvcdf(1-risklevel)*moments_e(i,T,opf.s,opf.S,ps;glob=glob)[2]                            )
    
    # variance constraint for generator
    nr_generators = length(mysys.inds[:gen][:buses])
    [ @constraint(opf.m, [0.01 ; moments_u(i, k, opf.u, opf.U, ps;glob = glob)[2]] in SecondOrderCone()) for i in 1:nr_generators for k in collect(4:9) ]
    # [ @constraint(opf.m, moments_u(i,k,opf.u,opf.U,ps;glob=glob)[2] <= 0.01) for i in [2,4] for k in union(1:5,8:12) ]
    return opf
end

function plot_certain(mu::AbstractVector, T::Int, bus::Int, sign::String; col::String = "#228B22", figsize::Tuple = (2, 2), linewidth = 1, fontdict = Dict("color" => "black", "size" => fontsize, "family" => "serif"), layout_pad = 0, labelsize = 8)
    tspan = 1:T

    PyPlot.figure(num = bus, figsize = figsize)
    ax = gca()

    majorformatter = matplotlib.ticker.FormatStrFormatter("%.1f")
    ax.yaxis.set_major_formatter(majorformatter)
    PyPlot.grid(true)

    PyPlot.plot(tspan, mu, color = col, linewidth = linewidth)

    PyPlot.xlabel(latexstring("t"), fontdict = fontdict)
    PyPlot.ylabel(latexstring(sign,"d_{$bus}(t)"), fontdict = fontdict)
    # title("Uncertain disturbance at bus $bus")
    mx = matplotlib.ticker.MultipleLocator(2) # Define interval of minor ticks
    ax.xaxis.set_minor_locator(mx) # Set interval of minor ticks
    PyPlot.xticks(1:2:T)
    PyPlot.setp(ax.get_yticklabels(), size = labelsize)
    PyPlot.setp(ax.get_xticklabels(), size = labelsize)
    PyPlot.axis("tight")
    PyPlot.xlim(minimum(tspan), maximum(tspan))
    # grid(b=true,which="minor")
    bus == 15 && PyPlot.yticks([3.2,3.3,3.4,3.5])
    PyPlot.tight_layout(pad = layout_pad)
end

function plot_uncertain(μ::AbstractVector, L::AbstractMatrix, V::AbstractMatrix, sign::String, T::Int, bus::Int; col::String = "#228B22", figsize::Tuple = (2, 2), linewidth = 1, fontdict = Dict("color" => "black", "size" => fontsize, "family" => "serif"), layout_pad = 0, labelsize = 8, Nscens::Int = 10)
    tspan = 1:T
    scens = repeat(μ, 1, Nscens) + L * randn(T, Nscens)
    V = L * L'

    PyPlot.figure(num = bus, figsize = figsize)
    ax = gca()

    majorformatter = matplotlib.ticker.FormatStrFormatter("%.1f")
    ax.yaxis.set_major_formatter(majorformatter)
    PyPlot.grid(true)  

    PyPlot.plot(tspan, scens, linestyle = "dotted", color = col, linewidth = linewidth)
    PyPlot.fill_between(tspan, μ - 3 * sqrt.(diag(V)), μ + 3 * sqrt.(diag(V)), color = col, alpha = 0.25, linewidth = 0)
    PyPlot.plot(tspan, μ, linestyle = "solid", color = col, linewidth = linewidth)

    PyPlot.xlabel(latexstring("t"), fontdict = fontdict)
    PyPlot.ylabel(latexstring(sign*"d_{$bus}(t)"), fontdict = fontdict)
    # title("Uncertain disturbance at bus $bus")
    mx = matplotlib.ticker.MultipleLocator(2) # Define interval of minor ticks
    ax.xaxis.set_minor_locator(mx) # Set interval of minor ticks
    PyPlot.xticks(1:2:T)
    PyPlot.setp(ax.get_yticklabels(), size = labelsize)
    PyPlot.setp(ax.get_xticklabels(), size = labelsize)
    PyPlot.axis("tight")
    PyPlot.xlim(minimum(tspan), maximum(tspan))
    # grid(b=true,which="minor")
    PyPlot.tight_layout(pad = layout_pad)
end