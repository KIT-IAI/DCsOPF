export find_load_buses,
        get_unc_stor_sets,
        load_cost_comptime,
        formula_nr_decision_variables,
        balancing

function find_load_buses(ref_::Dict{Symbol, Any}; max_load::Bool=false)

    bus_loads = ref_[:bus_loads]
    loads = ref_[:load]
    
    loads_ = Dict()
    for (load_bus,load_idx) in bus_loads
        if ! isempty(load_idx)
            idx = load_idx[1]
            load = loads[idx]["pd"]
            loads_[load_bus] = load
        end
    end
    if max_load
        loads_ = sort(collect(loads_), by=x->-x[2]) # max load first
    end
    return loads_
end

function get_unc_stor_sets(case_nr::String, nr_unc::Int64, nr_stor::Int64; max_load::Bool=true)
    data = PowerModels.parse_file("casefiles/case" * string(case_nr) * ".m")
    ref_ = PowerModels.build_ref(data)[:it][:pm][:nw][0]
    loads_ = find_load_buses(ref_, max_load=max_load)
    loads_buses_ = [bus for (bus, load) in loads_]
    uncertainties = loads_buses_[1:nr_unc]
    storages = collect(setdiff(keys(ref_[:bus]), uncertainties))[1:nr_stor] # random location!
    uncertainties_subsets = [uncertainties[1:n] for n in 1:length(uncertainties)]
    storages_subsets = [storages[1:n] for n in 1:length(storages)]

    return uncertainties_subsets, storages_subsets
end

function formula_nr_decision_variables(n_unc::Int64, n_stor::Int64, n_gen::Int64,T::Int64)
    return (n_gen + n_stor) * (T + n_unc * T * (T+1) / 2)
end

function load_cost_comptime(dir)

    costs = Dict()
    comp_times = Dict()
    
    try
        costs = load(dir*"/costs.jld2")["costs"]
        comp_times = load(dir*"/comp_time.jld2")["comp_time"]
    catch
    end
    return costs, comp_times
end

function balancing(mysys::PowerSystem, optimizer; glob::Bool=false, storage::Bool=false, variance::Bool=false, text::String="Start balancing")
    # optimizer_with_attributes(Mosek.Optimizer,  "Threads" => 2)
    # options: Mosek.MSK_IINF_INTPNT_NUM_THREADS (interior-point algorithm), Mosek.MSK_IPAR_NUM_THREADS (optimizer)
    opf = buildDCsOPF_core(mysys, optimizer; storage = storage, glob = glob)
    if storage && variance
        addSpecificConstraints_Storage_WithVariance!(mysys, opf; glob = glob)
    elseif storage
        addSpecificConstraints_Storage!(mysys, opf; glob = glob)
    else
        addSpecificConstraints_NoStorage!(mysys, opf; glob = glob)
    end
    addAllChanceConstraints!(mysys, opf; glob = glob)
    optimize!(opf.m)
    return opf
end

function plot_complexity(d::Dict{Any,Any}, mysys::PowerSystem, case_nr::String, ylab::String, directory::String; type::String="cases", szenario::Int64=3, risk_level::Symbol=:risk_05, loc::Bool=true, glob::Bool=false, plot_dec_var::Bool=false, save::Bool=false)
    @assert loc || glob
    @assert type in ["cases", "risk levels"]
    gr()
    p = Plots.plot()
    x = 1:length(keys(d))
    ts = 2
    fs = 6
    lst_nr_unc_stor = sort([(parse(Int,split(n,"+")[1]),parse(Int,split(n,"+")[2])) for n in sort(collect(keys(d)))])
    # risk_levels = [:risk_025, :risk_05, :risk_10]
    risk_levels = collect(keys(d[collect(keys(d))[1]]))
    if type == "cases" # szenarios (no storage, storage, storage + var. constr.)
        y_loc = loc ? [[d[string(string(n_unc),"+",string(n_stor))][risk_level][:local][s] for (n_unc, n_stor) in lst_nr_unc_stor] for s in [1,2,3]] : []
        y_glob = glob ? [[d[string(string(n_unc),"+",string(n_stor))][risk_level][:global][s] for (n_unc, n_stor) in lst_nr_unc_stor] for s in [1,2,3]] : []
        labels_loc = loc ? ["No storage (loc)" "Storage (loc)" "Storage + var. constr. (loc)"] : []
        labels_glob = glob ? ["No storage (glob)" "Storage (glob)" "Storage + var. constr. (glob)"] : []
    elseif type == "risk levels"
        y_loc = loc ? [[d[string(string(n_unc),"+",string(n_stor))][r][:local][szenario] for (n_unc, n_stor) in lst_nr_unc_stor] for r in risk_levels] : []
        y_glob = glob ? [[d[string(string(n_unc),"+",string(n_stor))][r][:global][szenario] for (n_unc, n_stor) in lst_nr_unc_stor] for r in risk_levels] : []
        labels_loc = loc ? ["2,5% (loc)" "5% (loc)" "10% (loc)"] : []
        labels_glob = glob ? ["2,5% (glob)" "5% (glob)" "10% (glob)"] : []
    end
    y = [y_loc; y_glob]
    title = replace(string(string(case)," ",ylab), "_" => "-")
    #labels = loc ? (glob ? hcat(labels_loc, labels_glob) : labels_loc) : labels_glob
    xticks = (x, lst_nr_unc_stor)

    if plot_dec_var
        #  title = title,
        plot!(p, x, y_loc, label = labels_loc, ls=:solid, thickness_scaling=ts, right_margin=16mm, bottom_margin=3mm, xlabel="number of uncertainties and storages", ylabel=ylab, show=false, xticks=xticks, xtickfontsize=fs-1, ytickfontsize=fs, xlabelfontsize=fs, ylabelfontsize=fs, legendfontsize=fs)
        plot!(p, x, y_glob, label = labels_glob, ls=:dash, thickness_scaling=ts, xtickfontsize=fs-1, ytickfontsize=fs, xlabelfontsize=fs, ylabelfontsize=fs, legendfontsize=fs)
        nr_gen, T= mysys.inds[:N][:gen], mysys.horizon
        # uncertainties = mysys.inds[:dist][:unc][:buses]
        lst_nr_dec_var = [formula_nr_decision_variables(n_unc, n_stor, nr_gen, T) for (n_unc, n_stor) in lst_nr_unc_stor]
        scatter!(Plots.twinx(), lst_nr_dec_var, lw=3, lc="red", label="number of decision variables", ylabel="number of decision variables", show=true, xticks=xticks, legend=:right, xtickfontsize=fs-1, ytickfontsize=fs, xlabelfontsize=fs, ylabelfontsize=fs, legendfontsize=fs)
    else
        # title = title, 
        plot!(p, x, y_loc, label = labels_loc, ls=:solid, thickness_scaling=ts, right_margin=16mm, bottom_margin=3mm, xlabel="number of uncertainties and storages", ylabel=ylab, show=false, xticks=xticks, xtickfontsize=fs-1, ytickfontsize=fs+2, xlabelfontsize=fs, ylabelfontsize=fs+2, legendfontsize=fs)
        plot!(p, x, y_glob, label = labels_glob, ls=:dash, thickness_scaling=ts, legend_position=:left, xtickfontsize=fs-1, ytickfontsize=fs+2, xlabelfontsize=fs, ylabelfontsize=fs+2, legendfontsize=fs)
    end
    if save
        png(directory*"/complexity_"*ylab*"_"*type)
    end
end