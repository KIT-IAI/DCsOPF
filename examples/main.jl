using DCsOPF, PyPlot, JuMP, MosekTools, PowerModels, MAT, LinearAlgebra, PowerModelsParsing, StatsFuns, JLD2, FileIO

cd("examples")
include("../src/helperfunctions.jl")
include("casefiles/buildsystem_all.jl")

##############################################################################

# Setup
PyPlot.matplotlib.rc("text", usetex = true) # allow tex rendering
PyPlot.matplotlib.rc("font", family = "serif")
fontsize = 12 # 9
fontdict = Dict("color" => "black", "size" => fontsize, "family" => "serif")
labelsize = 7
figsize = (2.75, 1.75) # (width, height) in inch
layout_pad = 0.
linewidth = 0.75

# plot parameters
plot_curves = true
plot_scenarios_in_one = true
plot_scenarios_individually = true
plot_ub = true # plot upper bound
plot_ub_and_noub = true
load_old_costs_times = true
export_results = true
plot_curves_global = false

# covariance matrix
cov_mat = "CovMatrix"

# optimization
loc_opt = true
glob_opt = true

# risk levels 
risk_levels_dict = Dict(:risk_025 => 0.025, :risk_05 => 0.05, :risk_10 => 0.1)
risk_levels = [:risk_05]
# risk_levels = [:risk_025, :risk_05, :risk_10]

##############################################################################

# CASES
# Struct "CaseParams" inputs:
# case number, uncertainties, storages, deterministic (yes: covariance matrix = zeros), artificial uncertainties, ...
# ... load factor, disturbance factor, uncertainty factor, storage upper bound, maximum change of ±γ % for the generators (sine amplitude)
# get_unc/stor_sets creats several sets of uncertainties/storages of the form [[..],..,[..]]

params = Dict()

# case 5
unc_set, stor_set = get_unc_stor_sets("5", 3, 2)
stor_set = append!(stor_set, [last(stor_set)])
params["5_artificial"] = CaseParams("5", unc_set, stor_set, false, true, 1., 0.5, 8., 6, 0.8)
params["5_volatile"] = CaseParams("5", unc_set, stor_set, false, false, 1., 1., 10., 10, 0.8)

# case 39
unc_set, stor_set = get_unc_stor_sets("39", 8, 8)
# unc_set, stor_set = [], [[19,23,22,29,2,25]]
# params["39_artificial"]   = CaseParams("39", unc_set, stor_set, false, true, 1., 2., 20., 6, 0.2)
params["39_volatile"]   = CaseParams("39", unc_set, stor_set, false, false, 1., 1., 10., 100, 0.2) # false, false, 1./1.5, 1.5, 10./1.2-1.7, 5, 0.4

# case 57
unc_set, stor_set = get_unc_stor_sets("57", 8, 8)
# unc_set, stor_set = [[12,8,9,6,1,16,17,3,47,18,15]], [[50,53,49,51,12]]
# unc_set, stor_set = [[12,8,9,6,1,16,17]], [[50,53,49,51,12]]
# params["57_artificial"]  = CaseParams("57", unc_set, stor_set, false, true, 1., 1., 1., 4, 0.2)
params["57_volatile"]  = CaseParams("57", unc_set, stor_set, false, false, 1., 1., 20., 100, 0.8)
params["57_volatile"]  = CaseParams("57", unc_set, stor_set, false, false, 1., 0.5, 20., 100, 0.8) # for 4
params["57_volatile"]  = CaseParams("57", unc_set, stor_set, false, false, 2., 0.5, 20., 100, 0.8) # for 4

# case 118
unc_set, stor_set = get_unc_stor_sets("118", 8, 8)
unc_set, stor_set =  [[59,116,90,80,54,42,15,49,56,60]], [[11,62,27,78,112]]
params["118_volatile"] = CaseParams("118", unc_set, stor_set,false, false, 1., 1., 1., 4, 0.2)

# case 300
unc_set, stor_set = get_unc_stor_sets("300", 8, 8)
# unc_set, stor_set = 
params["300_volatile"] = CaseParams("300", unc_set, stor_set, false, false, 1., 40., 40., 100, 0.2)

##############################################################################
#### OPTIMIZATION
##############################################################################

dir_results = "results"
isdir(dir_results) ? nothing : mkdir(dir_results)

##############################################################################
## Cases
##############################################################################

for case in ["300_volatile"]
# for case in ["5_volatile"]
    # case = "39_volatile"

    close("all") # close all figures

    # parameters
    case_nr, uncertainties, storages, local_opt, global_opt, det, artificial_unc, load_factor, dist_factor, unc_factor, stor_ub, gen_change = [getfield(params[case], key) for key ∈ fieldnames(typeof(params[case]))]

    dir_case = dir_results*"/case"*case
    isdir(dir_case) ? nothing : mkdir(dir_case)

    # load (and then overwrite) or initialize cost and computation time
    costs, comp_times = load_cost_comptime(dir_case)

    # uncertainties, storages = [first(uncertainties)], [first(storages)]
    # uncertainties, storages = [last(uncertainties)], [last(storages)]
    # uncertainties, storages = [[4,2]], [[5,1]]
    if case_nr == "300"
        uncertainties, storages = [uncertainties[5], uncertainties[6], uncertainties[7], uncertainties[8]], [storages[5], storages[6], storages[7], storages[8]]
    end
    # uncertainties, storages = [uncertainties[5]], [storages[5]]

##############################################################################
## Uncertainties
##############################################################################

    for (unc, stor) in zip(uncertainties, storages)

        nr_unc_stor = string(length(unc))*"+"*string(length(stor)) # e.g. "3+2" with 3 uncertainties and 2 storages

        dir_unc_stor = dir_case*"/"*nr_unc_stor
        isdir(dir_unc_stor) ? nothing : mkdir(dir_unc_stor)
        dir_curves_loads = dir_case*'/'*nr_unc_stor*"/Curves_loads"
        isdir(dir_curves_loads) ? nothing : mkdir(dir_curves_loads)

        costs_ = Dict()
        comp_times_ = Dict()
        
##############################################################################
## Risk level
##############################################################################
        for risk in risk_levels
            # risk = :risk_05

            dir_curves = dir_unc_stor*"/Curves"
            isdir(dir_curves) ? nothing : mkdir(dir_curves)
            dir_exportedResults = dir_unc_stor*"/ExportedResults/"
            isdir(dir_exportedResults) ? nothing : mkdir(dir_exportedResults)
    
##############################################################################
## Build system
##############################################################################

            # build buildSystem
            global mysys, case_data = buildSystem(case_nr, unc, stor, params[case]; ε=risk_levels_dict[risk])
            display(checkGenerationCapacity(mysys))

            # save (un)certain load curves
            [my_savefig(dir_curves_loads, i, i, "/uncertain_disturbance") for i in mysys.inds[:dist][:unc][:buses]]
            [my_savefig(dir_curves_loads, i, i, "/certain_disturbance") for i in mysys.inds[:dist][:cert][:buses]]

##############################################################################
## Optimization
##############################################################################

            # balancing
            # @elapsed stores the computation time
            display("#########################################")
            display("################ START ##################")
            display("#########################################")
            ## NO STORAGE
            display(string("NoStorage, Case: ", case_nr, ", nr_unc_stor: ", nr_unc_stor, ", risk level: ", risk))
            t_loc1 = @elapsed DCsOPF_NoStorage_local = balancing(mysys, Mosek.Optimizer, glob=false, storage=false, text="No Storages + local balancing")
            t_glob1 = @elapsed DCsOPF_NoStorage_global = balancing(mysys, Mosek.Optimizer, glob=true, storage=false, text="No Storages + global balancing")
            
            ## STORAGE, NO VARIANCE
            display(string("Storage, Case: ", case_nr, ", nr_unc_stor: ", nr_unc_stor, ", risk level: ", risk))
            t_loc2 = @elapsed DCsOPF_Storage_local = balancing(mysys, Mosek.Optimizer, glob=false, storage=true, text="Storages + local balancing + no variance constraint")
            t_glob2 = @elapsed DCsOPF_Storage_global = balancing(mysys, Mosek.Optimizer, glob=true, storage=true, text="Storages + global balancing + no variance constraint")
            
            ## STORAGE, VARIANCE
            display(string("Storage and Var, Case: ", case_nr, ", nr_unc_stor: ", nr_unc_stor, ", risk level: ", risk))
            t_loc3 = @elapsed DCsOPF_Storage_local_withVariance = balancing(mysys, Mosek.Optimizer, glob=false, storage=true, variance=true, text="Storages + local balancing + variance constraint")
            t_glob3 = @elapsed DCsOPF_Storage_global_withVariance = balancing(mysys, Mosek.Optimizer, glob=true, storage=true, variance=true, text="Storages + global balancing + variance constraint")

            # export results
            try export_results
                    opf_dict = export_opf(DCsOPF_NoStorage_local, mysys, dir_exportedResults, false, false, "DCsOPF_NoStorage_local")
                    matwrite(string(dir_exportedResults, "Results_", "DCsOPF_NoStorage_local", ".mat"), opf_dict)
                    opf_dict = export_opf(DCsOPF_NoStorage_global, mysys, dir_exportedResults, true, false, "DCsOPF_NoStorage_global")
                    matwrite(string(dir_exportedResults, "Results_", "DCsOPF_NoStorage_global", ".mat"), opf_dict)
                    opf_dict = export_opf(DCsOPF_Storage_local, mysys, dir_exportedResults, false, true, "DCsOPF_Storage_local")
                    matwrite(string(dir_exportedResults, "Results_", "DCsOPF_Storage_local", ".mat"), opf_dict)
                    opf_dict = export_opf(DCsOPF_Storage_global, mysys, dir_exportedResults, true, true, "DCsOPF_Storage_global")
                    matwrite(string(dir_exportedResults, "Results_", "DCsOPF_Storage_global", ".mat"), opf_dict)
                    opf_dict = export_opf(DCsOPF_Storage_local_withVariance, mysys, dir_exportedResults, false, true, "DCsOPF_Storage_local_withVariance")
                    matwrite(string(dir_exportedResults, "Results_", "DCsOPF_Storage_local_withVariance", ".mat"), opf_dict)
                    opf_dict = export_opf(DCsOPF_Storage_global_withVariance, mysys, dir_exportedResults, true, true, "DCsOPF_Storage_global_withVariance")
                    matwrite(string(dir_exportedResults, "Results_", "DCsOPF_Storage_global_withVariance", ".mat"), opf_dict)
            catch
                display("Export did not work.")
            end

            # computation times 
            comp_times__ = Dict()
            comp_times__[:local] = [t_loc1, t_loc2, t_loc3]
            comp_times__[:global] = [t_glob1, t_glob2, t_glob3]
            comp_times_[risk] = comp_times__

            # costs
            costs__ = Dict()
            costs__[:local] = [ cost(DCsOPF_NoStorage_local), cost(DCsOPF_Storage_local), cost(DCsOPF_Storage_local_withVariance) ]
            costs__[:global] = [ cost(DCsOPF_NoStorage_global), cost(DCsOPF_Storage_global), cost(DCsOPF_Storage_global_withVariance) ]
            costs_[risk] = costs__

##############################################################################
## Plots: curves over horizon
##############################################################################

            ### PLOTTING
            if plot_curves && risk == :risk_05 # only plot for one risk level

                try
                    close("all") # close all figures

                    plotted_gens = 1:mysys.inds[:N][:gen]
                    plotted_stor = 1:mysys.inds[:N][:stor]
                    plotted_line = 1:mysys.inds[:N][:line]

                    # plot each scenario in own plot
                    if plot_scenarios_individually
                        dir = dir_curves * "/Single"
                        isdir(dir) ? nothing : mkdir(dir)

                        # plot generation, change of generation, storage injection, storge state & line flows
                        plot_scenario(DCsOPF_NoStorage_local, mysys, 1000, plotted_gens, plotted_stor, plotted_line; col="#ff0000", α=0.15, figsize=figsize, linewidth=linewidth, glob=false, plot_ub=plot_ub, stor=false, var=false)
                        plot_scenario(DCsOPF_Storage_local, mysys, 2000, plotted_gens, plotted_stor, plotted_line, col="#0000ff", α=0.15, figsize=figsize, linewidth=linewidth, glob=false, plot_ub=plot_ub, stor=true, var=false)
                        plot_scenario(DCsOPF_Storage_local_withVariance, mysys, 3000, plotted_gens, plotted_stor, plotted_line, col="#228B22", α=0.15, figsize=figsize, linewidth=linewidth, glob=false, plot_ub=plot_ub, stor=true, var=true)

                        # save
                        save_scenario(dir, DCsOPF_NoStorage_local, mysys, 1000, plotted_gens, plotted_stor, plotted_line)
                        save_scenario(dir, DCsOPF_Storage_local, mysys, 2000, plotted_gens, plotted_stor, plotted_line)
                        save_scenario(dir, DCsOPF_Storage_local_withVariance, mysys, 3000, plotted_gens, plotted_stor, plotted_line)
                    end

                    # plot all scenarios in one plot
                    if plot_scenarios_in_one
                        dir = dir_curves * "/Together"
                        isdir(dir) ? nothing : mkdir(dir)

                        # plot generation, change of generation, storage injection, storge state & line flows
                        plot_scenario(DCsOPF_NoStorage_local, mysys, 4000, plotted_gens, plotted_stor, plotted_line, col="#ff0000", α=0.15, figsize=figsize, linewidth=linewidth, glob=false, plot_ub=plot_ub, stor=false, var=false)
                        plot_scenario(DCsOPF_Storage_local, mysys, 4000, plotted_gens, plotted_stor, plotted_line, col="#0000ff", α=0.15, figsize=figsize, linewidth=linewidth, glob=false, plot_ub=plot_ub, stor=true, var=false)
                        plot_scenario(DCsOPF_Storage_local_withVariance, mysys, 4000, plotted_gens, plotted_stor, plotted_line, col="#228B22", α=0.15, figsize=figsize, linewidth=linewidth, glob=false, plot_ub=plot_ub, stor=true, var=true)

                        # save
                        save_scenario(dir, DCsOPF_Storage_local, mysys, 4000, plotted_gens, plotted_stor, plotted_line)

                        if plot_ub_and_noub
                            # plot generation, change of generation, storage injection, storge state & line flows
                            plot_scenario(DCsOPF_NoStorage_local, mysys, 8000, plotted_gens, plotted_stor, plotted_line, col="#ff0000", α=0.15, figsize=figsize, linewidth=linewidth, glob=false, plot_ub=!plot_ub, stor=false, var=false)
                            plot_scenario(DCsOPF_Storage_local, mysys, 8000, plotted_gens, plotted_stor, plotted_line, col="#0000ff", α=0.15, figsize=figsize, linewidth=linewidth, glob=false, plot_ub=!plot_ub, stor=true, var=false)
                            plot_scenario(DCsOPF_Storage_local_withVariance, mysys, 8000, plotted_gens, plotted_stor, plotted_line, col="#228B22", α=0.15, figsize=figsize, linewidth=linewidth, glob=false, plot_ub=!plot_ub, stor=true, var=true)

                            # save
                            save_scenario(dir, DCsOPF_Storage_local, mysys, 8000, plotted_gens, plotted_stor, plotted_line)
                        end
                    end

                catch
                    display("Plots did not work.")
                end

            end # plots

        end # risk levels

        comp_times[nr_unc_stor] = comp_times_
        costs[nr_unc_stor] = costs_

        save(dir_case * "/costs.jld2", "costs", costs)
        save(dir_case * "/comp_time.jld2", "comp_time", comp_times)

##############################################################################
## Plots for each uncertainty/storage set
##############################################################################

        try
            plot_bars_(costs, dir_unc_stor, risk_levels_dict, nr_unc_stor, "costs", string("Total costs ",L"\$","/hr"), offset_text=10000.0, local_opt=true, global_opt=true, fignum=999)
            plot_bars_(comp_times, dir_unc_stor, risk_levels_dict, nr_unc_stor, "computation time", "Total computation time (s)", offset_text=0.01, local_opt=true, global_opt=true, fignum=888)
        catch
            println("Bar plots did not work.")
        end

    end # uncertainties, storages

##############################################################################
## Plot for one case
##############################################################################

    using Plots, GR, Measures

    try
        # costs
        plot_complexity(costs, mysys, case, "costs", dir_case, type="cases", risk_level=:risk_05, loc=true, glob=true, plot_dec_var=true, save=true)
        plot_complexity(costs, mysys, case, "costs", dir_case, type="risk levels", szenario=3, loc=true, glob=true, plot_dec_var=true, save=true)
        
        # computation times
        plot_complexity(comp_times, mysys, case, "computation time", dir_case, type="cases", risk_level=:risk_05, loc=true, glob=true, plot_dec_var=true, save=true)
        plot_complexity(comp_times, mysys, case, "computation time", dir_case, type="risk levels", szenario=3, loc=true, glob=true, plot_dec_var=true, save=true)
    catch
        println("Complexity plots did not work.")
    end

end # cases

##############################################################################
## Plot for all cases
##############################################################################

# all cases (time, costs)
# plot_complexity_cases(["5", "39", "57", "118", "300"], :risk_05, :local, "time", false, :identity, 3, dir_results, title=false, save=true)
# plot_complexity_cases(["5", "39", "57", "118", "300"], :risk_05, :local, "cost", false, :identity, 3, dir_results, title=false, save=true)
# plot_complexity_cases(["5", "39", "57", "118", "300"], :risk_05, :local, "time", false, :log, 3, dir_results, title=false, save=true)
# plot_complexity_cases(["5", "39", "57", "118", "300"], :risk_05, :local, "cost", false, :log, 3, dir_results, title=false, save=true)
# plot_complexity_cases(["5", "5_deterministic", "39", "39_deterministic", "57", "57_deterministic", "118", "300"], :risk_05, :local, "time", true, :identity, 3, dir_results, title=false, save=true)
# plot_complexity_cases(["5", "5_deterministic", "39", "39_deterministic", "57", "57_deterministic", "118", "300"], :risk_05, :local, "cost", true, :identity, 3, dir_results, title=false, save=true)
# plot_complexity_cases(["5", "5_deterministic", "39", "39_deterministic", "57", "57_deterministic", "118", "300"], :risk_05, :local, "time", true, :log, 3, dir_results, title=false, save=true)
# plot_complexity_cases(["5", "5_deterministic", "39", "39_deterministic", "57", "57_deterministic", "118", "300"], :risk_05, :local, "cost", true, :log, 3, dir_results, title=false, save=true)