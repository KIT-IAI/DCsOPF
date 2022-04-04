export cost,
        momentsOverHorizon,
        plot_unit,
        plot_map,
        plot_unit_2in1,
        plot_bars_,
        plot_complexity_cases,
        plot_complexity,
        plot_scenario,
        save_scenario,
        # save_costs,
        my_savefig,
        export_opf,
        export_opf_uncertainties,
        export_opf_load,
        export_opf_

cost(opf::OPF) = opf.fval(opf.m)

function map2name(x::Symbol)
    if x == :s || x == :e
        return :stor
    elseif x == :u || x == :Δu
        return :gen
    elseif x == :pl
        return :line
    end
end

function momentsOverHorizon(i::Int64, opf::OPF, ps::PowerSystem;kind::Symbol = :undef,glob::Bool = false)
    T = ps.horizon; tspan = 1:T

    lb, ub = ps.con[:val][kind][:min][i,tspan], ps.con[:val][kind][:max][i,tspan]
    λlb, λub = ps.con[:λ][kind][:min][i,tspan], ps.con[:λ][kind][:max][i,tspan]
    kind in [:Δu,:e] ? (tspan = 2:T) : ()
    μ, σ = zeros(length(tspan)), zeros(length(tspan))
    for (j, t_) in enumerate(tspan)
        # display((j,t_))
        μ[j], VAR = moments(i, t_, opf, ps;eval = true,kind = kind,glob = glob)
        σ[j] = norm(VAR, 2)
    end

    if kind == :Δu
        μ, σ = [0.; μ], [0.; σ]
    elseif kind == :e
        μ, σ = [ ps.stor[i][:ic][:μ]; μ ], [ ps.stor[i][:ic][:σ²]; σ ]
    end
    tspan = 1:T

    return Dict(:tspan => tspan, :μ => μ, :σ => σ, :lb => lb, :ub => ub, :λlb => λlb, :λub => λub)
end

# plotting

function plot_unit(i::Int64,opf::OPF,ps::PowerSystem;fignum::Int64 = 1,kind::Symbol = :undef,plot_lb::Bool = true,plot_ub::Bool = false,
    col::String = "blue",colshade::String = "#eeeeee",glob::Bool = false,alpha::Float64 = 0.1,
    lstyle::String = "solid",lwidth::Float64 = 0.75, figsize::Tuple = (2.75, 1.75),
    fontdict::Dict = Dict("color" => "black", "size" => 9, "family" => "serif"),
    labelsize = 7,
    layout_pad = 0)
    @assert kind in [:u,:Δu,:s,:e,:pl] "invalid keyword value for 'kind'"
    @assert i in 1:ps.inds[:N][map2name(kind)] "Specified unit $i does not exist (at most $(ps.inds[:N][map2name(kind)]) units)"

    vals = momentsOverHorizon(i, opf, ps;kind = kind,glob = glob)
    # μ, σ, lb, ub, λlb, λub = vals[:μ], vals[:σ], vals[:lb], vals[:ub], vals[:λlb], vals[:λub]
    # tspan = vals[:tspan]
    f = figure(num = fignum, figsize = figsize)
    bus = kind == :pl ? i : ps.inds[map2name(kind)][:buses][i]
    ylab = kind == :Δu ? latexstring("\\Delta u_{$bus}(t)") : latexstring("$(string(kind))_{$bus}(t)")
    _plot(vals,f,kind,ylab;plot_lb = plot_lb,plot_ub = plot_ub,col = col,colshade = colshade,alpha = alpha,lstyle = lstyle,lwidth = lwidth,
        fontdict = fontdict, labelsize = labelsize, layout_pad = layout_pad)

end

function _plot(d::Dict,fig::PyPlot.Figure,kind::Symbol,lab;plot_lb::Bool = true,plot_ub::Bool = true,
            col = "blue",colshade = "#eeeeee",alpha = 0.1,
            lstyle = "solid", lwidth::Float64 = 1.,
            fontdict::Dict = Dict("color" => "black", "size" => 9, "family" => "serif"),
            labelsize = 7,
            layout_pad = 0.)
    @assert 0 < alpha <= 1
    # println("PRINT BEGIN")
    # println("mu, lambdalb, lambdaub, sigma")
    # println(d[:μ])
    # println(d[:λlb])
    # println(d[:λub])
    # println(d[:σ])
    # println("Fill between")
    # println(d[:μ] - d[:λlb] .* d[:σ])
    # println(d[:μ] + d[:λub] .* d[:σ])
    # println("PRINT END")
    fill_between(d[:tspan], d[:μ] - d[:λlb] .* d[:σ], d[:μ] + d[:λub] .* d[:σ], color = colshade, alpha = alpha, linewidth = lwidth, linestyle = lstyle)
    grid(true); 
    plot(d[:tspan], d[:μ], color = col, linewidth = lwidth)
    plot_ub ? (plot(d[:tspan], d[:ub], color = "black", linestyle = "dashed")) : ()
    plot_lb ? (plot(d[:tspan], d[:lb], color = "black", linestyle = "dashed")) : ()
    xlabel(latexstring("t"), fontdict = fontdict)
    ylabel(lab, fontdict = fontdict)
    xlim(minimum(d[:tspan]), maximum(d[:tspan]))
    axis("tight")
    tight_layout(pad = layout_pad)
    ax = gca()
    majorformatter = matplotlib.ticker.FormatStrFormatter("%.1f")
    ax.yaxis.set_major_formatter(majorformatter)
    setp(ax.get_yticklabels(), size = labelsize)
    setp(ax.get_xticklabels(), size = labelsize)
    mx = matplotlib.ticker.MultipleLocator(2) # Define interval of minor ticks
    ax.xaxis.set_minor_locator(mx) # Set interval of minor ticks
    xticks(1:2:maximum(d[:tspan]))
    xlim([minimum(d[:tspan]),maximum(d[:tspan])])

end

function OLD_plot_unit(i::Int64,opf::OPF,ps::PowerSystem;fignum::Int64 = 1,kind::Symbol = :undef,plot_lb::Bool = true,plot_ub::Bool = true,
            col::String = "blue",colshade::String = "#eeeeee",glob::Bool = false,alpha::Float64 = 0.1,
            lstyle::String = "solid",lwidth::Float64 = 0., figsize::Tuple = (2.75, 1.75),
            fontdict::Dict = Dict("color" => "black", "size" => 9, "family" => "serif"),
            labelsize = 7,
            layout_pad = 0)
    @assert kind in [:u,:Δu,:s,:e,:pl] "invalid keyword value for 'kind'"
    @assert i in 1:ps.inds[:N][map2name(kind)] "Specified unit $i does not exist (at most $(ps.inds[:N][map2name(kind)]) units)"

    vals = momentsOverHorizon(i, opf, ps;kind = kind,glob = glob)
    println("Vals, mu and sigma:")
    println(vals[:μ])
    println(vals[:σ])
    # μ, σ, lb, ub, λlb, λub = vals[:μ], vals[:σ], vals[:lb], vals[:ub], vals[:λlb], vals[:λub]
    # tspan = vals[:tspan]
    f = figure(num = fignum, figsize = figsize)
    bus = kind == :pl ? i : ps.inds[map2name(kind)][:buses][i]
    ylab = kind == :Δu ? latexstring("\\Delta u_{$bus}(t)") : latexstring("$(string(kind))_{$bus}(t)")
    _plot(vals,f,kind,ylab;plot_lb = plot_lb,plot_ub = plot_ub,col = col,colshade = colshade,alpha = alpha,lstyle = lstyle,lwidth = lwidth,
            fontdict = fontdict, labelsize = labelsize, layout_pad = layout_pad)

end

function OLD_plot(d::Dict,fig::PyPlot.Figure,kind::Symbol,lab;plot_lb::Bool = true,plot_ub::Bool = true,
                    col = "blue",colshade = "#eeeeee",alpha = 0.1,
                    lstyle = "solid", lwidth::Float64 = 1.,
                    fontdict::Dict = Dict("color" => "black", "size" => 9, "family" => "serif"),
                    labelsize = 7,
                    layout_pad = 0.)
    println("PLOT")
    @assert 0 < alpha <= 1
    println("d, mu, lambdalb, lambdaub, sigma")
    println(d[:μ])
    println(d[:λlb])
    println(d[:λub])
    println(d[:σ])
    println("Color")
    println(col)
    println("Fill between")
    println(d[:μ] - d[:λlb] .* d[:σ])
    println(d[:μ] + d[:λub] .* d[:σ])
    plot(d[:tspan], d[:μ], color = col, linewidth = lwidth)
    plot(d[:tspan], d[:μ] .-2, color='y', linewidth=lwidth)
    fill_between(d[:tspan], d[:μ] - d[:λlb] .* d[:σ], d[:μ] + d[:λub] .* d[:σ], color = colshade, alpha = alpha, linewidth = lwidth, linestyle = lstyle)
    grid(true);
    plot_ub ? (plot(d[:tspan], d[:ub], color = "black", linestyle = "dashed")) : ()
    plot_lb ? (plot(d[:tspan], d[:lb], color = "black", linestyle = "dashed")) : ()
    # mu = [i for i in d[:μ]]
    # plot(d[:tspan], mu, color = "black", linewidth = lwidth)
    xlabel(latexstring("t"), fontdict = fontdict)
    ylabel(lab, fontdict = fontdict)
    xlim(minimum(d[:tspan]), maximum(d[:tspan]))
    axis("tight")
    tight_layout(pad = layout_pad)
    ax = gca()
    majorformatter = matplotlib.ticker.FormatStrFormatter("%.1f")
    ax.yaxis.set_major_formatter(majorformatter)
    setp(ax.get_yticklabels(), size = labelsize)
    setp(ax.get_xticklabels(), size = labelsize)
    mx = matplotlib.ticker.MultipleLocator(2) # Define interval of minor ticks
    ax.xaxis.set_minor_locator(mx) # Set interval of minor ticks
    xticks(1:2:maximum(d[:tspan]))
    xlim([minimum(d[:tspan]),maximum(d[:tspan])])

end

function plot_map(opf::OPF, ps::PowerSystem;fignum::Int64 = 1,kind::Symbol = :undef,exclude::Vector{Int64} = [-1],α::Float64 = 0.,scalefun::Function = x->x)
    @assert 0 <= α <= 1
    @assert kind in [:u,:Δu,:s,:e,:pl] "invalid keyword value for 'kind'"
    # scale     -->     scaling of relative values
    # α         -->     offset for plots, i.e. anything that is lower than
    #                   α*scale will not be displayed
    scale = 100
    n = ps.inds[:N][map2name(kind)]
    T = ps.horizon
    t = repeat(collect(1:T)', n, 1)
    u = repeat(collect(1:n), 1, T)
    z = zeros(n, T)
    col = repeat(["blue"], n, T)
    for i = 1:n
        vals = momentsOverHorizon(i, opf, ps;kind = kind)
        pmax = ps.con[:val][kind][:max][i,:]
        z[i,:] = (1 .- (pmax - (vals[:μ] + vals[:λub] .* vals[:σ])) ./ pmax) * scale
    end
    z_ = copy(z)
    # find all indices that should not be displayed
    # notice that z is going to be HEIGHT of the bars
    s = findall(x->x <= α * scale, z)
    z = z .- α * scale
    z[s] .= 0
    if exclude != [-1]
        z[exclude,:] = 0
        col[exclude,:] = "red"
    end


    f = figure(num = fignum, figsize = (10, 10))
    grid(true)
    # scatter(vec(t),vec(u),s=vec(z).^2,c=vec(z),cmap="Oranges")
    # scatter(vec(t),vec(u),s=0.5*vec(z).^2.5,c=vec(z),cmap="seismic")
    scatter(vec(t), vec(u), s = scalefun.(vec(z)), c = vec(z), cmap = "seismic")
    xticks(1:T)
    yticks(1:n)
    xlabel(latexstring("t"))
    ylabel(latexstring("$(kind)_i(t)"))
    # colorbar(orientation="horizontal")

    f = figure(num = fignum + 1, figsize = (10, 10))
    grid(true)
    # PyPlot.bar3d(vec(t),vec(u),zeros(n*T),0.75*ones(n*T),0.75*ones(n*T),vec(z))
    PyPlot.bar3d(vec(t), vec(u), α * scale * ones(n * T), 0.75 * ones(n * T), 0.75 * ones(n * T), vec(z), color = vec(col), alpha = 0.45)
    # scatter(vec(t),vec(u),s=vec(z),c=vec(z),cmap="Oranges")
    # xticks(1:T)
    # yticks(1:n)
    # zticks(scale*(α:0.1:1))
    xlabel(latexstring("t"))
    ylabel(latexstring("$(kind)_i(t)"))
    zlabel("Probability")
    zlim([α * scale,1.00 * scale])
    ylim([1,n])
    # colorbar(orientation="horizontal")
    ax = gca()
    ax[:view_init](35, -145)
    return f
end

function plot_unit_2in1(i::Int64,opf::OPF,ps::PowerSystem;fignum::Int64 = 1,kind::Symbol = :undef,col = "blue",figsize::Tuple = (2.75, 1.75),
    fontdict::Dict = Dict("color" => "black", "size" => 9, "family" => "serif"),
    labelsize = 7,
    layout_pad = 0,
    fontsize = 9)
    println(kind)
    @assert kind in [:u, :s]
    @assert i in 1:ps.inds[:N][:gen] "Specified unit $i does not exist (at most $(ps.inds[:N][map2name(:gen)]) units)"
    T = ps.horizon;
    tspan = 1:T
    if kind == :u
        x = :u; y = :Δu
    elseif kind == :s
        x = :s; y = :e
    end

    valsx = momentsOverHorizon(i, opf, ps;kind = x)
    μx, σx, lbx, ubx = valsx[:μ], valsx[:σ], valsx[:lb], valsx[:ub]
    λlbx, λubx = valsx[:λlb], valsx[:λub]
    tspan = valsx[:tspan]

    valsy = momentsOverHorizon(i, opf, ps;kind = y)
    μy, σy, lby, uby = valsy[:μ], valsy[:σ], valsy[:lb], valsy[:ub]
    λlby, λuby = valsy[:λlb], valsx[:λub]
    tspan = valsy[:tspan]

    fig = figure(num = fignum, figsize = figsize)
    
    # colorx="blue"
    colorx = col
    font1 = Dict("color" => colorx, "family" => "serif", "size" => fontsize)
    fill_between(tspan, μx - λlbx .* σx, μx + λubx .* σx, color = "#dddddd", alpha = 0.5)
    grid(true); 
    plot(tspan, μx, color = colorx)
    plot(tspan, ubx, color = "black", linestyle = "dashed")
    ax = gca()
    # if plt[:fignum_exists](fignum)==false
    xlabel(latexstring("t"))
    ylabel(latexstring(string(x) * "_{$i}(t)"), fontdict = font1)
    setp(ax[:get_yticklabels](), color = colorx, size = labelsize)
    # end
    # plot(tspan,lb,color="red",linestyle="dashed")


    # colory = "green"
    colory = col
    new_position = [0.06;0.06;0.77;0.91] # Position Method 2
    ax[:set_position](new_position)
    ax2 = ax[:twinx]() # Create another axis on top of the current axis
    font2 = Dict("color" => colory, "family" => "serif", "size" => fontsize)
    ylabel(latexstring(string(y) * "_{$i}(t)"), fontdict = font2)
    fill_between(tspan, μy - λlby .* σy, μy + λuby .* σy, color = "#dddddd", alpha = 0.5)
    plot(tspan, μy, color = colory, linestyle = "dotted")
    plot(tspan, uby, color = "black", linestyle = "dashed")
    plot(tspan, lby, color = "black", linestyle = "dashed")
    ax2[:set_position](new_position)
    # if plt[:fignum_exists](fignum)==false
    setp(ax2[:get_yticklabels](), color = colory, size = labelsize)
    xticks(tspan)
    # end
    # xlabel(latexstring("t"))
    axis("tight")
    tight_layout(pad = layout_pad)
    xlim(minimum(tspan), maximum(tspan))

    fig[:canvas][:draw]() # Update the figure

    return ax
end

function plot_complexity(d::Dict{Any,Any}, mysys::PowerSystem, case_nr::String, ylab::String, directory::String; type::String="cases", szenario::Int64=3, risk_level::Symbol=:risk_05, loc::Bool=true, glob::Bool=false, plot_dec_var::Bool=false, save::Bool=false)
    @assert loc || glob
    @assert type in ["cases", "risk levels"]
    gr()
    # figure(figsize = (10,5))
    display("test")
    p = Plots.plot(figsize = (10,5))
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

function plot_complexity_cases(cases::Array{String}, risk_level::Symbol, loc_glob::Symbol, type::String, deterministic::Bool, yscale::Symbol, szenario::Int64, directory::String; title::Bool=true, save::Bool=true)

    @assert type in ["cost", "time"]

    gr()
    p = Plots.plot()
    lst_nr_unc_stor = []
    lst_x = []
    figure(figsize = (5,5))

    # plot data
    for case in cases
        d = type=="cost" ? load(directory*"/case"*case*"/costs.jld2")["costs"] : load(directory*"/case"*case*"/comp_time.jld2")["comp_time"]
        lst_nr_unc_stor = union(lst_nr_unc_stor, sort([(parse(Int,split(n,"+")[1]),parse(Int,split(n,"+")[2])) for n in sort(collect(keys(d)))]))
        d_ = sort(d, by=x->parse(Int,split(x,"+")[1]))

        x = 1:length(d_)

        if case=="300"
            y = [d_[k][collect(keys(d_[k]))[1]][loc_glob][szenario] for k in keys(d_)]
        else
            y = [d_[k][risk_level][loc_glob][szenario] for k in keys(d_)]
        end

        @assert length(x) == length(y) string("Length of x and y do not match for case ", case)

        plot!(x, y, labels=replace(case, "deterministic" => "det"))

        lst_x = union(lst_x, x)
    end

    # layout
    title_ = ""
    if title
        title_ = "Complexity of different network sizes"
    end
    labels = cases
    xlabel= "number of (uncertainties, storages)"
    ylabel = type=="time" ? "Total computation time" : string("Total costs ")#,L"\$","/hr")
    xticks = (lst_x, lst_nr_unc_stor)
    fs = 6
    lw = 2
    plot!(p, label = labels, title = title_, yscale=yscale, xticks=xticks, ls=:solid, thickness_scaling=lw, left_margin=-5mm, right_margin=0mm, bottom_margin=-4mm, xlabel=xlabel, ylabel=ylabel, show=true, xtickfontsize=fs-1, ytickfontsize=fs, xlabelfontsize=fs, ylabelfontsize=fs, legendfontsize=fs)

    det = ""
    if deterministic
        det = "_deterministic"
    end

    if save
        png(directory*"/complexity_cases_"*type*"_"*string(yscale)*det)
    end

end

function plot_scenario(opf::OPF, mysys::PowerSystem, c::Int64, plot_gens::UnitRange{Int64}, plot_stor::UnitRange{Int64}, plot_lines::UnitRange{Int64}; col::String="#ff0000", α::Float64=0.15, figsize=figsize, linewidth=linewidth, glob::Bool=true, plot_ub::Bool=false, stor::Bool=true, var::Bool=false)
    lstyle = var ? "dashed" : "solid"
    # generation
    [ plot_unit(i, opf, mysys;fignum = 100 + c + i,figsize = figsize,kind = :u,col = col,colshade = col,plot_lb = false,plot_ub = plot_ub,alpha = α,
            glob = glob, lstyle = lstyle,lwidth = linewidth) for i in plot_gens ]
    # change in generation
    [ plot_unit(i, opf, mysys;fignum = 150 + c + i,figsize = figsize,kind = :Δu,col = col,colshade = col,plot_lb = false,plot_ub = plot_ub,alpha = α,
            glob = glob, lstyle = lstyle,lwidth = linewidth) for i in plot_gens ]
    if stor
        # storage injections
        [ plot_unit(i, opf, mysys;fignum = 600 + c + i,figsize = figsize,kind = :s,col = col,colshade = col,plot_lb = false,plot_ub = plot_ub,alpha = α,
                glob = glob, lstyle = lstyle,lwidth = linewidth) for i in plot_stor ]
        # storage states
        [ plot_unit(i, opf, mysys;fignum = 650 + c + i,figsize = figsize,kind = :e,col = col,colshade = col,plot_lb = false,plot_ub = plot_ub, alpha = α,
                glob = glob, lstyle = lstyle,lwidth = linewidth) for i in plot_stor ]
    end
    # branches
    [ plot_unit(i, opf, mysys; fignum = 700 + c + i, kind=:pl,alpha = α,col=col,plot_lb=false,colshade=col,glob=false) for i in plot_lines ]

end

function plot_szenario_OLD(scenario::Int64, temp, mysys, plot_gens, plot_stor, col, α, global_opt, plot_ub, plot_curves_single, plot_curves_together)
    
    func, c, nrs = choose_scenario(scenario)

    # plot
    if plot_curves_single && plot_ub
        fignums = nrs .+ c
        func(temp, mysys, fignums, plot_gens, plot_stor, col, α, global_opt, plot_ub)
    end
    if plot_curves_single
        fignums = nrs .+ 4000 .+ c
        func(temp, mysys, fignums, plot_gens, plot_stor, col, α, global_opt, false)
    end
    if plot_curves_together && plot_ub
        fignums = nrs .+ 4000
        func(temp, mysys, fignums, plot_gens, plot_stor, col, α, global_opt, plot_ub)
    end
    if plot_curves_together
        fignums = nrs .+ 8000
        func(temp, mysys, fignums, plot_gens, plot_stor, col, α, global_opt, false)
    end
end

function plot_bars_(dict, directory::String, risk_levels, nr_unc_stor::String, type::String, y_title::String; offset_text::Float64=2000.0, local_opt::Bool = true, global_opt::Bool = false, fignum::Int64 = 999, save::Bool=true, figsize::Tuple = (6, 2), layout_pad = 0) 
    close(fignum)
    # setup
    width = 0.07
    tspan = 1:3
    offsets = [-0.3, 0, 0.3]
    # offset_text = 2000
    # cmax = maximum(d[:risk_025][:global])
    figure(num = fignum, figsize = figsize)
    ax = gca()

    d = dict[nr_unc_stor]

    i = 1
    min_y = minimum(d[collect(keys(d))[1]][:local])
    max_y = maximum(d[collect(keys(d))[1]][:local])

    risks = [r for r in [:risk_025, :risk_05, :risk_10] if r in keys(d)] # sorted risks

    for r in risks
        pos = 0
        d_loc = d[r][:local]
        d_glob = d[r][:global]
        if local_opt 
            pos = collect(tspan).+ offsets[i] .- width/2
            PyPlot.bar(pos, d_loc, width = width, color="#1f77b4", label = "Local")
            min_y = minimum(d_loc) < min_y ? minimum(d_loc) : min_y
            max_y = maximum(d_loc) > max_y ? maximum(d_loc) : max_y
        end
        if global_opt
            pos = collect(tspan).+ offsets[i] .+ width/2
            PyPlot.bar(pos, d_glob, width = width, color="#ff7f0e", label = "Global")
            min_y = minimum(d_glob) < min_y ? minimum(d_glob) : min_y
            max_y = maximum(d_glob) > max_y ? maximum(d_glob) : max_y
        end
        lab = string(risk_levels[r]*100) * L"\%"
        annotate(lab, (tspan[1] + offsets[i], maximum([d_loc[1], d_glob[1]])+offset_text))
        annotate(lab, (tspan[2] + offsets[i], maximum([d_loc[2], d_glob[2]])+offset_text))
        annotate(lab, (tspan[3] + offsets[i], maximum([d_loc[3], d_glob[3]])+offset_text))
        i = i+1
    end

    # layout
    ax.legend(["Local","Global"])
    PyPlot.grid("on", axis="y")
    ylabel(y_title)
    PyPlot.ylim(0.8 * min_y, 1.2 * max_y)
    tight_layout(pad = layout_pad)
    PyPlot.xticks(ticks=[1,2,3], labels=["c1 - no storage", "c2 - storage", "c3 - storage and constr."])

    save ? my_savefig(directory, fignum, fignum, "bars_"*type) : nothing
end

function plot_bars_percent(dict, directory::String, risk_levels, nr_unc_stor::String, type::String, y_title::String; offset_text::Float64=2000.0, local_opt::Bool = true, global_opt::Bool = false, fignum::Int64 = 999, save::Bool=true, figsize::Tuple = (6, 2), layout_pad = 0, yscale::String="linear", upper_space::Float64=1.2) 
    close(fignum)
    # setup
    width = 0.07
    tspan = 1:3
    offsets = [-0.3, 0, 0.3]
    figure(num = fignum, figsize = figsize)
    ax = gca()

    d = dict[nr_unc_stor]

    ref_val = d[:risk_05][:local][2] # 5% risk, local opt., storage (no variance constraints)

    i = 1
    min_y = 0 # minimum(d[collect(keys(d))[1]][:local])
    max_y = 0 # maximum(d[collect(keys(d))[1]][:local])

    risks = [r for r in [:risk_025, :risk_05, :risk_10] if r in keys(d)] # sorted risks

    for r in risks
        pos = 0
        y_loc = d[r][:local]
        y_glob = d[r][:global]
        if local_opt 
            pos = collect(tspan).+ offsets[i] .- width/2
            y_loc = y_loc / ref_val .- 1
            PyPlot.bar(pos, y_loc, width = width, color="#1f77b4", label = "Local")
            min_y = minimum(y_loc) < min_y ? minimum(y_loc) : min_y
            max_y = maximum(y_loc) > max_y ? maximum(y_loc) : max_y
        end
        if global_opt
            pos = collect(tspan).+ offsets[i] .+ width/2
            y_glob = y_glob / ref_val .- 1
            PyPlot.bar(pos, y_glob, width = width, color="#ff7f0e", label = "Global")
            min_y = minimum(y_glob) < min_y ? minimum(y_glob) : min_y
            max_y = maximum(y_glob) > max_y ? maximum(y_glob) : max_y
        end
        lab = string(risk_levels_[r]*100) * L"\%"
        annotate(lab, (tspan[1] + offsets[i], maximum([y_loc[1], y_glob[1]])+offset_text))
        annotate(lab, (tspan[2] + offsets[i], maximum([y_loc[2], y_glob[2]])+offset_text))
        annotate(lab, (tspan[3] + offsets[i], maximum([y_loc[3], y_glob[3]])+offset_text))
        i = i+1
    end

    display(min_y)
    display(max_y)
    # layout
    ax.legend(["Local","Global"])
    PyPlot.grid("on", axis="y")
    PyPlot.yscale(yscale)
    ylabel(y_title)
    PyPlot.ylim(0.8 * min_y, upper_space * max_y)
    tight_layout(pad = layout_pad)
    PyPlot.xticks(ticks=[1,2,3], labels=["c1 - no storage", "c2 - storage", "c3 - storage and constr."])

    save ? my_savefig(directory, fignum, fignum, "bars_"*type*"_percent_"*yscale) : nothing
end

function plot_costs_(d::Dict;fignum::Int64 = 999,figsize::Tuple = (8, 5), layout_pad = 0) 
    width = 0.35
    tspan = 1:3
    # cmax = maximum(d[:global])
    figure(num = fignum, figsize = figsize)        
    bar(collect(tspan).- width/2, d[:local], width = width, label = "Local")
    bar(collect(tspan).+ width/2, d[:global], width = width, label = "Global")
    legend()
    ylabel("Total cost / hour")
    ylim(0.95 * maximum(d[:global]), 590000)
    tight_layout(pad = layout_pad)
    xticks([1,2,3])
end

# realizations

function makeL(X, i::Int64, t::Int64, n::Int64)
    Y = zeros(t, t)
    [Y[k,j] = X[i,k,j,n] for k = 1:t for j = 1:k]
    return Y
end

function getRealization_su_i(i::Int64, x, X, Ξ, xsym;glob::Bool = false)
    @assert xsym in [:u :s]
    @assert i >= 1
    nu, T = size(x)
    ndunc, T_, N = size(Ξ)
    @assert T == T_
    gen = Dict()
    if glob == false
            # local balancing
        u  = repeat(x[i,:], 1, N) + sum(makeL(X, i, T, n) * Ξ[n,:,:] for n = 1:ndunc)
        Δu = [0. ; diff(u)]
    else
            # global balancing
        error("not yet implemented")
    end
    return Dict(xsym => u, Symbol(:Δ, xsym) => Δu)
end

function getRealization_su(x, X, Ξ, xsym;glob::Bool = false)
    gen = Dict()
    nu, T = size(x)
    [ gen[i] = getRealization_su_i(i, x, X, Ξ, xsym;glob = glob) for i = 1:nu]
    return gen
end

function getRealization_e(S::Dict, ps::PowerSystem)
    ns = ps.inds[:N][:stor]
    @assert length(S) == ns
    T, Nscen = size(S[1][:s])
    @assert T == ps.horizon
    E = Dict()
    for i = 1:ns
        e, h = zeros(T, Nscen), ps.stor[i][:h]
        @assert abs(ps.stor[i][:ic][:σ²]) < 1e-10 "uncertain initial condition not supported"
        e[1,:] = ps.stor[i][:ic][:μ]
        [ e[t,:] = e[t - 1,:] - h * S[i][:s][t - 1,:] for t = 2:T ]
        Δe = [ zeros(1, Nscen); diff(e) ]
        E[i] = Dict(:e => e, :Δe => Δe)
    end
    return E
end

# getRealization_ppl(d::Dict,ps::PowerSystem) = getRealization_ppl(d[:dist][:cert],d[:dist][:unc],d[:gen],d[:stor],ps)

function getRealization_ppl(cert::Dict, unc::Dict, gen::Dict, stor::Dict, ps::PowerSystem)
    T, nbus, nline = ps.horizon, ps.inds[:N][:bus], ps.inds[:N][:line]
    T_, Nscen = size(cert[1])
    @assert T == T_
    # pow = Dict()
    # for i=1:nbus
    #     p = zeros(T,Nscen)
    #     # generated power
    #     [ p += gen[j][:u] for j=1:ps.inds[:N][:gen] ]
    #     # certain disturbance
    #
    #     # uncertain disturbance
    #
    #     # storage
    # end

    pow = zeros(nbus, T, Nscen)
    # add generated power for every bus
    [ pow[ps.inds[:gen][:buses][i],:,:] += gen[i][:u] for i = 1:ps.inds[:N][:gen]]
    # add certain disturbance for every bus
    [ pow[ps.inds[:dist][:cert][:buses][i],:,:] -= cert[i] for i = 1:ps.inds[:N][:dist][:cert]]
    # add uncertain disturbance for every bus
    [ pow[ps.inds[:dist][:unc][:buses][i],:,:] -= unc[i] for i = 1:ps.inds[:N][:dist][:unc]]
    # add storage power for every bus
    [ pow[ps.inds[:stor][:buses][i],:,:] += stor[i][:s] for i = 1:ps.inds[:N][:stor]]

    pl = zeros(nline, T, Nscen)
    [ pl[:,t,n] = ps.ptdf * pow[:,t,n] for t = 1:T for n = 1:Nscen ]
    d = Dict()
    [ d[i] = pl[i,:,:] for i = 1:nline ]

    return pow, d

    # [ d[i] = pow[i,:,:] for i=1:nbus ]
    # return d
end

function getRealization(opf::OPF, ps::PowerSystem, N::Int64 = 1;glob::Bool = false)
    T = ps.horizon
    nu, ndcert, ndunc, nline = ps.inds[:N][:gen], ps.inds[:N][:dist][:cert], ps.inds[:N][:dist][:unc], ps.inds[:N][:line]
    Ξ = randn(ndunc, T, N)

    unc = Dict()
    [ unc[i] = -(repeat(ps.dist[:unc][i][:μ], 1, N) + ps.dist[:unc][i][:Σ] * Ξ[i,:,:]) for i = 1:ndunc ]
    cert = Dict()
    [ cert[i] = -repeat(ps.dist[:cert][i], 1, N) for i = 1:ndcert ]
    gen = getRealization_su(getvalue(opf.u), getvalue(opf.U), Ξ, :u;glob = glob)
    stor = getRealization_su(getvalue(opf.s), getvalue(opf.S), Ξ, :s;glob = glob)
    # display(stor)
    ener = getRealization_e(stor, ps)

    pow, pl = getRealization_ppl(cert, unc, gen, stor, ps)

    Dict(:dist => Dict(:unc => unc, :cert => cert),
         :gen => gen,
         :stor => Dict(:power => stor, :energy => ener),
         :pline => pl)
    # (unc,cert,gen,stor,ener)
end

# save and export

function my_savefig_OLD(case_nr::String, n::Int64, horizon::Int64, fignum::Int64; prefix::String = "myfile",fileformat::String = "jpg", subfolder::String = "")
    figure(fignum)
    gca()
    dir = "graphics/results/case"*case_nr*"/" * subfolder * "/"
    isdir(dir) ? nothing : mkdir(dir) 
    fname = dir * prefix * "_" * string(fignum) * "." * fileformat
    savefig(fname, format = fileformat, pad_inches = 0)
end

function my_savefig(dir::String, n::Int64, fignum::Int64, prefix::String; fileformat::String="jpg")
    PyPlot.figure(fignum)
    PyPlot.gca()
    directory = dir * "/" * prefix * "_" * string(fignum) * "." * fileformat
    PyPlot.savefig(directory, format = fileformat, pad_inches = 0)
end

# function save_comp_time(directory, d)
#     # e.g. loaded: d["costs"]["57"][:risk_025][:local]
#     save(directory * "/comp_time.jld2", "comp_time", d)
# end

function save_scenario(dir::String, opf::OPF, mysys::PowerSystem, c::Int64, plot_gens::UnitRange{Int64}, plot_stor::UnitRange{Int64}, plot_lines::UnitRange{Int64})

    save_u  = plot_gens .+ (100 + c)
    [ my_savefig(dir, mysys.inds[:gen][:buses][i - 100-c], i, "gen_u") for i in save_u]

    save_Δu = plot_gens .+(150 + c)
    [ my_savefig(dir, mysys.inds[:gen][:buses][i - 150-c], i, "gen_delta_u") for i in save_Δu]
    
    if opf.stor
        save_s  = plot_stor .+ (600 + c)
        [ my_savefig(dir, mysys.inds[:stor][:buses][i - 600-c], i, "storage_s") for i in save_s]
        save_e  = plot_stor .+ (650 + c)
        [ my_savefig(dir, mysys.inds[:stor][:buses][i - 650-c], i, "storage_e") for i in save_e]
    end

    save_pl = plot_lines .+ (700 + c)
    [ my_savefig(dir, i, i, "/line") for i in save_pl]
end

# function save_costs(directory, d)
#     # e.g. loaded: d["costs"]["57"][:risk_025][:local]
#     save(directory * "/costs.jld2", "costs", d)
# end

function export_opf_(opf::OPF, ps::PowerSystem, name::String, type::Symbol = :undef, kind::Symbol = :undef, glob::Bool=false)
    @assert kind in [:u,:e,:pl] "invalid keyword value for 'kind'"
    d = Dict()
    mean = Dict()
    var = Dict()
    ub = Dict()
    if kind == :u || kind == :e # for gen/stor take bus nunber
        buses = ps.inds[type][:buses]
    elseif kind == :pl # for lines simply enumerate
        buses = collect(1:ps.inds[:N][type])
    end
    for i in 1:ps.inds[:N][type]
        bus = buses[i]
        vals = momentsOverHorizon(i, opf, ps; kind = kind, glob = glob)
        mean[string("bus_",bus)] = vals[:μ]
        var[string("bus_",bus)] = vals[:σ]
        ub[string("bus_",bus)] = vals[:ub]
    end
    d["mean"] = mean
    d["var"] = var
    d["ub"] = ub
    return d
end

function export_opf_load(ps::PowerSystem)
    d = Dict()
    mean = Dict()
    var = Dict()
    for (i,bus) in enumerate(ps.inds[:dist][:unc][:buses])
        mean[string("bus_", bus)] = ps.dist[:unc][i][:μ]
        var[string("bus_", bus)] = diag(ps.dist[:unc][i][:Σ_full])
    end
    for (i,bus) in enumerate(ps.inds[:dist][:cert][:buses])
        mean[string("bus_", bus)] = ps.dist[:cert][i]
    end
    d["mean"] = mean
    d["var"] = var
    return d
end

function export_opf_uncertainties(ps::PowerSystem)
    return ps.inds[:dist][:unc][:buses]
end

function export_opf(opf::OPF, ps::PowerSystem, directory::String, glob::Bool = false, storage::Bool = false, save_name::String = "DCsOPF_NoName")
    opf_dict = Dict()
    names = ["lines", "storages", "generators"]
    types = [:line, :stor, :gen]
    kinds = [:pl, :e, :u]
    for i in 1:3
        if ~(kinds[i] == :e && ~storage)
            opf_dict[names[i]] = export_opf_(opf, ps, names[i], types[i], kinds[i], glob)
        end
    end
    opf_dict["load"] = export_opf_load(ps)
    opf_dict["uncertainties"] = export_opf_uncertainties(ps)

    return opf_dict
    # matwrite(string(directory, "Results_", save_name, ".mat"), opf_dict)
end