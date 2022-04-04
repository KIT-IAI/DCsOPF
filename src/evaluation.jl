function evalCC(opf::OPF, s::PowerSystem, i::Int64, t::Int64, x::Symbol, kind::Symbol)
    addCC!(opf, s, i, t, x, kind;eval = true)
end

function evalConstraints(opf::OPF, ps::PowerSystem)
    u, U = getvalue(opf.u), getvalue(opf.U)
    T = ps.horizon
    nu, nl = ps.inds[:N][:gen], ps.inds[:N][:line]
    x = [:u, :Δu, :pl]
    N = [nu,  nu,  nl]
    k1 = [1,   2,  1]
    bnd = [:lb, :ub]

    cUmax = [ evalCC(opf, ps, i, k, :u, :ub) for i = 1:nu, k = 1:T ]
    cUmin = [ evalCC(opf, ps, i, k, :u, :lb) for i = 1:nu, k = 1:T ]
    d = Dict()
    m = Dict()
    for (i, x_) in enumerate(x)
        d[x_], m[x_] = Dict(), Dict()
        for bnd_ in bnd
            d[x_][bnd_] = [ evalCC(opf, ps, i, k, x_, bnd_) for i = 1:N[i], k = k1[i]:T ]
            m[x_][bnd_] = Dict(:val => maximum(d[x_][bnd_]), :ind => indmax(d[x_][bnd_]))
        end
    end
    d, m
end
