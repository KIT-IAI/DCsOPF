function buildCostSOC_perbus(x, X, quad, lin, nd::Int, bus::Int;glob::Bool = false)
    nu, T = size(x)
    # build component for mean
    y = [ x[bus,k] for k in 1:T ]
    # build components for Sigma
    r = Int(T * (T + 1) / 2)
    Y = Vector{JuMP.VariableRef}(undef, nd * r)
    for i = 1:nd
        rows = ((i - 1) * r + 1):(i * r)
        Y[rows] = 
        if glob
            # global balancing
            [ X[bus,k,l] for k in 1:T for l in 1:k ]
        else
            # local balancing
            [ X[bus,k,l,i] for k in 1:T for l in 1:k ]
        end
    end
    h = [ 0.5 * lin[bus] * ones(T); spzeros(nd * r) ]
    e = ones(T + nd * r)
    H = quad[bus] * e
    Hsq, Hsqi = sqrt(quad[bus]) * e, 1 / sqrt(quad[bus]) * e

    return [y; Y], h, H, Hsq, Hsqi
end

function buildCostSOC(x, X, quad, lin, nd;glob::Bool = false)
    nu, T = size(x)
    r = Int(T * (T + 1) / 2)
    ni = T + nd * r # number of decision variables per bus
    n = nu * ni   # number of decision variables in total
    Y, h, H, Hsq, Hsqi = Vector{JuMP.VariableRef}(undef, n), spzeros(n), zeros(n), zeros(n), zeros(n)
    for i in 1:nu
        rows = (i - 1) * ni + 1:(i * ni)
        Y[rows], h[rows], H[rows], Hsq[rows], Hsqi[rows] = buildCostSOC_perbus(x, X, quad, lin, nd, i;glob = glob)
    end
    Y, h, H, Hsq, Hsqi
end