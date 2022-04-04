export OPF, PowerSystem, CaseParams

struct OPF
    m::JuMP.Model
    u::Matrix{JuMP.VariableRef}
    U::JuMP.Containers.SparseAxisArray
    s::Matrix{JuMP.VariableRef}
    S::JuMP.Containers.SparseAxisArray
    con::Dict
    fval::Function
    stor::Bool
end

struct PowerSystem
    horizon::Int64
    inds::Dict
    dist::Dict
    stor::Dict
    ptdf::Matrix{Float64}
    con::Dict
    cost::Dict
end

struct CaseParams
    case_nr::String
    uncertainties::Vector{Vector{Int64}}
    storages::Vector{Vector{Int64}}
    local_opt::Bool
    global_opt::Bool
    deterministic::Bool
    artificial_unc::Bool
    load_factor::Float64
    dist_factor::Float64
    unc_factor::Float64
    stor_ub::Int64
    γ::Float64
end