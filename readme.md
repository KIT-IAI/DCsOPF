# Analytical Uncertainty Propagation for Multi-Period Stochastic Optimal Power Flow

*How can we combine stochastic DC optimal power flow with Gaussian processes?*

This repository provides both the finished [paper](docs/main.pdf) and a Julia package to tackle this question.



## Developing with Julia

Getting started with Julia: [getting started](https://docs.julialang.org/en/v1/manual/getting-started/), [workflow tips](https://docs.julialang.org/en/v1/manual/workflow-tips/), or [notable differences from other languages](https://docs.julialang.org/en/v1/manual/noteworthy-differences/)

To develop a Julia package the package `Revise` (`using Revise`) is helpful, as it allows changes in a file to be applied directly, so that you don't have to restart Julia. Use by `using Revise`.

Many questions are answered in Julias [discourse forum](https://discourse.julialang.org/).



## Installation

__The code should be self-contained and run out of the box.__

First, install [Julia](https://julialang.org/) (we used Julia 1.3).


### Packages

Required packages: `PyPlot`, `JuMP`, `MosekTools`, `PowerModels`, `MAT`, `LinearAlgebra`, `StatsFuns`.
To install packages 
- start Julia
- enter `]` to switch the the package manager. Your console should read `(v1.3) pkg>` (depends on your Julia version)
- enter `add PyPlot JuMP MosekTools PowerModels MAT StatsFuns`.

Additionally, the code requires two self-written Julia packages that need to be registered locally.
Note that this very repository *is* a package, namely the Julia Package `DCsOPF`.

- clone the repo for [PowerModelsParsing](https://iai-vcs.iai.kit.edu/advancedcontrol/code/PowerModelsParsing)  to `<your-path>/PowerModelsParsing`
- clone *this* repo to `<your-path>/DCsOPF`
- start Julia
- enter `]` to switch the the package manager. Your console should read `(v1.3) pkg>` (depends on your Julia version)
- enter `dev <your-path>/PowerModelsParsing`
- enter `dev <your-path>/DCsOPF`.

For more information on adding unregistered packages, see the [Julia package manager documentation](https://julialang.github.io/Pkg.jl/v1/managing-packages/#Adding-unregistered-packages-1).


### Mosek license

The code needs Mosek to solve the optimization problem.
You can obtain a license [here](https://www.mosek.com/license/request/personal-academic/).
The Julia package `MosekTools` will complain if it can't find the license file, and it should tell you where to put it by default.
In case of any issues, consult the [MosekTools documentation](https://juliapackages.com/p/mosektools)



## Run & modify the code


### Run the code
To run the code simply execute the `main.jl` file from top to bottom. 
At the top there are parameters for setup of plots, then there are different cases, and then the main pipeline.


### Modify the code
We can modify the code by adding a new case (`casefiles/case.m`) and by adding new wind data (`data/CovMatrix.mat`).

To add a new case (network) the MatPower casefile `case.m` with the network topology and parameters is needed, as well as the PTDF (power transfer distribution factor) file. The wind data (covariance matrix) can be applied to any case.

| Description | Example | Note |
| --- | --- | --- |
| [Matpower](https://matpower.org/) case file | [`casefiles/case5.m`](examples/casefiles/case5.m) | Has to be of a specific [caseformat](https://matpower.org/docs/ref/matpower5.0/caseformat.html) e.g. [IEEE case5](https://matpower.org/docs/ref/matpower5.0/case5.html) |
| Power transfer distribution factor matrix | [`casefiles/case5_PTDF.mat`](examples/casefiles/case5_PTDF.mat) | The PTDF can be extracted out of `case.m` with MATLAB using `mpc = loadcase('<case>.m'); ptdf = makePTDF(mpc); save('<case>_PTDF.mat'), 'ptdf');`. |
| GP wind data (covariance matrix) | [data/CovMatrix.mat](examples/data/CovMatrix_artificial.mat) | The GP (Gaussian process) is created with GPR (Gaussian process regression), e.g. with the [GPflow](https://github.com/GPflow/GPflow.git) package. |


#### Gaussian processes (covariance matrices)

The wind forecast is given in form of a covariance matrix, e.g. [`CovMatrix_artificial.mat`](examples/data/CovMatrix_artificial.mat), that contains the mean and variance of the forecast. If we forecast over `n` time steps, then the files are:
- `mu_post`: mean (size `nx1`)
- `Lpost`: variance (size `nxn`, lower triangular matrix)
Letting `Σ = Lpost * Lpost'`, then `(mu_post, Σ)` is a Gaussian process over `n` time steps.

#### Add case in main.jl

To test the code with your own network add the parameters for another case with the following code:

```julia
params["<case>"] = CaseParams("<case number>", # case number, e.g. "5"
unc_set, # wind buses (uncertainties), e.g. [[1,2],[1,2,3]]
stor_set, # storage busses, e.g. [[4],[4,5]]
true, # executed local optimization
true, # execute global optimization
false, # solve deterministic problem (covariance matricies are zero)
false, # use artificial wind power curves
1., # multiply load with factor
1., # multiply wind injection with factor
10., # multiply uncertainty (variance) of wind with factor
10, # storage upper bound
0.8) # amplitude of load modelled as sine curve
```

where each load is a dictionary with two fields `:μ` and `:Σ`. This is where you can insert your GP information. Just make sure that the numerical range is about the same. Also note that loads are *negative injections*, hence the minus sign.

### GPR

For the Gaussian process regression we use [GPflow](https://github.com/GPflow/GPflow.git). These files are not currently for display.

### Data

The wind data is taken from [Zenodo](https://zenodo.org/record/4682697#.YksQ2OdCTmG) containing hourly data from 2014 to 2021. 

## Questions?
__Write to us!__ 