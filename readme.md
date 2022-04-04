# Stochastic DC optimal power flow using Gaussian processes

How can we combine stochastic DC optimal power flow with Gaussian processes?

This repository provides both a [paper draft](docs/main.pdf) and a Julia package to tackle this question.

__The code should be self-contained and run out of the box (modulo installing missing/required packages).__

## Developing with Julia

The Julia manual has plenty to offer, e.g. [getting started](https://docs.julialang.org/en/v1/manual/getting-started/), [workflow tips](https://docs.julialang.org/en/v1/manual/workflow-tips/), or [notable differences from other languages](https://docs.julialang.org/en/v1/manual/noteworthy-differences/).

My workflow has become something like this:

- While in the experimental phase of coding, I just code stuff down, and run it in the Julia repl via `include("myfile.jl")`.
- Once the code base is more mature, I put it in a package such that I can load it conveniently via `using MyPackage`. Since Julia is pre-compiled, any change to a Julia package requires a restart of Julia to become effective -- in theory. There is (the amazing) package [`Revise.jl`](https://juliapackages.com/p/revise). This lets you do the following:
    - start Julia
    - enter `using Revise`
    - enter `using MyPackage`
    - try things out and realize that a function in `MyPackage` needs to be changed.
    - apply these changes
    - in the repl you *do not have to quit Julia*, but the changes are effective immediately
- I use VSCodium with the Julia plugin, there is something similar for Atom
- For any questions, I recommend posting them in the [discourse forum](https://discourse.julialang.org/)

## How to run the code?

__The code should be self-contained and run out of the box.__

First, install [Julia](https://julialang.org/).
__This code has been tested only for Julia 1.3.__

### Install packages

Make sure you have the following official packages installed: `PyPlot`, `JuMP`, `MosekTools`, `PowerModels`, `MAT`, `LinearAlgebra`, `StatsFuns`.
To do so

- start Julia,
- enter `]` to switch the the package manager. Your console should read `(v1.3) pkg>` (depends on your Julia version),
- enter `add PyPlot JuMP MosekTools PowerModels MAT StatsFuns`.

Additionally, the code requires two self-written Julia packages that need to be registered locally.
Note that this very repository *is* a package, namely the Julia Package `DCsOPF`.

- clone the repo for [PowerModelsParsing](https://iai-vcs.iai.kit.edu/advancedcontrol/code/PowerModelsParsing)  to `<your-path>/PowerModelsParsing`
- clone *this* repo to `<your-path>/DCsOPF`
- start Julia
- enter `]` to switch the the package manager. Your console should read `(v1.3) pkg>` (depends on your Julia version)
- enter `dev <your-path>/PowerModelsParsing`
- enter `dev <your-path>/DCsOPF`.

__For more information on adding unregistered packages, see the [Julia package manager documentation](https://julialang.github.io/Pkg.jl/v1/managing-packages/#Adding-unregistered-packages-1).__

### Obtain a Mosek license

The code needs Mosek to solve the optimization problem.
You can obtain a license [here](https://www.mosek.com/license/request/personal-academic/).
The Julia package `MosekTools` will complain if it can't find the license file, and it should tell you where to put it by default.
In case of any issues, consult the [MosekTools documentation](https://juliapackages.com/p/mosektools).

### Execute a file

- switch to `<your-path>/DCsOPF/examples`
- start Julia
- enter `include("main_case39_v01.jl")`

Quite likely, you get a warning reading

```julia
julia> include("main_case39_v01.jl")
[ Info: Precompiling DCsOPF [9f832317-d226-4fc3-8cf0-f6d45d985ec6]
┌ Warning: Package DCsOPF does not have PowerModelsParsing in its dependencies:
│ - If you have DCsOPF checked out for development and have
│   added PowerModelsParsing as a dependency but haven't updated your primary
│   environment's manifest file, try `Pkg.resolve()`.
│ - Otherwise you may need to report an issue with DCsOPF
└ Loading PowerModelsParsing into DCsOPF from project dependency, future warnings for DCsOPF are suppressed.
```

This is because our self-written package `DCsOPF` uses the other self-written package `PowerModelsParsing`.
I couldn't yet figure our how to circumvent the warning.

## How to modify the code

As a rule of thumb, I would suggest we don't change the working prototype that much.
Instead, we only change where the GP information comes in.
Still, here's some background information:

There a few file inputs required to set up and solve a problem:

| Description | Origin | Used in |
| --- | --- | --- |
| [Matpower](https://matpower.org/) case file | [`casefiles/case39_v01.m`](examples/casefiles/case39_v01.m) | [`casefiles/case39_v01_buildsystem.jl`](examples/casefiles/case39_v01_buildsystem.jl#L8)
| Power transfer distribution factor matrix | `casefiles/PTDF_case.mat` | [`casefiles/case39_v01_buildsystem.jl`](examples/casefiles/case39_v01_buildsystem.jl#L10) |
| __GP load data__ | [data/LoadData.mat.mat](examples/data/LoadData.mat.mat) | [`casefiles/case39_v01_buildsystem.jl`](examples/casefiles/case39_v01_buildsystem.jl#L123)

Besides this information, there are several other modelling assumptions that have to be made, e.g.

- where to add storages,
- where to add uncertain loads,
- how to choose constraints

This is all taken care of in [`casefiles/case39_v01_buildsystem.jl`](examples/casefiles/case39_v01_buildsystem.jl#L6).

### Gaussian process information

The file [`LoadData.mat`](examples/data/LoadData.mat) contains two data structures

```matlab
  Name           Size              Bytes  Class     Attributes

  Lpost        200x200            320000  double              
  mu_post      200x1                1600  double  
```

where `mu_post` is a vector of `200` values, and `Lpost` is a lower triangular matrix. Letting `Σ = Lpost * Lpost'`, then `(mu_post, Σ)` is a Gaussian process over `200` time steps.

In the code, __only the information about `Lpost` is being used__. You can see this tracing back the few lines [from here](examples/casefiles/case39_v01_buildsystem.jl#L133): The information for `Σ` is taken from the file `LoadData.mat`, but the information about the mean is essentially taken from the function [`my_mean`](examples/casefiles/case39_v01_buildsystem.jl#L74).
This is totally hacky and totally for just playing around -- __it's exactly here where you come in!__

## To-Dos

- if you want to insert your GP information, then I suggest you modify the function [`initializeUncertainty`](examples/casefiles/case39_v01_buildsystem.jl#L78).
- currently, there are a total of 7 uncertain loads,

```julia
julia> mysys.dist[:unc]
Dict{Any,Any} with 7 entries:
  7 => Dict{Symbol,Array{Float64,N} where N}(:μ=>[-2.81, -2.88273, -2.9505, -3.0087, -3.05335, -3.08143, -3.091, -3.08143, -3.05335, -3.0087, -2.9505, -2.88273],:Σ=>[-0.00866991 -0.0 … -0.0 -0.0; -0.0175542 -0.00200176 … -0.0 -0.0; … ; -0.150397 -0.0964019 … -…
  4 => Dict{Symbol,Array{Float64,N} where N}(:μ=>[-6.8, -6.976, -7.14, -7.28083, -7.3889, -7.45683, -7.48, -7.45683, -7.3889, -7.28083, -7.14, -6.976],:Σ=>[-0.00866991 -0.0 … -0.0 -0.0; -0.0175542 -0.00200176 … -0.0 -0.0; … ; -0.150397 -0.0964019 … -0.00027119…
  2 => Dict{Symbol,Array{Float64,N} where N}(:μ=>[-5.22, -5.3551, -5.481, -5.58911, -5.67207, -5.72421, -5.742, -5.72421, -5.67207, -5.58911, -5.481, -5.3551],:Σ=>[-0.00866991 -0.0 … -0.0 -0.0; -0.0175542 -0.00200176 … -0.0 -0.0; … ; -0.150397 -0.0964019 … -0.…
  3 => Dict{Symbol,Array{Float64,N} where N}(:μ=>[-3.29, -3.37515, -3.4545, -3.52264, -3.57492, -3.60779, -3.619, -3.60779, -3.57492, -3.52264, -3.4545, -3.37515],:Σ=>[-0.00866991 -0.0 … -0.0 -0.0; -0.0175542 -0.00200176 … -0.0 -0.0; … ; -0.150397 -0.0964019 ……
  5 => Dict{Symbol,Array{Float64,N} where N}(:μ=>[-2.74, -2.81092, -2.877, -2.93375, -2.97729, -3.00466, -3.014, -3.00466, -2.97729, -2.93375, -2.877, -2.81092],:Σ=>[-0.00866991 -0.0 … -0.0 -0.0; -0.0175542 -0.00200176 … -0.0 -0.0; … ; -0.150397 -0.0964019 … -…
  6 => Dict{Symbol,Array{Float64,N} where N}(:μ=>[-1.39, -1.42598, -1.4595, -1.48829, -1.51038, -1.52426, -1.529, -1.52426, -1.51038, -1.48829, -1.4595, -1.42598],:Σ=>[-0.00866991 -0.0 … -0.0 -0.0; -0.0175542 -0.00200176 … -0.0 -0.0; … ; -0.150397 -0.0964019 ……
  1 => Dict{Symbol,Array{Float64,N} where N}(:μ=>[-5.0, -5.12941, -5.25, -5.35355, -5.43301, -5.48296, -5.5, -5.48296, -5.43301, -5.35355, -5.25, -5.12941],:Σ=>[-0.00866991 -0.0 … -0.0 -0.0; -0.0175542 -0.00200176 … -0.0 -0.0; … ; -0.150397 -0.0964019 … -0.000…
```

where each load is a dictionary with two fields `:μ` and `:Σ`. This is where you can insert your GP information. Just make sure that the numerical range is about the same. Also note that loads are *negative injections*, hence the minus sign.


# GP Regression (addition)

To the existing repository have been added files that calculate the LoadData.mat files (mean and covariance). The file is:

* main_pipeline.ipynb (and main_pipeline.py)

The python script is derived via the jupyter notebook via this command:

`jupyter nbconvert --to script main_pipeline.ipynb --template=simplepython.tpl --RegexRemovePreprocessor.patterns="['# hide']"`

The template ensures that only code cells are exported and the preprocessor ensure that code cells containing `# hide` are not exported. This way we get a clean python script.

## Repository
The original repository is:
* [GP Regression](https://iai-vcs.iai.kit.edu/uydye/gpregression)

This repository implements GPR on a dataset of the Reduced Turkish High Voltage Transmission Network.

The packages are `GPflow` and `GPyTorch`. We use both for comparison, will later probably move to __GPyTorch__. The code is leaned on the respective repositories' examples:
* [GPflow](https://github.com/GPflow/GPflow.git) (Tutorial: [regression.pct.py](https://github.com/GPflow/docs/blob/master/doc/source/notebooks/basics/regression.ipynb))
* [GPyTorch](https://github.com/cornellius-gp/gpytorch.git) (Tutorial: [Simple_GP_Regression.ipynb](https://github.com/cornellius-gp/gpytorch/blob/master/examples/01_Exact_GPs/Simple_GP_Regression.ipynb))

The code is in __Jupyter Notebooks__ only. 

### Packages

We mainly need the two main packages:
* [GPflow](https://github.com/GPflow/GPflow.git)
* [GPyTorch](https://github.com/cornellius-gp/gpytorch.git)
of which GPflow uses `Tensorflow`.

For optional data visualisation in Jupyter Lab we use `ipywidgets`.

### Data

The data is from the reduced Turkish high voltage transmision system with 116 nodes and 190 transmission lines. It contains hourly demand measurements of 2016.

![TurkishTransmissionNetwork](src/data/Grid.png)

The data is provided with the following paper and can be downloaded from the same page: 
* [A real data set for a 116-node power transmission system](https://data.mendeley.com/datasets/dv3vjnwwf9/1)
This data (`Data.rar`) contains also many information about the grid. The relevant file is `Counties.xlsx` containing the sheet `hourly demand of conties`. From that sheet there is a short excerpt in this repository with which the code runs:
* [Test data from ADANA Aladag](src/data/Test_County_ADANA_Aladag.xlsx)

### Results (Test Data)

We tried both packages on the test data set. Becaues GPR does not work for large data sets we sampled few (e.g. 10, not yet firmly representative) points and performed GPR on it. The resulst vary grandly as the number of points chosen is relatively small.

Test data:

![Test Data Sampled](pics/test_data.png)

#### GPflow

Performing GPR with GPflow on the test data gives:

![Test Results GPflow](pics/results_test_gpflow.png) 

#### GPyTorch

Performing GPR with GPyTorch on the test data gives:

![Test Results GPyTorch](pics/results_test_gpytorch.png)