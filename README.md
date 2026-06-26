# tn_with_julia_Galicia
Notes based on previous notes from @lukyluket and Glen Evenbly (tensors.net) and are transcripted into julia by @abouco and @starsfordummies. They follow a three part tutorial: 
1. Linear Algebra + Tensor contraction. 
2. Basic tensor networks: MPS + iMPS. 
3. iTEBD

## Basics to start with Jullia

Setup: Installing Julia. Recommended way is via juliaup: https://julialang.org/downloads/

tl;dr is 

```
curl -fsSL https://install.julialang.org | sh
```

then close/reopen shell 

After starting julia ( `julia` from console ), you can enter `Pkg` mode to install packages by pressing `]` 

The packages required to install are, and below the code you need to write to install them:

1. First of all enter into `Pkg` mode by pressing `]` or type into the terminal `using Pkg`
2. You can then install all the Packages required by typing `add ...` where `...` corresponds to the name of the package if you are in the `Pkg` mode, or if you have typed `using Pkg` into the julia terminar you can type `Pkg.add(...)`.
3. Below there is an exhaustive list of all the packages required to run these notebooks. Remember you also need to have installed jupyter notebooks. 

- IJulia.jl -> This package will allow you to run the julia kernel in a jupyter notebook.

If you are in the `Pkg` mode you need to write `add IJulia` else you need to write `Pkg.add("IJulia")`

- Plots.jl -> We will use it to plot our results

If you are in the `Pkg` mode you need to write `add Plots` else you need to write `Pkg.add("Plots")`

- TensorOperations.jl -> we will mainly use its function `ncon`

If you are in the `Pkg` mode you need to write `add TensorOperations` else you need to write `Pkg.add("TensorOperations")`

- OMEinsum.jl -> we will use `ein` and the optimization of tensor contraction path

If you are in the `Pkg` mode you need to write `add OMEinsum` else you need to write `Pkg.add("OMEinsum")`

ITensors.jl -> Package that makes your life easier with finite MPS

If you are in the `Pkg` mode you need to write `add ITensors` else you need to write `Pkg.add("ITensors")`

ITensorMPS.jl -> Package that makes your life easier with finite MPS

If you are in the `Pkg` mode you need to write `add ITensorMPS` else you need to write `Pkg.add("ITensorMPS")`

ITransverse.jl -> Package that makes your life easier with the transverse contraction

You need to write `Pkg.add(url="https://github.com/starsfordummies/ITransverse.jl.git")`

Observers.jl

If you are in the `Pkg` mode you need to write `add Observers` else you need to write `Pkg.add("Observers")`
