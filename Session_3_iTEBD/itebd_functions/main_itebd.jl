using LinearAlgebra
using ITensors, ITensorMPS
using KrylovKit
using ITransverse

import Observers
using Observers: observer

# extended for SASBCell in unit_cell.jl
import ITensors: siteinds, siteind
import ITensorMPS: linkinds

using Plots

include("expH_ising.jl")
include("unit_cell.jl")
include("cell_checks.jl")
include("contract_cell.jl")
include("pretty_prints.jl")
include("environments.jl")
include("ortho_checks.jl")
include("tebd.jl")
include("expvals.jl")


function main_tebd(Nt::Int; ortho_every::Int=10, dt::Number=0.1)

    cutoff = 1e-10
    maxdim = 128

    # IsingParams is a struct that take Jxx, gperp and hpara and for the 
    # Ising model then this can be used to inicialize the function the default are (1.,0.4,0.)
    mp = IsingParams()
    init_state = up_state
    
    c = init_c(init_state)

    sitesAB = siteinds(c)

    hAB = exphAB_ising(mp, dt, sitesAB)

    Za = op(sitesAB[1], "Z")

    obs = observer(
        "energy"  => (; phi)   -> real(scalar(apply(hAB, phi) * dag(phi))),
        "zeta"    => (; phi)   -> real(scalar(apply(Za, phi) * dag(phi))),
        "norm"    => (; phi)   -> real(scalar(phi * dag(phi))),
        "entropy" => (; workS) -> entanglement_entropy(workS),
        "chi"     => (; workS) -> dim(workS, 1)
    )

    c = itebd_vidal(Nt, [hAB, hAB], obs, c; cutoff, maxdim, ortho_every)

    return (; c, obs)
end

dt = 0.01
TT = 500
Tf = TT*dt
resu_itebd = main_tebd(TT; dt)

# Check the evolved cell is consistent / in canonical form
c = resu_itebd.c
check_inds_sasbcell(c)
check_ortho(c)

@show entanglement_entropy(c, "AB")
@show compute_expvals(c, "Z")

let
plot((dt:dt:Tf),resu_itebd.obs.zeta)
xlabel!("Time")
ylabel!("<Z>")
end
let
plot((dt:dt:Tf),resu_itebd.obs.chi)
xlabel!("Time")
ylabel!("χ")
end
let
    plot((dt:dt:Tf),resu_itebd.obs.entropy, label="S_VN")
    plot!((dt:dt:Tf),0.25*(dt:dt:Tf).+0.2, label="Volume law")
    
    xlabel!("Time")
    ylabel!("S_VN")
end

let
    plot((dt:dt:Tf),resu_itebd.obs.energy.-resu_itebd.obs.energy[1])
    xlabel!("Time")
    ylabel!("ΔE")
end

let
plot((dt:dt:Tf),resu_itebd.obs.entropy, label="S_VN")
plot!((dt:dt:Tf),0.25*(dt:dt:Tf).+0.2, label="Volume law")

xlabel!("Time")
ylabel!("S_VN")
end