using ITensors, ITensorMPS
using JLD2
using Plots

using ITransverse

function main_cone()

    JXX = 1.0  
    hz = 0.7
    gx = 0.0 # 0.5

    dt = 0.1

    nbeta = 0

    optimize_op = vX
    init_state = up_state

    cutoff = 1e-10
    maxdim = 128
    direction = :left
    alg = "densitymatrix"

    truncp = (;cutoff, maxdim, direction, alg)

    Nsteps = 80

    #time_sites = siteinds("S=3/2", 1)

    mp = IsingParams()
    #tp = tMPOParams(dt, expH_ising_murg, mp, nbeta, init_state, Id)
    tp = tMPOParams(mp; dt, scheme=SymSVD(), nbeta, init_state)
    b = FoldtMPOBlocks(tp)
    c0 = init_cone(b)

    #@info length(c0)

    cone_params = ConeParams(;truncp, opt_method=:sym, optimize_op)

    cp = DoCheckpoint(
        "cp_cone_xx.jld2";
        params=Dict("tparams" => tp, "cparams" => cone_params),
        save_at = [],
        f_obs = (
            SVN = s -> vn_entanglement_entropy(s.R),
            S_SVD = s -> generalized_svd_vn_entropy(s.L, s.R),
            overlap = s -> overlap_noconj(s.L, s.R),
            expvals = s -> ITransverse.compute_expvals(s.L, s.R, ["Z","X"], s.b),
            chis = s -> maxlinkdim(s.R)
        ),
        f_savestate = (
            L = s -> s.L,
            R = s -> s.R,
            b = s -> s.b
        )
    )

    psi, psiR, checkpt = run_cone(c0, b, cone_params, cp, Nsteps)

    return  psi, psiR, checkpt

end

psi, psiR, checkpt_cone = main_cone()
