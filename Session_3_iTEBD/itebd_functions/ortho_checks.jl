function non_ortho(t::ITensor)
    @assert ndims(t) == 2 
    return (norm(matrix(t)-I(size(t,1))))
end

# Helper: contract a 3-tensor with its conjugate over the bond index `lind`
# (i.e. trace over physical + one link), return the deviation from identity.
function _lenv_dev(phi::ITensor, lind::Index, sitesAB)
    rind = uniqueinds(phi, lind, sitesAB)
    env = phi * replaceind(dag(phi), rind..., sim.(rind)...)
    return non_ortho(env)
end

function _renv_dev(phi::ITensor, lind::Index)
    env = phi * replaceind(dag(phi), lind, sim(lind))
    return non_ortho(env)
end

"""
    right_isometry_dev(T::ITensor, right_ind::Index) -> Float64

Frobenius-norm deviation from right-isometry for a rank-3 tensor `T` with two
virtual indices and one physical index (tagged "Site").

Contracts T·T† over (`right_ind`, phys) and returns ‖result − I‖ on the
remaining (left) virtual index.  Zero iff T is a right isometry.
"""
function right_isometry_dev(T::ITensor, right_ind::Index)
    @assert ndims(T) == 3 "Expected rank-3 tensor, got rank $(ndims(T))"
    phys     = only(inds(T, "Site"))
    left_ind = only(uniqueinds(T, right_ind, phys))
    return _renv_dev(T, left_ind)
end

"""
    vidal_deviations(c::SASBCell) -> NamedTuple

Return the four Frobenius-norm deviations from the single-site Vidal isometry
conditions.  All four are zero iff `c` is in proper Γ–Λ (Vidal) canonical form.

    sA  : ‖ (Λ_BA · Γ_A)† (Λ_BA · Γ_A) − I ‖    (left-isometry, inter-cell BA bond)
    As  : ‖ (Γ_A · Λ_AB)† (Γ_A · Λ_AB) − I ‖    (left-isometry, intra-cell BA bond)
    sB  : ‖ (Λ_AB · Γ_B)† (Λ_AB · Γ_B) − I ‖    (left-isometry, AB bond)
    Bs  : ‖ (Γ_B · Λ_BA)† (Γ_B · Λ_BA) − I ‖    (left-isometry, intra-cell AB bond)
"""
function vidal_deviations(c::SASBCell)
    sitesAB = siteinds(c)

    phi, lind = contract_chain(c,[:sBA, :A])
    δsA = _lenv_dev(phi, lind, sitesAB)

    phi, lind = contract_chain(c,[:A, :sAB])
    δAs = _renv_dev(phi, lind)

    phi, lind = contract_chain(c,[:sAB,:B])
    δsB = _lenv_dev(phi, lind, sitesAB)

    phi, lind = contract_chain(c,[:B,:sBA])
    δBs = _renv_dev(phi, lind)

    return (sA=δsA, As=δAs, sB=δsB, Bs=δBs)
end

"""
    is_vidal_form(c::SASBCell; tol=1e-10) -> Bool

Return `true` if all four single-site Vidal isometry conditions hold within `tol`.
"""
function is_vidal_form(c::SASBCell; tol=1e-9)
    devs = vidal_deviations(c)
    return all(d -> d < tol, devs)
end

function check_ortho(c::SASBCell; tol=1e-9)
    sitesAB = siteinds(c)

    # single-site Vidal isometry deviations
    devs = vidal_deviations(c)

    # two-site deviations
    lenv, lind = contract_chain(c,[:sBA,:A,:sAB,:B])
    lenv = lenv * replaceind(dag(lenv), uniqueinds(lenv, lind, sitesAB)..., sim(uniqueinds(lenv, lind, sitesAB)...))
    δsAsB = non_ortho(lenv)

    renv, lind = contract_chain(c,[:B,:sBA,:A,:sAB])
    renv = renv * replaceind(dag(renv), lind, sim(lind))
    δBsAs = non_ortho(renv)

    lenv, lind = contract_chain(c,[:sAB,:B,:sBA,:A])
    lenv = lenv * replaceind(dag(lenv), uniqueinds(lenv, lind, sitesAB)..., sim(uniqueinds(lenv, lind, sitesAB)...))
    δsBsA = non_ortho(lenv)

    renv, lind = contract_chain(c,[:A,:sAB,:B,:sBA])
    renv = renv * replaceind(dag(renv), lind, sim(lind))
    δAsBs = non_ortho(renv)

    all_devs = (sA=devs.sA, As=devs.As, sB=devs.sB, Bs=devs.Bs,
                sAsB=δsAsB, BsAs=δBsAs, sBsA=δsBsA, AsBs=δAsBs)

    if all(d -> d < tol, all_devs)
        @info "SASBCell ortho ✓  (all deviations < $tol)"
    else
        @warn "SASBCell non-ortho deviations:"
        @info "  [.A]>   = $(devs.sA)"
        @info "  <[A.]   = $(devs.As)"
        @info "  [.B]>   = $(devs.sB)"
        @info "  <[B.]   = $(devs.Bs)"
        @info "  [.A.B]> = $δsAsB"
        @info "  <[B.A.] = $δBsAs"
        @info "  [.B.A]> = $δsBsA"
        @info "  <[A.B.] = $δAsBs"
    end

    return all_devs
end
