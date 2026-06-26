
# ── high-level loop ───────────────────────────────────────────────────────────
#
# Applies the gate sequence `gates[1..K]` for each of `Nt` physical timesteps,
# alternating step_AB / step_BA by position in the sequence (odd → AB, even → BA).
# The observer is called once per physical step with keyword arguments `c` and `step`.
#
# Single gate (backwards compat):  wrapped as [gate, gate]  →  AB + BA per step.
#
# Examples:
#   tebd(Nt, [gate, gate],                     obs!, c0)  # 1st-order
#   tebd(Nt, [gate_half, gate_full, gate_half], obs!, c0)  # 2nd-order Strang

function itebd_vidal(Nt::Int, gate::ITensor, observer!, c0::SASBCell; kwargs...)
    itebd_vidal(Nt, [gate, gate], observer!, c0; kwargs...)
end

function itebd(Nt::Int, gates::AbstractVector{<:ITensor}, observer!, c::SASBCell; cutoff=1e-8, maxdim=128)
    sites = siteinds(c)
    
    for step = 1:Nt
        for (k, gate) in enumerate(gates)
            if isodd(k)
                c = step_AB(c, gate; cutoff, maxdim)
            else
                c = step_BA(c, gate; cutoff, maxdim)
            end
        end
        # Pass the same kwargs as tebd_hastings so observers are interchangeable
        phi, _, _ = contract_chain(c, [:sBA, :A, :sAB, :B])
        workS = c.sAB
        Observers.update!(observer!; phi, workS, step)
    end
    return c
end


"""
    _standardize_cell(A, sAB, B, sBA) -> (A, sAB, B, sBA)

Relabel the four shared bond indices of a freshly-built cell to the standard
SASBCell tag convention (`Link,sA | Link,As | Link,sB | Link,Bs`). The SVD
inside the step generates link indices with default tags, so this keeps the
cell tags consistent across steps (and keeps `check_inds_sasbcell` happy).
"""
function _standardize_cell(A::ITensor, sAB::ITensor, B::ITensor, sBA::ITensor)
    function retag(t1, t2, tag)
        old = commonind(t1, t2)
        new = settags(old, tag)
        return replaceind(t1, old => new), replaceind(t2, old => new)
    end
    sBA, A   = retag(sBA, A,   "Link,sA")
    A,   sAB = retag(A,   sAB, "Link,As")
    sAB, B   = retag(sAB, B,   "Link,sB")
    B,   sBA = retag(B,   sBA, "Link,Bs")
    return A, sAB, B, sBA
end

function itebd_vidal(Nt::Int, gates::AbstractVector{<:ITensor}, observer!,
    c::SASBCell;
    cutoff=1e-8, maxdim=128, ortho_every::Int=0)
    sites = siteinds(c)

    for step in 1:Nt

        # read the (evolving) cell into work tensors at the start of each step
        s_ext = c.sBA
        workS = c.sAB
        workL = c.A
        workR = c.B

        for gate in gates
            # One-body gates are trivial
            if ndims(gate) <= 2
                if hascommoninds(workL, gate)
                    workL = apply(gate, workL)
                elseif hascommoninds(workR, gate)
                    workR = apply(gate, workR)
                else
                    error("Don't know where to apply this gate? $(inds(gate))")
                end

            else # for 2body we truncate (and swap)

                phi, lind = build_phi(s_ext, workL, workS, workR)

                thetaAB = apply(gate, phi)

                siteL = only(commoninds(sites, workL))

                leftinds = IndexSet(lind, siteL)

                u,s_new,vd = svd(thetaAB, leftinds;
                cutoff, maxdim) # TODO adapt SVD tags , lefttags, righttags)

                s_new = normalize(s_new)
                s_inv = inv.(s_ext)

                workR = replaceind(s_inv, commonind(s_inv, workR) => lind) * u
                workL = vd * s_inv

                workS = s_ext
                s_ext = noprime(s_new)

                @assert ndims(workL) == 3 "$(inds(workL))"
                @assert ndims(workR) == 3 "$(inds(workR))"
                @assert ndims(workS) == 2 "$(inds(workS))"
            end
        end

        # After an even (AB,BA) gate sequence the work tensors map back to a
        # standard-layout cell:  A←workL, sAB←workS, B←workR, sBA←s_ext.
        A, sAB, B, sBA = _standardize_cell(workL, workS, workR, s_ext)
        c = SASBCell(A, sAB, B, sBA)

        # periodic full re-orthogonalization of the unit cell
        if ortho_every > 0 && step % ortho_every == 0
            c = orthogonalize_cell(c)
        end

        phi, _ = build_phi(c.sBA, c.A, c.sAB, c.B)
        Observers.update!(observer!; phi, workS = c.sAB, step)
    end

    return c
end


# ── single bond steps ─────────────────────────────────────────────────────────

function step_AB(c, hAB; cutoff=1e-8, maxdim=500)

    @info "AB update"

    #  -Bs - sasbs - sA 
    sasbs, lind, _ = contract_chain(c, [:sBA, :A, :sAB, :B, :sBA])

    thetaAB = apply(hAB, sasbs)

    # Update two-site WF
    #         [pA]   [pB]
    # [Bs]-sBA--A-sAB-B--sBA-[sA]

    s_ext = c.sBA
    u,s_new,vd = svd(thetaAB, (lind,siteind(c,"A")); cutoff, maxdim, lefttags="As,Link", righttags="sB,Link")

    s_new = normalize(s_new)
    s_inv = inv.(s_ext)
    
    Anew = replaceind(s_inv, commonind(s_inv, c.B), lind) * u
    Bnew = vd * s_inv

    # Update 
    return SASBCell(Anew, s_new, Bnew, s_ext)

end

function step_BA(c, hBA; cutoff=1e-20, maxdim=500)

    @info "BA update"

    #          [pB]  [pA]
    # [As]-sAB--B-sBA-A--sAB-[sB]
    sbsas, lind, _ = contract_chain(c, [:sAB, :B, :sBA, :A, :sAB])

    thetaBA = apply(hBA, sbsas)
    u,s_new,vd = svd(thetaBA, (lind, siteind(c,"B")); cutoff, maxdim,lefttags="Bs,Link", righttags="sA,Link")

    s_new = normalize(s_new)

    s_ext = c.sAB
    s_inv = inv.(s_ext)
    Bnew = u * replaceind(s_inv, commonind(s_inv, c.A), lind)
    Anew = vd * s_inv
    #@info "chi Snew = $(dims(s_new))"
    return SASBCell(Anew, s_ext, Bnew, s_new)
end
