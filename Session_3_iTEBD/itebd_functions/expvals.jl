
function compute_expvals(c::SASBCell, obs = "X")

    sitesAB = siteinds(c)

    sas, _, _ = contract_chain(c, [:sBA, :A, :sAB])
    evsas = sas * op(sitesAB[1], obs)
    evsas *= dag(sas) * delta( sitesAB[1], sitesAB[1]')
    evsas = scalar(evsas)

    sbs, _, _ = contract_chain(c, [:sAB, :B, :sBA])
    evsbs = sbs * op(sitesAB[2], obs)
    evsbs *= dag(sbs) * delta( sitesAB[2], sitesAB[2]')
    evsbs = scalar(evsbs)


    sasbs, _, _ = contract_chain(c, [:sBA, :A, :sAB, :B, :sBA])
    evAB = sasbs * op(sitesAB[1], obs)
    evAB = evAB * dag(sasbs) * delta( sitesAB[1], sitesAB[1]')
    evAB = scalar(evAB)

    if abs(evsas - evAB) > 1e-8
        @warn "different evs? check can form "
    elseif abs(evsas - evsbs) > 1e-8
        @warn "different evs? check trasl inv"
    end

    # TODO return one, all, mean ? 
    return evAB

end

"""
rdm_1site(c, which="A") assuming we are in canonical form,
builds 1-site reduced density matrix for site A or B.
"""
function rdm_1site(c::SASBCell, which::String="A")
    if which == "A" || which == "a"
        phi, _, _ = contract_chain(c, [:sBA, :A, :sAB])
        s = siteind(c, "A")
    else
        phi, _, _ = contract_chain(c, [:sAB, :B, :sBA])
        s = siteind(c, "B")
    end
    return phi * prime(dag(phi), s)   # trace over bond indices, keep s open
end

"""
2-site reduced density matrix for the A-B unit cell, assuming canonical form 
"""
function rdm_2site(c::SASBCell)
    phi, _, _ = contract_chain(c, [:sBA, :A, :sAB, :B])
    sA, sB = siteinds(c)
    return phi * prime(dag(phi), sA, sB)
end

""" 
Computes RDM by partial tracing everything to the right of site `which` (=A,B)
"""
function rdm_left(c::SASBCell, which="A")
    phi, lind, rind = if which == "A"
        contract_chain(c, [:sBA, :A, :sAB, :B, :sBA])
    else
        contract_chain(c, [:sBA, :A, :sAB, :B, :sBA])
    end

    return phi * prime(prime(dag(phi), which), lind), combiner(lind, inds(phi, which))

end



"""
    entanglement_entropy(c, bond="AB")

Von Neumann entropy S = -Σᵢ λᵢ² log λᵢ² across the `bond` cut ("AB" or "BA").

In Vidal form the singular values stored in sAB / sBA ARE the Schmidt coefficients,
so this is O(χ) with no tensor contraction needed.
"""
function entanglement_entropy(c::SASBCell, bond::String="AB")
    s_tensor = bond == "AB" ? c.sAB : c.sBA
    entanglement_entropy(s_tensor)
end

function entanglement_entropy(s_tensor::ITensor)
    λ = diag(matrix(s_tensor))
    λ2 = λ .^ 2
    λ2 ./= sum(λ2)   # guard against numerical drift from imperfect normalization
    return -sum(x -> x > 1e-15 ? x * log(x) : 0.0, λ2)
end
