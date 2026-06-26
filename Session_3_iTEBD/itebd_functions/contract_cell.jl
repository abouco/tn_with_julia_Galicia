
# ---------------------------------------------------------------------------
# General chain contraction
# ---------------------------------------------------------------------------

"""
    _piece(c::SASBCell, p::Symbol) -> (tensor, left_ind, right_ind)

Return the ITensor and its left/right bond indices for piece `p` of the
periodic chain (physical left-to-right order):

    ... —sBA— A —sAB— B —sBA— A —sAB— B— ...

Valid symbols: `:sBA`, `:A`, `:sAB`, `:B`.
"""
function _cell_piece(c::SASBCell, p::Symbol)
    if p === :sBA
        return c.sBA, commonind(c.B, c.sBA), commonind(c.sBA, c.A)
    elseif p === :A
        return c.A,   commonind(c.sBA, c.A),  commonind(c.A, c.sAB)
    elseif p === :sAB
        return c.sAB, commonind(c.A, c.sAB),  commonind(c.sAB, c.B)
    elseif p === :B
        return c.B,   commonind(c.sAB, c.B),  commonind(c.B, c.sBA)
    else
        error("Unknown piece: $p. Expected one of: :sBA, :A, :sAB, :B")
    end
end

"""
contract_chain(c::SASBCell, pieces) -> (phi, left_ind)
Convention: we prime the leftmost ind of the first tensor of the list
"""
function contract_chain(c::SASBCell, pieces::AbstractVector{Symbol})

    @assert 0 < length(pieces) < 6 "bad list of elements"

    _, _, r = _cell_piece(c, pieces[end])
    t1, l, _ = _cell_piece(c, pieces[1])

    phi  = prime(t1, l)
    for p in pieces[2:end]
        phi *= first(_cell_piece(c, p))
    end

    return phi, l', r

end

function build_phi(s_ext::ITensor, workL::ITensor, workS::ITensor, workR::ITensor)
    lind = only(uniqueinds(s_ext, workL))
    phi = prime(s_ext, lind) * workL 
    phi *= workS
    phi *= workR
    phi *= s_ext

    return phi, lind'
end


function compute_tm(ket::ITensor, leftind::Index, rightind::Index)

    @assert hasind(ket, leftind)
    @assert hasind(ket, rightind)

    bra = prime(dag(ket), 2, (leftind, rightind))

    tm = ket*bra 

    cL = combiner(leftind,  leftind'',  tags="left")
    cR = combiner(rightind, rightind'', tags="right")

    tm *= cL 
    tm *= cR

    return tm, cL, cR
end

function compute_tm(cell, pieces::AbstractVector{Symbol})
    compute_tm(contract_chain(cell, pieces)...)
end
