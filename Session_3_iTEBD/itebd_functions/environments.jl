""" Given a transfer matrix and combiners, builds its dominant 
eigenvector and diagonalizes it. Returns eigen struct F and vectorized environment """
function diag_env(tm::ITensor, cL::ITensor, cR::ITensor; 
         direction::Symbol=:right, normalize::Bool=true, kwargs...)

    (dir_ind, open_ind) = direction == :left ? (cL, cR) : (cR, cL)

    # dominant vector
    dlam, dvec, info = dominant_eigenvectors(tm, combinedind(dir_ind); ishermitian=false)
   
    @assert only(inds(dvec[1])) == combinedind(open_ind) " ? $(only(inds(dvec[1]))) != $(open_ind)"

    dv = dvec[1] * dag(open_ind)

    if normalize
        dv = dv / tr(matrix(dv))
    end

    # at least for hermitian, F.Vt and F.V are the same right eigenvector, Vt has the "left" index 
    # F.Vt * F.D * dag(F.V) ≈ m    #true
    ITransverse.isapproxherm(matrix(dv)) || @warn "diag_env: dominant environment not hermitian"
    F = eigen(dv, inds(dv)...; ishermitian=true, cutoff=1e-8, kwargs...)

    return F, dv
end

function diag_env(cell, pieces::AbstractVector{Symbol}; kwargs...)
    diag_env(compute_tm(cell, pieces)...; kwargs...)
end

function sqrt_env(cell, pieces::AbstractVector{Symbol}; direction::Symbol=:right, kwargs...)

    F, dv = diag_env(cell, pieces; direction, kwargs...)
    X = sqrt(complex(F.D)) * F.Vt
    Xinv = complex(F.D).^(-0.5) * dag(F.Vt)

    # debugging 
    Xinv_dag =  complex(F.D).^(-0.5) * (F.V)

    #@show inds(dv)
    #@show inds(Xinv)
    if !isid(dv * Xinv * Xinv_dag)
        @warn "Not good inverse ?"
    end

    return X, Xinv
end


function orthogonalize_cell(c::SASBCell; kwargs...)
   
    # Center in AB 

    # Left 
    X, Xinv = noprime.(sqrt_env(c, [:sAB, :B, :sBA, :A], direction=:left))
    # Right
    Y, Yinv = noprime.(sqrt_env(c, [:B, :sBA, :A, :sAB]))

    new_S =  X * c.sAB * Y
    @assert ndims(new_S) == 2 
    U,new_sAB,Vd = svd(new_S, commonind(new_S, Xinv); lefttags="Link,As", righttags="Link,sB") # todo check no trunc here ? 

    Xinv *= U
    Yinv *= Vd

    new_A = c.A * Xinv
    new_B = Yinv * c.B

    c = SASBCell(new_A, normalize(new_sAB), new_B, c.sBA)
    
    # Center in BA
    
    # Left 
    X, Xinv = noprime.(sqrt_env(c, [:sBA, :A, :sAB, :B], direction=:left))
    # Right
    Y, Yinv = noprime.(sqrt_env(c, [:A, :sAB, :B, :sBA]))

    new_S =  X * c.sBA * Y
    U,new_sBA,Vd = svd(new_S, commonind(new_S, Xinv); lefttags="Link,Bs", righttags="Link,sA") # todo check no trunc here ? 

    Xinv *= U
    Yinv *= Vd

    new_B = c.B * Xinv
    new_A = c.A * Yinv

    c = SASBCell(new_A, c.sAB, new_B, normalize(new_sBA))

    return normalize_cell(c)

end

function normalize_cell(c::SASBCell)

    sAB = normalize(c.sAB)
    sBA = normalize(c.sBA)

    sAs = sBA * c.A * sAB
    sBs = sAB * c.B * sBA

    An = c.A/norm(sAs)
    Bn = c.B/norm(sBs)

    return SASBCell(An, sAB, Bn, sBA)

end

