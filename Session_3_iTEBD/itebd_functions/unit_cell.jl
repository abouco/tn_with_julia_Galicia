
""" Unit cell (Gamma-Lambda form): 

          |       |            |       |
[...][--sBA--A--sAB--B--][--sBA--A--sAB--B--][...]

and we have
      |           |                        |            |
--sBA-A-  =   -|ALeft>--    and      -sAB--B-  =   --<Bright|-

 |               |                     |                 |
-A--sAB-  =   -<ARight|--    and      -B-sBA--  =   --|Bleft>-

Alternatively: Orthogonal form, eg. Left everywhere

       |       |           |       |
--sBA--A--sAB--B-- =   -|ALeft>--|Bleft>-

But more interesting for us: mixed form so we are in ortho centre

                 |       |           |       |           |       |
    [...][--sBA--A--sAB--B--][--sBA--A--sAB--B--][--sBA--A--sAB--B--][...]


            |       |        |       |            |       |
 =  [...][--|AL>--|BL>--][--|AL>--|BL>--][-sBA][--A--sAB--B--][...]


            |       |        |       |            |       |
 =  [...][--|AL>--|BL>--][--|AL>--|BL>--][-sBA][--<AR--sAB--B--][...]

"""

""" Unit cell, A sAB B sBA """
struct SASBCell
    A::ITensor
    sAB::ITensor
    B::ITensor
    sBA::ITensor

    function SASBCell(A::ITensor, sAB::ITensor, B::ITensor, sBA::ITensor)
        @assert ndims(A) == 3 "A? $(inds(A))"
        @assert ndims(B) == 3 "B? $(inds(B))"
        @assert ndims(sAB) == 2 "sAB? $(inds(sAB))"
        @assert ndims(sBA) == 2 "sBA? $(inds(sBA))"
        @assert isdiag(sAB) "sAB not diag?"
        @assert isdiag(sBA) "sBA not diag?"

        new(A,sAB,B,sBA)
    end
end

"""
Build a `SASBCell` from two tensors that share link indices in a periodic structure:

    [...] --i_BA-- A --i_AB-- B --i_BA-- [...]

`i_AB`: the shared index sitting between A (right) and B (left)  → becomes the sAB bond.
`i_BA`: the shared index sitting between B (right) and A (left)  → becomes the sBA bond.

Identity singular-value tensors (Λ = I) are inserted at both bonds, and all
link indices are relabeled to the standard SASBCell tag convention:

    Link,sA  (left of A)  |  Link,As  (right of A)
    Link,sB  (left of B)  |  Link,Bs  (right of B)
"""
function SASBCell(A::ITensor, B::ITensor, i_AB::Index, i_BA::Index)
    chi_AB = dim(i_AB)
    chi_BA = dim(i_BA)

    # New tagged indices
    link_as = Index(chi_AB, "Link,As")   # right of A / left  of sAB
    link_sb = Index(chi_AB, "Link,sB")   # right of sAB / left of B
    link_bs = Index(chi_BA, "Link,Bs")   # right of B / left  of sBA
    link_sa = Index(chi_BA, "Link,sA")   # right of sBA / left of A

    newA = replaceind(replaceind(A, i_AB, link_as), i_BA, link_sa)
    newB = replaceind(replaceind(B, i_AB, link_sb), i_BA, link_bs)

    sAB = diagITensor(ones(chi_AB), link_as, link_sb)
    sBA = diagITensor(ones(chi_BA), link_bs, link_sa)

    return SASBCell(newA, sAB, newB, sBA)
end

function linkinds(c::SASBCell)
    # As, Bs, sB, sA
    return commonind(c.A, c.sAB), commonind(c.sBA, c.B), commonind(c.sAB, c.B), commonind(c.A, c.sBA)
end


function siteinds(c::SASBCell)
    return uniqueind(c.A, c.sAB, c.sBA ), uniqueind(c.B, c.sAB, c.sBA)
end


""" TODO This is not the best.. """
function ITensors.siteind(c::SASBCell, which_ab::String)
    si =  which_ab == "A" || which_ab == "a" ? uniqueind(c.A, c.sAB, c.sBA ) : uniqueind(c.B, c.sAB, c.sBA)
    return si
end


function init_c(init_vec::Vector=ComplexF64[1,0])
    link_as = Index(1,"Link,As")
    link_sb = Index(1,"Link,sB")
    link_bs = Index(1,"Link,Bs")
    link_sa = Index(1,"Link,sA")

    phys_a = Index(2, "Site,S=1/2,A")
    phys_b = Index(2, "Site,S=1/2,B")

    # A - sAB - B 
    sAB = diagITensor(1., link_as, link_sb)
    sBA = diagITensor(1., link_bs, link_sa)
    A = ITensor(init_vec, (link_sa, phys_a, link_as))
    B = ITensor(init_vec, (link_sb, phys_b, link_bs))

    c = SASBCell(A, sAB, B, sBA)

    return c
end

"""
    random_cell(; chi=10, d=2) -> SASBCell

Construct a `SASBCell` with random Γ tensors and identity Λ weights.
Indices follow the standard tag convention:

    sBA(Link,sA|Link,Bs) — A(Link,sA|Site|Link,As) — sAB(Link,As|Link,sB) — B(Link,sB|Site|Link,Bs)
"""
function random_cell(; chi1::Int=10, chi2::Int=11, d::Int=2)
    link_sa = Index(chi1, "Link,sA")
    link_as = Index(chi2, "Link,As")
    link_sb = Index(chi2, "Link,sB")
    link_bs = Index(chi1, "Link,Bs")

    phys_a = Index(d, "Site,S=1/2,A")
    phys_b = Index(d, "Site,S=1/2,B")

    A   = random_itensor(ComplexF64, link_sa, phys_a, link_as)
    B   = random_itensor(ComplexF64, link_sb, phys_b, link_bs)
    sAB = normalize(diagITensor(ones(chi2), link_as, link_sb))
    sBA = normalize(diagITensor(ones(chi1), link_bs, link_sa))

    return SASBCell(A, sAB, B, sBA)
end
