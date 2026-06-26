
"""
    check_inds_sasbcell(c::SASBCell)

Verify that the indices of a SASBCell have the expected tag structure and
are properly shared between adjacent tensors.

Expected layout:  ... --sBA-- A --sAB-- B --sBA-- ...

    A   : Link,sA  (left, shared with sBA) | Site,A  | Link,As (right, shared with sAB)
    sAB : Link,As  (left, shared with A)   | Link,sB (right, shared with B)
    B   : Link,sB  (left, shared with sAB) | Site,B  | Link,Bs (right, shared with sBA)
    sBA : Link,Bs  (left, shared with B)   | Link,sA (right, shared with A)
"""
function check_inds_sasbcell(c::SASBCell)
    (;A, sAB, B, sBA) = c
    ok = true
    
    function _check_has(tensor, tag)
        idx = only(inds(tensor, tags=tag))  # throws if not exactly 1
        return idx
    end

    # --- tag presence ---
    try
        idx_sA_on_A   = _check_has(A, "Link,sA")
        idx_As_on_A   = _check_has(A, "Link,As")
                        _check_has(A, "Site")

        idx_As_on_sAB = _check_has(sAB, "Link,As")
        idx_sB_on_sAB = _check_has(sAB, "Link,sB")

        idx_sB_on_B   = _check_has(B, "Link,sB")
        idx_Bs_on_B   = _check_has(B, "Link,Bs")
                        _check_has(B, "Site")

        idx_Bs_on_sBA = _check_has(sBA, "Link,Bs")
        idx_sA_on_sBA = _check_has(sBA, "Link,sA")

        # --- shared (same Index object) between neighbours ---
        idx_As_on_A   == idx_As_on_sAB || (ok = false; @warn "Link,As: index on A ≠ index on sAB")
        idx_sB_on_sAB == idx_sB_on_B   || (ok = false; @warn "Link,sB: index on sAB ≠ index on B")
        idx_Bs_on_B   == idx_Bs_on_sBA || (ok = false; @warn "Link,Bs: index on B ≠ index on sBA")
        idx_sA_on_sBA == idx_sA_on_A   || (ok = false; @warn "Link,sA: index on sBA ≠ index on A")

    catch e
        @error "check_inds_sasbcell failed: $e"
        return false
    end

    ok && @info "SASBCell index tags OK"
    return ok
end

