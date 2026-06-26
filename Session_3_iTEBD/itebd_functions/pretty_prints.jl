# -- pretty print SASBCell


function _SASBCell_summary(io::IO, c::SASBCell)
    bond_l  = commonind(c.sBA, c.A)    # Link,sA  (left outer)
    bond_as = commonind(c.A,   c.sAB)  # Link,As  (A → sAB)
    bond_sb = commonind(c.sAB, c.B)    # Link,sB  (sAB → B)
    bond_bs = commonind(c.B,   c.sBA)  # Link,Bs  (B → sBA)

    # short label: strip "Link," prefix, append dim
    lbl(i) = replace(string(tags(i)), "Link," => "") * "($(dim(i)))"
    lbls = [lbl(bond_l), lbl(bond_as), lbl(bond_sb), lbl(bond_bs), lbl(bond_l)]

    # wire segment width: wide enough to center the longest label with 1 char to spare
    seg = max(5, maximum(length, lbls) + 1)
    w   = "-" ^ seg

    # Build row2 left-to-right, tracking positions as we go
    nodes = ["A", "sAB", "B", "sBA"]
    row2 = "  "
    col = 3  # 1-indexed position of next char

    wire_starts  = Int[]   # start col of each wire segment (5 total)
    node_centers = Int[]   # center col of each node (for site legs)

    for nd in nodes
        push!(wire_starts, col)
        row2 *= w
        col  += seg
        push!(node_centers, col + div(length(nd) - 1, 2))
        row2 *= nd
        col  += length(nd)
    end
    push!(wire_starts, col)   # trailing wire
    row2 *= w

    # midpoint of each wire segment
    mids = [ws + div(seg, 2) for ws in wire_starts]

    N    = length(row2) + 4
    row1 = fill(' ', N)   # site legs
    row3 = fill(' ', N)   # bond labels

    row1[node_centers[1]] = '|'   # above A
    row1[node_centers[3]] = '|'   # above B

    function place!(buf, s, center)
        start = center - div(length(s) - 1, 2)
        for (k, ch) in enumerate(s)
            idx = start + k - 1
            1 <= idx <= length(buf) && (buf[idx] = ch)
        end
    end

    for (l, m) in zip(lbls, mids)
        place!(row3, l, m)
    end

    println(io, "SASBCell:")
    println(io, rstrip(join(row1)))
    println(io, row2)
    print(io,   rstrip(join(row3)))
end

# compact show (e.g. inside containers, print())
function Base.show(io::IO, c::SASBCell)
    _SASBCell_summary(io, c)
end

# rich display: REPL, notebooks, text/plain
function Base.show(io::IO, ::MIME"text/plain", c::SASBCell)
    _SASBCell_summary(io, c)
end
