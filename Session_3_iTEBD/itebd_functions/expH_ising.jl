
function exphAB_ising(gz::Float64, dt::Number, s)
    A = s[1]
    B = s[2]
    gz = gz / 2
    hAB = op("X", A) * op("X", B) +
          gz * op("Z", A) * op("Id", B) + gz * op("Id", A) * op("Z", B)
    return exp(im * dt * hAB)
end

"""
    exphAB_ising(p::IsingParams, dt, s)

Two-site gate ``\\exp(i\\,dt\\,h_{AB})`` for the Ising Hamiltonian
``H = J\\sum XX - g_\\perp \\sum Z - h_\\parallel \\sum X``,
where the on-site fields are split evenly across the AB and BA bonds.
"""
function exphAB_ising(p::IsingParams, dt::Number, s)
    A = s[1]
    B = s[2]
    g2 = p.gperp / 2
    h2 = p.hpar  / 2
    hAB = p.Jtwo * op("X", A) * op("X", B) +
          g2 * op("Z", A) * op("Id", B) + g2 * op("Id", A) * op("Z", B) +
          h2 * op("X", A) * op("Id", B) + h2 * op("Id", A) * op("X", B)
    return exp(im * dt * hAB)
end


function exphAB_ising_rotated(gz::Float64, dt::Number, s)
    gate = exphAB_ising(gz, dt, s)
    t1 = sim(s[2]', tags="time")
    t2 = sim(s[2], tags="time")
    return replaceinds(gate, (s[1], s[2], s[1]', s[2]') => (t1', t1, t2', t2))
end

function exphAB_ising_rotated(p::IsingParams, dt::Number, s)
    gate = exphAB_ising(p, dt, s)
    t1 = sim(s[2]', tags="time")
    t2 = sim(s[2], tags="time")
    return replaceinds(gate, (s[1], s[2], s[1]', s[2]') => (t1', t1, t2', t2))
end


"""
    trotter2_gates_ising(gz, dt, s)       -> (gate_AB_half, gate_BA_full)
    trotter2_gates_ising(p::IsingParams, dt, s) -> (gate_AB_half, gate_BA_full)

Return the gate pair for one second-order (Strang) Trotter step:

    gate_AB_half = exp(im * dt/2 * h_AB)
    gate_BA_full = exp(im * dt   * h_BA)

Apply as:
    c = step_AB(c, gate_AB_half)
    c = step_BA(c, gate_BA_full)
    c = step_AB(c, gate_AB_half)

This gives O(dt³) local Trotter error vs O(dt²) for the first-order scheme.
"""
function trotter2_gates_ising(gz::Float64, dt::Number, s)
    gate_AB_half = exphAB_ising(gz, dt / 2, s)
    gate_BA_full = exphAB_ising(gz, dt,     reverse(s))
    return (gate_AB_half, gate_BA_full)
end

function trotter2_gates_ising(p::IsingParams, dt::Number, s)
    gate_AB_half = exphAB_ising(p, dt / 2, s)
    gate_BA_full = exphAB_ising(p, dt,     reverse(s))
    return (gate_AB_half, gate_BA_full)
end


# ── X | Z Trotter splitting ───────────────────────────────────────────────────
#
# H = H_X + H_Z  where
#   H_X = J ∑ XᵢXᵢ₊₁ + h∥ ∑ Xᵢ   (all terms commute → exact layer)
#   H_Z = g⊥ ∑ Zᵢ                 (all terms commute → exact layer)
#
# exp(dt·H_X) and exp(dt·H_Z) are each exact products of local gates;
# the only Trotter error comes from [H_X, H_Z] ≠ 0.
#
# 1st-order:  exp(dt·H_X) · exp(dt·H_Z)
# 2nd-order (ZXZ):  exp(dt/2·H_Z) · exp(dt·H_X) · exp(dt/2·H_Z)

"""
    expHX_bond_ising(p::IsingParams, dt, s)

Two-site gate ``\\exp(i\\,dt\\,h_X)`` for the X-part of the Ising Hamiltonian
on bond `(s[1], s[2])`:

    h_X = J · X_{s[1]} X_{s[2]} + (h∥/2) · (X_{s[1]} + X_{s[2]})

The h∥/2 factor splits the on-site X field evenly between the AB and BA bonds.
All terms in H_X commute, so the full layer ``\\prod_\\text{bonds} \\exp(dt\\,h_X)``
is exact (no intra-layer Trotter error).
"""
function expHX_bond_ising(p::IsingParams, dt::Number, s)
    A, B = s[1], s[2]
    h2 = p.hpar / 2
    HX = p.Jtwo * op("X", A) * op("X", B) +
         h2 * op("X", A) * op("Id", B) +
         h2 * op("Id", A) * op("X", B)
    return exp(im * dt * HX)
end


"""
    expHZ_site_ising(p::IsingParams, dt, i::Index)

Single-site gate ``\\exp(i\\,dt\\,g_\\perp Z_i)`` for the Z-part of the Ising
Hamiltonian on site `i`. Apply independently to every site in the unit cell.
All Z gates commute, so the full layer is exact.
"""
function expHZ_site_ising(p::IsingParams, dt::Number, i::Index)
    return exp(im * dt * p.gperp * op("Z", i))
end


"""
    trotter2_XZ_gates_ising(p::IsingParams, dt, s)
        -> (gate_X_full, gate_ZA_half, gate_ZB_half)

Return the three gates for one second-order ZXZ Trotter step:

    exp(dt/2 · H_Z) · exp(dt · H_X) · exp(dt/2 · H_Z)

where:
  - `gate_X_full`  = `expHX_bond_ising(p, dt, s)`      — full-step 2-site X gate
  - `gate_ZA_half` = `expHZ_site_ising(p, dt/2, s[1])` — half-step Z gate for site A
  - `gate_ZB_half` = `expHZ_site_ising(p, dt/2, s[2])` — half-step Z gate for site B

Apply as (per step):
    apply gate_ZA_half, gate_ZB_half to each site
    apply gate_X_full  to each bond
    apply gate_ZA_half, gate_ZB_half to each site
"""
function trotter2_XZ_gates_ising(p::IsingParams, dt::Number, s)
    gate_X_full  = expHX_bond_ising(p, dt,     s)
    gate_ZA_half = expHZ_site_ising(p, dt / 2, s[1])
    gate_ZB_half = expHZ_site_ising(p, dt / 2, s[2])
    return (gate_X_full, gate_ZA_half, gate_ZB_half)
end
