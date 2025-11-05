function compliance_and_sensitivity(model_type, ρ, dh, u, cell_values, penalty, E0, ν, Emin)
    dC_dρ = zeros(length(ρ))
    C = 0.0

    ke = zeros(getnbasefunctions(cell_values), getnbasefunctions(cell_values))

    for (cell_index, cell) in enumerate(CellIterator(dh))
        reinit!(cell_values, cell)
        fill!(ke, 0.0)

        ρ_local = ρ[cell_index]  # scalar density for this element

        # Element stiffness (for reference or use)
        assemble_cell!(model_type, ke, cell_values, ρ_local, 0.0, E0, ν, Emin)

        eldofs = celldofs(cell)
        ue = u[eldofs]

        # --- Element constitutive interpolation ---
        coeff_comp = Emin + (ρ_local ^ penalty) * (E0 - Emin)

        # --- Element sensitivity term ---
        coeff_dc = -penalty * (E0 - Emin) * (ρ_local ^ (penalty - 1))

        # --- Accumulate global compliance ---
        C += coeff_comp * (ue' * ke * ue)

        # --- Elemental sensitivity ---
        dC_dρ[cell_index] = coeff_dc * (ue' * ke * ue)
    end

    return C, dC_dρ
end
