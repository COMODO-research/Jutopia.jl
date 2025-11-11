function optimality_criteria(x, volfrac, dc, dv, element_volumes, filter_type, neighbors, weights, n)

    # optimality_criteria(x, volfrac, dc, dv, element_volumes, filter_type, neighbors, weights, n)
    l1 = 0.0
    l2 = 1e9
    move = 0.2
    total_volume = sum(element_volumes)
    xnew = copy(x)
    xPhys = copy(x)  # ✅ Ensure xPhys is always defined

    while (l2 - l1 > 1e-3)
        lmid = 0.5 * (l2 + l1)

        # OC update rule
        xnew = @. max(0.0, max(x - move, min(1.0, min(x + move, x * sqrt(-dc / (lmid * dv))))))
        # --- Filtering logic ---
        if filter_type == :sensitivity_filter
            xPhys = xnew
        elseif filter_type == :density_filter
            #xPhys = rho_filter(xnew, neighbors, weights, n)
            xPhys = density_filter(n, xnew , element_volumes, neighbors, weights)
        else
            error("the filter type $filter_type is incorrect")
        end

        # --- Volume constraint check ---
        if sum(xPhys .* element_volumes) - volfrac * total_volume > 0
            l1 = lmid
        else
            l2 = lmid
        end
    end

    return xnew, xPhys
end
