function omega_vec = calculate_aging_coefficients(E_cap_agg, C_bat, eta_dis, N_seg)
    A = 330.33; z = 0.552;
    Ea = 31500; R_gas = 8.314; Temp = 298;
    Arrhenius_Factor = A * exp(-Ea / (R_gas * Temp)); 
    E_rate_MWh = E_cap_agg / 1000; 
    omega_vec = zeros(N_seg, 1);
    for m = 1:N_seg
        h_upper = m / N_seg;
        h_lower = (m - 1) / N_seg;
        Q_upper = Arrhenius_Factor * (h_upper * E_rate_MWh)^z;
        Q_lower = Arrhenius_Factor * (h_lower * E_rate_MWh)^z;
        diff_Q = Q_upper - Q_lower;
        omega_vec(m) = (C_bat / (2 * eta_dis * E_rate_MWh)) * diff_Q;
    end
end