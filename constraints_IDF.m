function [c,ceq] = constraints_IDF(x)
    
    mach = x(1) * 0.80199380281831;
    h = x(2) * 1.18872e+04;
    c_root = x(3) * 11.3981;
    outboard_span = x(4) * 15.8602;
    outboard_taper_ratio = x(5) * 0.360667590375573;
    sweep_LE = x(6) * 31.5;
    a_upper = x(7:12) .* [0.1823, 0.1195, 0.1490, 0.2179, 0.1249, 0.2745];   % 6×1 vector
    a_lower = x(13:18) .* [-0.1809, -0.1261, -0.1357, -0.2602, -0.0796, 0.1409];  % 6×1 vector
    W_fuel = x(19) * 29274.3142414032 * 9.81;
    W_wing = x(20) * 10360.2 * 9.81;
    L_D_ratio_hat = x(21) * 16;
    W_aw = 115769.685758597 * 9.81;

    x_new = [ ...
    mach, h, c_root, outboard_span, outboard_taper_ratio, sweep_LE, ...
    a_upper(1), a_upper(2), a_upper(3), a_upper(4), a_upper(5), a_upper(6), ...
    a_lower(1), a_lower(2), a_lower(3), a_lower(4), a_lower(5), a_lower(6), ...
    W_fuel, W_wing, L_D_ratio_hat];

    [~, a, ~, rho] = atmosisa(h);
    v = a * mach;
    dynamicPressure = 0.5 * rho * v^2;
    drag_wingless = 29443.3155785225;
    dynamic_pressure_reference = 8859.44131748761;
    
    span_inboard = 7.9248;
    c_kink = c_root - span_inboard * cotd(90 - sweep_LE);
    c_tip = c_kink * outboard_taper_ratio;
    
    surface_inboard = span_inboard * (c_root + c_kink)/2;
    surface_outboard = outboard_span * (c_kink + c_tip) /2;
    
    wing_surface = 2 * (surface_outboard + surface_inboard);
    
    update_airfoil(a_upper,a_lower)
    
    % Loads discipline
    [spanwise_positions, lift_distribution, moment_distribution] = loads(x_new);
    % Structures discipline
    wing_weight = 9.81 * structures(x_new, spanwise_positions, lift_distribution, moment_distribution);
    
    % Aerodynamics discipline
    CLdes = (sqrt((W_wing + W_fuel + W_aw) * (W_wing + W_aw))) / (0.5 * rho * v^2 * wing_surface);
    aero = aerodynamics(CLdes, x_new, 1);
    spanwise_positions = aero.Wing.Yst;
    chord_distribution = aero.Wing.chord;
    L_D_ratio = (aero.CLwing * dynamicPressure * wing_surface) / (aero.CDwing * dynamicPressure * wing_surface + drag_wingless * dynamicPressure/dynamic_pressure_reference);
    
    % Performance discipline
    fuel_weight = performance(x_new);
    
    % Constants:
    
    W_fuel_ref = 29274.3142414032 * 9.81;
    W_wing_ref = 10360.2 * 9.81;
    L_D_ref = 16;
    wing_loading_ref = 552.38 * 9.81;
    rho_fuel = 0.817*10^3;
    f_tank = 0.93;
    v_fuel_ref = 38.528466644823280;
    
    
    % Wing loading constraint:
    wing_loading = (W_wing + W_fuel + W_aw) / wing_surface;
    c_wing_loading = (wing_loading - wing_loading_ref)/wing_loading_ref;
    
    % Fuel Volume constraint:
    
    v_fuel = W_fuel / 9.81 / f_tank / rho_fuel;
    V_tank = find_tank_volume(a_upper, a_lower, spanwise_positions, chord_distribution, outboard_span);
    c_fuel_volume = (v_fuel - V_tank)/v_fuel_ref;
    
    % Emissions constraint:
    
    c_emissions = (W_fuel - W_fuel_ref) / W_fuel_ref;
    
    % Consistency constraints
    
    c_w_fuel = abs(fuel_weight - W_fuel) / W_fuel_ref;
    c_w_wing = abs(wing_weight - W_wing) / W_wing_ref;
    c_L_D = abs(L_D_ratio - L_D_ratio_hat) / L_D_ref;
    
    %Feed back into optimizer
    
    c = [c_wing_loading, c_fuel_volume, c_emissions];
    ceq = [c_w_fuel, c_w_wing, c_L_D];

        % Names for debugging
    c_names   = {'wing_loading','fuel_volume','emissions'};
    ceq_names = {'fuel_consistency','wing_weight_consistency','L_D_consistency'};
    
    % Check inequality constraints
    for i = 1:length(c)
        if ~isfinite(c(i))
            fprintf('NaN/Inf in inequality constraint: %s\n', c_names{i});
        end
    end
    
    % Check equality constraints
    for i = 1:length(ceq)
        if ~isfinite(ceq(i))
            fprintf('NaN/Inf in equality constraint: %s\n', ceq_names{i});
        end
    end

end