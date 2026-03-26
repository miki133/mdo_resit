function [f,vararg] = Optim_IDFGauss_hybrid(x)
    
    mach = x(1) * 0.80199380281831;
    h = x(2) * 1.18872e+04;
    c_root = x(3) * 11.3981;
    outboard_span = x(4) * 15.8602;
    outboard_taper_ratio = x(5) * 0.360667590375573;
    sweep_LE = x(6) * 31.5;
    a_upper = x(7:12) .* [0.1823, 0.1195, 0.1490, 0.2179, 0.1249, 0.2745];   % 6×1 vector
    a_lower = x(13:18) .* [-0.1809, -0.1261, -0.1357, -0.2602, -0.0796, 0.1409];  % 6×1 vector
    W_fuel = x(19);
    W_wing = x(20);
    W_aw = 1;

    [~, a, ~, rho] = atmosisa(h);
    v = a * mach;
    dynamicPressure = 0.5 * rho * v^2;
    drag_wingless = 2.774713572457221e+4;
    dynamic_pressure_reference = 1.001172405760481e+4;

    span_inboard = 7.9248;
    c_kink = c_root - span_inboard * cotd(90 - sweep_LE);
    c_tip = c_kink * outboard_taper_ratio;

    surface_inboard = span_inboard * (c_root + c_kink)/2;
    surface_outboard = outboard_span * (c_kink + c_tip) /2;

    wing_surface = 2 * (surface_outboard + surface_inboard);

    update_airfoil(a_upper,a_lower)
    
    % Loads discipline
    [spanwise_positions, lift_distribution, moment_distribution] = loads(x);
    % Structures discipline
    wing_weight = structures(x, spanwise_positions, lift_distribution, moment_distribution);
    
    % Aerodynamics discipline
    CLdes = (sqrt((W_wing + W_fuel + W_aw) * (W_wing + W_aw))) / (0.5 * rho * v^2 * wing_surface);
    aero = aerodynamics(a_upper, a_lower, CLdes, mach, h, c_root, outboard_span, outboard_taper_ratio, sweep_LE, 1);
    spanwise_positions = aero.Wing.Yst;
    chord_distribution = aero.Wing.chord;
    L_D_ratio = (aero.CLwing * dynamicPressure * wing_surface) / (aero.CDwing * dynamicPressure * wing_surface + drag_wingless * dynamicPressure/dynamic_pressure_reference);
    
    % Performance discipline
    fuel_weight = performance(x);
    
    % Objective
    f = objective_mtow(W_fuel, W_wing);

    global couplings;
    
    couplings.W_fuel = fuel_weight;
    couplings.W_wing = wing_weight;
    couplings.L_D_ratio = L_D_ratio;
    couplings.wing_loading = (wing_weight + fuel_weight + W_aw) / wing_surface;
    couplings.V_tank = find_tank_volume(a_upper, a_lower, spanwise_positions, chord_distribution, outboard_span);
    
    vararg = {CL_design,mtow,counter,mtow_c, wing_weight, couplings.V_tank, couplings.Re, aero, L_D_ratio};

end


