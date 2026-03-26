function[spanwise_positions, lift_distribution, moment_distribution] = loads(x)

    
    h = x(2);
    c_root = x(3);
    outboard_span = x(4);
    outboard_taper_ratio = x(5);
    sweep_LE = x(6);
    a_upper = x(7:12);
    a_lower = x(13:18);
    W_fuel_hat = x(19);
    W_wing_hat = x(20);
    W_aw = 1;

    mtow_hat = W_aw + W_wing_hat + W_fuel_hat;
    CL_design = sqrt(mtow_hat * (mtow_hat - W_fuel_hat));

    [~, a, ~, rho] = atmosisa(h);
    v_MO = 251.563;
    mach = v_MO / a;
    [q3d_loads] = aerodynamics(a_upper, a_lower, CL_design, mach, h, c_root, outboard_span, outboard_taper_ratio, sweep_LE, 0);
       
    spanwise_positions = q3d_loads.Wing.Yst;
    lift_distribution = q3d_loads.Wing.ccl .* 0.5 .* rho .* v_MO^2;
    moment_distribution = q3d_loads.Wing.cm_c4 .* 0.5 .* rho .* v_MO^2 .* q3d_loads.Wing.chord;

end