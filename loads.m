function[spanwise_positions, lift_distribution, moment_distribution] = loads(x)

    c_root = x(3);
    outboard_span = x(4);
    outboard_taper_ratio = x(5);
    sweep_LE = x(6);
    W_fuel_hat = x(19);
    W_wing_hat = x(20);
    W_aw = 1114084.574055;
    v_MO = 251.563;
    [~, ~, ~, rho] = atmosisa(x(2));

    span_inboard = 7.9248;
    c_kink = c_root - span_inboard * cotd(90 - sweep_LE);
    c_tip = c_kink * outboard_taper_ratio;
    
    surface_inboard = span_inboard * (c_root + c_kink)/2;
    surface_outboard = outboard_span * (c_kink + c_tip) /2;
    
    lambda1 = c_kink/c_root;
    lambda2 = c_tip/c_kink;
    mac1 = (2/3) * c_root * ((1 + lambda1 + lambda1^2) / (1+lambda1));
    mac2 = (2/3) * c_kink * ((1 + lambda2 + lambda2^2) / (1+lambda2));
    mac_overall = (mac1*surface_inboard + mac2*surface_outboard)/(surface_outboard + surface_inboard);
    wing_surface = 2 * (surface_outboard + surface_inboard);

    mtow_hat = W_aw + W_wing_hat + W_fuel_hat;
    CL_design = 2.5 * mtow_hat / (0.5 * rho * v_MO^2 * wing_surface);

    [q3d_loads] = aerodynamics(CL_design, x, 0);
       
    spanwise_positions = q3d_loads.Wing.Yst;
    lift_distribution = q3d_loads.Wing.ccl .* 0.5 .* rho .* v_MO^2;
    moment_distribution = q3d_loads.Wing.cm_c4 .* 0.5 .* rho .* v_MO^2 .* q3d_loads.Wing.chord * mac_overall;
   
    % Ensure column vectors
    spanwise_positions = spanwise_positions(:);
    lift_distribution = lift_distribution(:);
    moment_distribution = moment_distribution(:);
    
    % Interpolate lift values at x = 0 and x = 1
    lift_0 = interp1(spanwise_positions, lift_distribution, 0, 'pchip', 'extrap');
    lift_1 = interp1(spanwise_positions, lift_distribution, 1, 'pchip', 'extrap');
    
    mom0 = interp1(spanwise_positions, moment_distribution, 0, 'pchip', 'extrap');
    mom1 = interp1(spanwise_positions, moment_distribution, 1, 'pchip', 'extrap');
    
    % Add them to the vectors
    spanwise_positions = [0; spanwise_positions; 1];
    lift_distribution = [lift_0; lift_distribution; lift_1];
    moment_distribution = [mom0; moment_distribution; mom1];

end