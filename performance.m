function[fuel_weight] = performance(x)

    h = x(2);
    mach = x(1);
    wing_weight = x(20);
    L_D_hat = x(21);
    [~, a] = atmosisa(h);
    v = a * mach;
    weight_wingless = 115769.685758597 * 9.81;
    
    vc_ref = 236.644; % m/s
    h_ref = 11887.2; % m
    ct_ref = 1.8639e-4;
    range_ref = 2907 * 1000;

    eta = exp( -((v - vc_ref)^2) / (2*70^2) -((h - h_ref)^2) / (2*2500^2));

    ct = ct_ref ./ eta;
    
    W_start_end_ratio = exp(range_ref * ct / v / L_D_hat);
    c = 1 - 0.938 * (1/W_start_end_ratio);

    fuel_weight = c * (wing_weight + weight_wingless) / (1 - c);

end