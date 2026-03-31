function[V_tank] = find_tank_volume(a_upper, a_lower, spanwise_positions, chord_distribution, outboard_span)


    X_vect = linspace(0,1,103)'; 
    [upper,lower] = D_airfoil2(a_upper,a_lower,X_vect);
    
    xu = upper(:,1);
    yu = upper(:,2);
    
    xl = lower(:,1);
    yl = lower(:,2);
    
    % Spar locations
    x1 = 0.2;
    x2 = 0.6;
    
    x = linspace(x1, x2, 1000);
    
    yu_i = interp1(xu, yu, x, 'pchip');
    yl_i = interp1(xl, yl, x, 'pchip');
    
    thickness = yu_i - yl_i;
    area_between_spars_normalized = trapz(x, thickness);
    span_inboard = 7.9248;
    span = span_inboard + outboard_span;
    last_pos = 0.85 * span; % Tanks end at 85% of span

    % Refined spanwise grid
    y_fine = linspace(0, last_pos, 500);

    % Interpolate chord distribution
    chord_fine = interp1(spanwise_positions, chord_distribution, y_fine, 'pchip');

    % Compute volume integral
    volume_integrand = area_between_spars_normalized .* chord_fine.^2;

    V_tank = trapz(y_fine, volume_integrand);

end