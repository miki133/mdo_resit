function[aero_results] = aerodynamics(a_upper, a_lower, CL, mach, h, c_root, outboard_span, outboard_taper_ratio, sweep_LE, visc)

    % Wing planform geometry 
    wing_dihedral = 6; %deg
    span_inboard = 7.9248;
    twist = 3;
    span = span_inboard + outboard_span;

    c_kink = c_root - span_inboard * cotd(90 - sweep_LE);
    z_kink = span_inboard*sind(wing_dihedral);
    y_kink = span_inboard * cosd(wing_dihedral);
    x_kink = span_inboard * tand(sweep_LE);

    c_tip = c_kink * outboard_taper_ratio;
    x_tip = span * tand(sweep_LE);
    y_tip = span * cosd(wing_dihedral);
    z_tip = span*sind(wing_dihedral);
    

    %                x    y     z   chord(m)    twist angle (deg) 
    AC.Wing.Geom = [0     0     0     c_root         twist; % root
                    x_kink*1.00001  y_kink   z_kink     c_kink      twist; % kink
                    x_tip  y_tip   z_tip     c_tip         twist]; % tip
    
    % Wing incidence angle (degree)
    AC.Wing.inc  = 4.25;   
                
                
    % Airfoil coefficients input matrix
    AC.Wing.Airfoils   = [a_upper a_lower;
                          a_upper a_lower;
                          a_upper a_lower];
                      
    AC.Wing.eta = [0;span_inboard/span;1];  % Spanwise location of the airfoil sections

    AC.Visc  = visc;    % 0 for inviscid and 1 for viscous analysis
    
    AC.Aero.MaxIterIndex = 150;    %Maximum number of Iteration for the
                                    %convergence of viscous calculation
    [~, a, ~, rho, nu] = atmosisa(h);

    surface_inboard = span_inboard * (c_root + c_kink)/2;
    surface_outboard = outboard_span * (c_kink + c_tip) /2;

    lambda1 = c_kink/c_root;
    lambda2 = c_tip/c_kink;
    mac1 = (2/3) * c_root * ((1 + lambda1 + lambda1^2) / (1+lambda1));
    mac2 = (2/3) * c_kink * ((1 + lambda2 + lambda2^2) / (1+lambda2));
    mac_overall = (mac1*surface_inboard + mac2*surface_outboard)/(surface_outboard + surface_inboard);

    % Flight Condition
    AC.Aero.V     = mach * a;            % flight speed (m/s)
    AC.Aero.rho   = rho;         % air density  (kg/m3)
    AC.Aero.alt   = h;             % flight altitude (m)
    AC.Aero.Re    = mach * a * mac_overall / nu;        % reynolds number (bqased on mean aerodynamic chord)
    AC.Aero.M     = mach;           % flight Mach number 
    AC.Aero.CL    = CL;          % lift coefficient - comment this line to run the code for given alpha%
    %AC.Aero.Alpha = 0;             % angle of attack -  comment this line to run the code for given cl 
    
    aero_results = Q3D_solver(AC);
    
    global couplings;
    couplings.Re = mach * a * mac_overall / nu;

end
