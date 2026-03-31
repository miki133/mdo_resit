wfuel_ref = 31477.7845000000 * 9.81;
vfuel_init = 54.53977968; % in cubic meters, obtained from reference fuel weight and divided by the fuel density from assignment
h = 11887.2;
[~, a, ~, rho] = atmosisa(h);
Mach = 236.644 / a;
vc_ref = 236.644; % m/s
v_MO = 251.563;
ct_ref = 1.8639e-4;
W_fuel = 29274.3142414032 * 9.81;
W_TO = 156489 * 9.81;
v = Mach * a;  % m/s
taper = 0.207;
visc = 0;
span = 47.57 / 2;
inboard_span = 7.9248;
sweep_LE = 34.6;
reference_range = 3260 * 1000; % in meters

c_root = 11.3981;
c_tip = c_root * 0.207;
c_kink = c_root - inboard_span * cotd(90 - sweep_LE);

%Performance
h_ref = 11887.2; % m
[~, a_ref, ~, rho_ref, nu_ref] = atmosisa(h_ref);
mach_ref = 236.644 / a_ref;

outboard_span = span - inboard_span;
outboard_taper_ratio = c_tip / c_kink;

% Upper surface coefficients
a_upper_1 = 0.1823;
a_upper_2 = 0.1195;
a_upper_3 = 0.1490;
a_upper_4 = 0.2179;
a_upper_5 = 0.1249;
a_upper_6 = 0.2745;

% Lower surface coefficients
a_lower_1 = -0.1809;
a_lower_2 = -0.1261;
a_lower_3 = -0.1357;
a_lower_4 = -0.2602;
a_lower_5 = -0.0796;
a_lower_6 = 0.1409;

w_fuel_ref = 31477.7845 * 9.81;
w_wing_ref = 11445 * 9.81;
l_d_ref = 16;
v_fuel_ref = 38.528466644823280;
mtow_ref = 156489 * 9.81;

x = [mach_ref, h_ref, c_root, outboard_span, outboard_taper_ratio, sweep_LE,...
    a_upper_1, a_upper_2, a_upper_3, a_upper_4, a_upper_5, a_upper_6, a_lower_1,...
    a_lower_2, a_lower_3, a_lower_4, a_lower_5, a_lower_6, w_fuel_ref, w_wing_ref, l_d_ref];
a_upper = x(7:12);
a_lower = x(13:18);

span = 47.57 / 2;
span_inboard = 7.9248;
c_kink = c_root - span_inboard * cotd(90 - sweep_LE);
c_tip = c_root * taper;

outboard_span = span - span_inboard;
outboard_taper_ratio = c_root * taper / c_kink;
surface_inboard = span_inboard * (c_root + c_kink)/2;
surface_outboard = outboard_span * (c_kink + c_tip) /2;

lambda1 = c_kink/c_root;
lambda2 = c_tip/c_kink;
mac1 = (2/3) * c_root * ((1 + lambda1 + lambda1^2) / (1+lambda1));
mac2 = (2/3) * c_kink * ((1 + lambda2 + lambda2^2) / (1+lambda2));
mac_overall = (mac1*surface_inboard + mac2*surface_outboard)/(surface_outboard + surface_inboard);

wing_surface = 2 * (surface_outboard + surface_inboard);

CL = 2.5 * W_TO / (0.5 * rho * v_MO^2 * wing_surface);

aero = aerodynamics(CL, x, visc);
x_spanwise = aero.Wing.Yst ./ span;
load_factor = 1;

lift_distribution = load_factor .* aero.Wing.ccl .* 0.5 .* rho .* v_MO.^2;
moment_distribution = load_factor .* aero.Wing.cm_c4 .* 0.5 .* rho .* v_MO.^2 .* aero.Wing.chord .* mac_overall;

% Ensure column vectors
x_spanwise = x_spanwise(:);
lift_distribution = lift_distribution(:);
moment_distribution = moment_distribution(:);

% Interpolate lift values at x = 0 and x = 1
lift_0 = interp1(x_spanwise, lift_distribution, 0, 'pchip', 'extrap');
lift_1 = interp1(x_spanwise, lift_distribution, 1, 'pchip', 'extrap');

mom0 = interp1(x_spanwise, moment_distribution, 0, 'pchip', 'extrap');
mom1 = interp1(x_spanwise, moment_distribution, 1, 'pchip', 'extrap');

% Add them to the vectors
x_spanwise = [0; x_spanwise; 1];
lift_distribution = [lift_0; lift_distribution; lift_1];
moment_distribution = [mom0; moment_distribution; mom1];

write_init(156489, 156489-29274.3142414032, c_root, outboard_span, outboard_taper_ratio, sweep_LE)
write_load(x_spanwise, lift_distribution, moment_distribution)
EMWET B767_EMWET

fid = fopen('B767_EMWET.weight','r');
line = fgetl(fid);
    fclose(fid);

wing_weight = sscanf(line, 'Wing total weight(kg) %f');
disp(wing_weight)

eta = exp( -((v - vc_ref)^2) / (2*70^2) -((h - h_ref)^2) / (2*2500^2) );

ct = ct_ref ./ eta;

W_end_start_ratio = (1 - W_fuel/W_TO) / 0.938;

range_ref = (v / ct_ref) * l_d_ref * log(1 / W_end_start_ratio);

V_tank_ref = find_tank_volume(a_upper, a_lower, aero.Wing.Yst, aero.Wing.chord, outboard_span);

W_aw = W_TO - wing_weight*9.81 - wfuel_ref;

CL = sqrt(W_TO*(W_TO-wfuel_ref)) / (0.5 * rho * v^2 * wing_surface);
aero = aerodynamics(CL, x, 1);

cd_aw = aero.CLwing/16 - aero.CDwing;
drag_aw = cd_aw * 0.5 * rho * vc_ref^2 * wing_surface;
%w_fuel_ref = 38.5285 * 0.817 * 1000;

