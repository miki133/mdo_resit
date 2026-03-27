wfuel_ref = 44559 * 9.81;
l_d_ref = 16;
vfuel_ref = 54.53977968; % in cubic meters, obtained from reference fuel weight and divided by the fuel density from assignment
h = 11887.2;
[~, a, ~, rho] = atmosisa(h);
Mach = 236.644 / a;
vc_ref = 236.644; % m/s
v_MO = 251.563;
h_ref = 11887.2; % m
ct_ref = 1.8639e-4;
W_fuel = 44559 * 9.81;
W_TO = 156489 * 9.81;
v = Mach * a;  % m/s
c_root = 11.3981;
sweep_LE = 31.5;
taper = 0.207;
visc = 0;
a_upper = [    0.1823
    0.1195
    0.1490
    0.2179
    0.1249
    0.2745]';
a_lower = [   -0.1809
   -0.1261
   -0.1357
   -0.2602
   -0.0796
    0.1409]';

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

CL = 2.5 * W_TO / (0.5 * rho * v^2 * wing_surface);

aero = aerodynamics(a_upper, a_lower, CL, Mach, h, c_root, outboard_span, outboard_taper_ratio, sweep_LE, visc);
x_spanwise = aero.Wing.Yst;
load_factor = 1;
lift_distribution = load_factor .* aero.Wing.ccl .* 0.5 .* rho .* v_MO.^2;
moment_distribution = load_factor .* aero.Wing.cm_c4 .* 0.5 .* rho .* v_MO.^2 .* aero.Wing.chord .* mac_overall;

write_init(156489, 156489-44559, c_root, outboard_span, outboard_taper_ratio, sweep_LE)
write_load(x_spanwise, lift_distribution, moment_distribution)
EMWET B767_EMWET

fid = fopen('B767_EMWET.weight','r');
line = fgetl(fid);
    fclose(fid);

wing_weight = sscanf(line, 'Wing total weight(kg) %f');


eta = exp( -((v - vc_ref)^2) / (2*70^2) -((h - h_ref)^2) / (2*2500^2) );

ct = ct_ref ./ eta;

W_end_start_ratio = (1 - W_fuel/W_TO) / 0.938;

range_ref = (v / ct) * l_d_ref * log(1 / W_end_start_ratio);

V_tank_ref = find_tank_volume(a_upper, a_lower, x_spanwise, aero.Wing.chord, outboard_span);

W_aw = W_TO - wing_weight*9.81 - wfuel_ref;

%cd_aw = aero.CLwing/16 - aero.CDwing;
%drag_aw = cd_aw * 0.5 * rho * vc_ref^2 * wing_surface;
%total_lift = 2 * trapz(x_spanwise, lift_distribution);