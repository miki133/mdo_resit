%% Main Optimization file for B767-300 wing
clear all
close all
clc

global iterHist fvalHist conHist
fvalHist = [];   % Objective value per iteration
iterHist = [];   % Iteration counter
conHist  = [];   % rows = iterations, columns = constraints

%Initial values and constants:

span = 47.57 / 2;
inboard_span = 7.9248;
sweep_LE = 31.5;
reference_range = 7445040; % in meters

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

sweep_lb = 24.49 / sweep_LE;
sweep_ub = 41.87 / sweep_LE;

outboard_span_lb = 10.08 / outboard_span;
outboard_span_ub = 18.08 / outboard_span;

mach_h_lb = 0.9;
mach_h_ub = 1.1;

c_root_lb = 0.9;
c_root_ub = 1.1;

outboard_taper_ratio_lb = 0.25 / outboard_taper_ratio;
outboard_taper_ratio_ub = 0.4 / outboard_taper_ratio;

%bounds
lb = [ ...
    mach_h_lb, mach_h_lb , c_root_lb, outboard_span_lb, outboard_taper_ratio_lb, sweep_lb, ...
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];

ub = [ ...
    mach_h_ub, mach_h_ub, c_root_ub, outboard_span_ub, outboard_taper_ratio_ub, sweep_ub, ...
    1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5];

x0 = [mach_ref, h_ref, c_root, outboard_span, outboard_taper_ratio, sweep_LE, a_upper_1, a_upper_2, a_upper_3, a_upper_4, a_upper_5, a_upper_6, a_lower_1, a_lower_2, a_lower_3, a_lower_4, a_lower_5, a_lower_6]; % initial design variables vector
normal_vector = x0;
x0 = x0 ./ normal_vector; % normalize first 

% Options for the optimization
options.Display         = 'iter-detailed';
options.Algorithm       = 'sqp';
options.FunValCheck     = 'off';
options.DiffMinChange   = 1e-2;         % Minimum change while gradient searching
options.DiffMaxChange   = 1e-1;         % Maximum change while gradient searching
options.TolCon          = 1e-3;         % Maximum difference between two subsequent constraint vectors [c and ceq]
options.TolFun          = 1e-3;         % Maximum difference between two subsequent objective value
options.TolX            = 1e-4;         % Maximum difference between two subsequent design vectors
options.MaxIter         = 50;           % Maximum iterations
options.OutputFcn       = @outfun;

tic;
[x,FVAL,EXITFLAG,OUTPUT] = fmincon(@(x) Optim_IDFGauss_hybrid(x),x0,[],[],[],[],lb,ub,@(x) constraints(x),options);
toc;

[f,vararg] = Optim_IDFGauss_hybrid(x);
[CL_design, mtow, counter, mtow_c, wing_weight, V_tank, Re, aero, L_D_ratio] = vararg{:};

f = -f * reference_range;
x_final = x .* normal_vector;
mach = x_final(1);
h = x_final(2);
c_root = x_final(3);
outboard_span = x_final(4);
outboard_taper_ratio = x_final(5);
sweep_LE = x_final(6);
a_upper = x_final(7:12);
a_lower = x_final(13:18);

[~, a, ~, rho, nu] = atmosisa(h);
v = mach * a;
vc_ref = 236.644; % m/s
h_ref = 11887.2; % m
ct_ref = 1.8639e-4;
eta = exp( -((v - vc_ref)^2) / (2*70^2) -((h - h_ref)^2) / (2*2500^2) );
ct = ct_ref ./ eta;

aero_CDI = aerodynamics(a_upper, a_lower, CL_design, mach, h, c_root, outboard_span, outboard_taper_ratio, sweep_LE, 0);
CDi = aero_CDI.CDiwing;

span_inboard = 7.9248;
c_kink = c_root - span_inboard * cotd(90 - sweep_LE);
c_tip = c_kink * outboard_taper_ratio;
lambda1 = c_kink/c_root;
lambda2 = c_tip/c_kink;
surface_inboard = span_inboard * (c_root + c_kink)/2;
surface_outboard = outboard_span * (c_kink + c_tip) /2;
wing_surface = 2 *(surface_outboard + surface_inboard);
mac1 = (2/3) * c_root * ((1 + lambda1 + lambda1^2) / (1+lambda1));
mac2 = (2/3) * c_kink * ((1 + lambda2 + lambda2^2) / (1+lambda2));
mac_overall = (mac1*surface_inboard + mac2*surface_outboard)/(surface_outboard + surface_inboard);
reynolds_ref = mach_ref * a_ref * mac_overall / nu_ref;

alpha = aero.Alpha;

drag_wingless_ref = 2.774713572457221e+4;
dynamic_pressure_reference = 1.001172405760481e+4;
dynamic_pressure = 0.5 * rho * v^2;
drag_wingless = drag_wingless_ref * dynamic_pressure / dynamic_pressure_reference;
cd_wingless = drag_wingless / wing_surface / dynamic_pressure;
cd_wingless_ref = drag_wingless_ref / 283.3 / dynamic_pressure_reference;

weight_wingless_plusfuel = mtow - wing_weight;
wing_loading = mtow / wing_surface;
AR = ( 2*span_inboard + 2*outboard_span)^ 2 / wing_surface;

figure;
plot(iterHist, fvalHist, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Iteration');
ylabel('Objective Function Value');
title('fmincon Convergence History');

tolCon = options.TolCon;

figure; hold on; grid on;

nCon = size(conHist,2);
colors = lines(nCon);

for i = 1:nCon
    plot(iterHist, conHist(:,i), ...
        'LineWidth',1.3, ...
        'Color',colors(i,:));
end

yline( tolCon,'k--','TolCon');
yline(-tolCon,'k--','HandleVisibility','off');

xlabel('Iteration');
ylabel('Constraint value');
title('Constraint Convergence History');

legend(arrayfun(@(i) sprintf('Constraint %d',i), ...
       1:nCon,'UniformOutput',false), ...
       'Location','best');

function stop = outfun(x, optimValues, state)

    global iterHist fvalHist conHist
    persistent nCon

    stop = false;
    reference_range = 7445040;

    switch state
        case 'init'
            iterHist = [];
            fvalHist = [];
            conHist  = [];
            nCon     = [];

        case 'iter'
            iterHist(end+1,1) = optimValues.iteration;
            fvalHist(end+1,1) = -optimValues.fval * reference_range;

            [c, ~] = constraints(x);

            if isempty(nCon)
                nCon = numel(c);
                conHist = zeros(0,nCon);
            end

            conHist(end+1,:) = c(:).';

        case 'done'
            % nothing needed
    end
end
