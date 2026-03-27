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

w_fuel_ref = 31477.7845 * 9.81;
w_wing_ref = 11445 * 9.81;
l_d_ref = 16;
v_fuel_ref = 38.528466644823280;
mtow_ref = 156489 * 9.81;

%bounds
lb = [ ...
    mach_h_lb, mach_h_lb , c_root_lb, outboard_span_lb, outboard_taper_ratio_lb, sweep_lb, ...
    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5];

ub = [ ...
    mach_h_ub, mach_h_ub, c_root_ub, outboard_span_ub, outboard_taper_ratio_ub, sweep_ub, ...
    1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5];
 % initial design variables vector
x0 = [mach_ref, h_ref, c_root, outboard_span, outboard_taper_ratio, sweep_LE,...
    a_upper_1, a_upper_2, a_upper_3, a_upper_4, a_upper_5, a_upper_6, a_lower_1,...
    a_lower_2, a_lower_3, a_lower_4, a_lower_5, a_lower_6, w_fuel_ref, w_wing_ref, l_d_ref];

normal_vector = x0;
x0 = x0 ./ normal_vector; % normalize first 

% Options for the optimization
options.Display         = 'iter-detailed';
options.Algorithm       = 'sqp';
options.FunValCheck     = 'off';
options.DiffMinChange   = 1e-6;         % Minimum change while gradient searching
options.DiffMaxChange   = 5e-2;         % Maximum change while gradient searching
options.TolCon          = 1e-6;         % Maximum difference between two subsequent constraint vectors [c and ceq]
options.TolFun          = 1e-6;         % Maximum difference between two subsequent objective value
options.TolX            = 1e-6;         % Maximum difference between two subsequent design vectors
options.MaxIter         = 50;           % Maximum iterations
options.OutputFcn       = @outfun;

tic;
[x,FVAL,EXITFLAG,OUTPUT] = fmincon(@(x) Optim_IDFGauss_hybrid(x),x0,[],[],[],[],lb,ub,@(x) constraints_IDF(x),options);
toc;

[f,vararg] = Optim_IDFGauss_hybrid(x);
[wing_loading, V_tank, CLdes, Re, dynamicPressure] = vararg{:};

mtow_final = -f * mtow_ref;
x_final = x .* normal_vector;
mach = x_final(1);
h = x_final(2);
c_root = x_final(3);
outboard_span = x_final(4);
outboard_taper_ratio = x_final(5);
sweep_LE = x_final(6);
a_upper = x_final(7:12);
a_lower = x_final(13:18);
w_fuel = x(19);
w_wing = x(20);
L_D_ratio = x(21);

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
    global iterHist fvalHist conHist ceqHist
    persistent nCon nCeq

    stop = false;
    reference_mtow = 156489 * 9.81;

    switch state
        case 'init'
            iterHist = [];
            fvalHist = [];
            conHist  = [];
            ceqHist  = [];
            nCon     = [];
            nCeq     = [];

        case 'iter'
            iterHist(end+1,1) = optimValues.iteration;
            fvalHist(end+1,1) = -optimValues.fval * reference_mtow;

            [c, ceq] = constraints_IDF(x);

            % ── Initialize sizes on first iteration ───────────────────────
            if isempty(nCon)
                nCon    = numel(c);
                nCeq    = numel(ceq);
                conHist = zeros(0, nCon);
                ceqHist = zeros(0, nCeq);
            end

            conHist(end+1,:)  = c(:).';
            ceqHist(end+1,:)  = ceq(:).';   % consistency constraint history

        case 'done'
            % nothing needed
    end
end