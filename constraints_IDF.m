function [c,ceq] = constraints_IDF(x)

global couplings;

% Constants:

W_fuel_ref = 1;
W_wing_ref = 1;
L_D_ref = 1;
wing_loading_ref = 552.38 * 9.81;
rho_fuel = 0.817*10^3;
f_tank = 0.93;
v_fuel_ref = 1;

% Guesses:

W_fuel_hat = x(19);
W_wing_hat = x(20);
L_D_hat = x(21);

% Disciplines from global couplings:

W_fuel = couplings.W_fuel;
W_wing = couplings.W_wing;
L_D_ratio = couplings.L_D_ratio;
wing_loading = couplings.wing_loading;
V_tank = couplings.V_tank;

% Wing loading constraint:

c_wing_loading = (wing_loading - wing_loading_ref)/wing_loading_ref;

% Fuel Volume constraint:

v_fuel = W_fuel / f_tank / rho_fuel;
c_fuel_volume = (v_fuel - V_tank)/v_fuel_ref;

% Emissions constraint:

c_emissions = (W_fuel - W_fuel_ref) / W_fuel_ref;

% Consistency constraints

c_w_fuel = abs(W_fuel - W_fuel_hat) / W_fuel_ref;
c_w_wing = abs(W_wing - W_wing_hat) / W_wing_ref;
c_L_D = abs(L_D_ratio - L_D_hat) / L_D_ref;

%Feed back into optimizer

c = [c_wing_loading, c_fuel_volume, c_emissions];
ceq = [c_w_fuel, c_w_wing, c_L_D];

end