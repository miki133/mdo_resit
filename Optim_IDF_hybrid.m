function [f,vararg] = Optim_IDF_hybrid(x)
    
    W_fuel = x(19);
    W_wing = x(20);
    w_fuel_ref = 29274.3142414032 * 9.81;
    w_wing_ref = 11445 * 9.81;
    % Objective
    f = - objective_mtow(W_fuel * w_fuel_ref, W_wing * w_wing_ref);
    
    vararg = {};

end


