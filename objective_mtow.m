function[mtow] = objective_mtow(W_fuel, W_wing)

mtow_ref = 156489 * 9.81;
W_aw = 115769.685758597 * 9.81;
mtow = (W_fuel + W_wing + W_aw) / mtow_ref;

end