function[mtow] = objective_mtow(W_fuel, W_wing)

mtow_ref = 156489 * 9.81;
W_aw = 1114084.574055;
mtow = (W_fuel + W_wing + W_aw) / mtow_ref;

end