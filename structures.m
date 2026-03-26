function[wing_weight] = structures(x, spanwise_positions, lift_distribution, moment_distribution)


    c_root = x(3);
    outboard_span = x(4);
    outboard_taper_ratio = x(5);
    sweep_LE = x(6);
    
    W_fuel = x(19);
    W_wing = x(20);
    W_aw = 1;
    mtow = W_fuel + W_wing + W_aw;
    mzf = mtow - W_fuel;
    
  
    write_init(mtow, mzf, c_root, outboard_span, outboard_taper_ratio, sweep_LE)
    write_load(spanwise_positions, lift_distribution, moment_distribution)

    EMWET B767_EMWET

    fid = fopen('B767_EMWET.weight','r');
    line = fgetl(fid);
        fclose(fid);

    wing_weight = sscanf(line, 'Wing total weight(kg) %f');

end