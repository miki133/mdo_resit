%%%_____Routine to write the input file for the EMWET procedure________% %%
function write_init(MTOW, MZF, c_root, outboard_span, outboard_taper_ratio, sweep_LE)

    wing_dihedral = 6; %deg
    span_inboard = 7.9248;
    span = span_inboard + outboard_span;
    
    c_kink = c_root - span_inboard * cotd(90 - sweep_LE);
    z_kink = span_inboard*sind(wing_dihedral);
    y_kink = span_inboard * cosd(wing_dihedral);
    x_kink = span_inboard * tand(sweep_LE);

    x_tip = span * tand(sweep_LE);
    y_tip = span * cosd(wing_dihedral);
    z_tip = span*sind(wing_dihedral);
    c_tip = c_kink * outboard_taper_ratio;

    nz_max      =    2.5;   
    spar_front  =    0.2;
    spar_rear   =    0.6;
    ftank_start =    0.1;
    ftank_end   =    0.85;
    eng_num     =    1;
    eng_ypos    =    span_inboard / span;
    eng_mass    =    7954;         %kg

    E_al        =    7.1E10;       %N/m2
    rho_al      =    2800;         %kg/m3
    Ft_al       =    295E6;        %N/m2
    Fc_al       =    295E6;        %N/m2
    pitch_rib   =    0.5;          %[m]
    eff_factor  =    0.96;         %Depend on the stringer type, Ztype in this case
    Airfoil     =    'updated_airfoil';


    section_num =    3;
    airfoil_num =    3;

    surface_inboard = span_inboard * (c_root + c_kink)/2;
    surface_outboard = outboard_span * (c_kink + c_tip) /2;

    wing_surf   =   2 * (surface_outboard + surface_inboard);
    
    fid = fopen( 'B767_EMWET.init','wt');
    fprintf(fid, '%g %g \n', MTOW, MZF);
    fprintf(fid, '%g \n',nz_max);
    
    fprintf(fid, '%g %g %g %g \n',wing_surf,span * 2,section_num,airfoil_num);
    
    fprintf(fid, '0 %s \n',Airfoil);
    fprintf(fid, '%g %s \n',eng_ypos, Airfoil);
    fprintf(fid, '1 %s \n',Airfoil); 

    fprintf(fid, '%g %g %g %g %g %g \n',c_root,0,0,0,spar_front,spar_rear);
    fprintf(fid, '%g %g %g %g %g %g \n',c_kink,x_kink,y_kink,z_kink,spar_front,spar_rear); % kink
    fprintf(fid, '%g %g %g %g %g %g \n',c_tip,x_tip,y_tip,z_tip,spar_front,spar_rear);
    
    fprintf(fid, '%g %g \n',ftank_start,ftank_end);
    
    fprintf(fid, '%g \n', eng_num);
    fprintf(fid, '%g  %g \n', eng_ypos,eng_mass);
    
    fprintf(fid, '%g %g %g %g \n',E_al,rho_al,Ft_al,Fc_al);
    fprintf(fid, '%g %g %g %g \n',E_al,rho_al,Ft_al,Fc_al);
    fprintf(fid, '%g %g %g %g \n',E_al,rho_al,Ft_al,Fc_al);
    fprintf(fid, '%g %g %g %g \n',E_al,rho_al,Ft_al,Fc_al);
    
    fprintf(fid,'%g %g \n',eff_factor,pitch_rib);
    fprintf(fid,'0 \n');
    
    fclose(fid);

end