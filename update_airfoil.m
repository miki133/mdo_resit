function update_airfoil(a_upper,a_lower)
% Generates airfoil x y coordinates using Bernstein polynomials for
% optimization loop

    X_vect = linspace(0,1,103)';      %points for evaluation along x-axis
    [Xtu,Xtl] = D_airfoil2(a_upper,a_lower,X_vect);
    
    fid = fopen('updated_airfoil.dat','wt');

    for i = 1:length(X_vect)
        fprintf(fid, '%g %g\n', Xtu(i, 1), Xtu(i, 2));
    end

    for i = 1:length(X_vect)
        fprintf(fid, '%g %g\n', Xtl(i, 1), Xtl(i, 2));
    end

    fclose(fid);


end
