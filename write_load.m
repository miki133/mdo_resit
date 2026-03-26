function write_load(x_spanwise, lift_distribution, moment_distribution)

    fid = fopen( 'B767_EMWET.load','wt');

    for i = 1:length(x_spanwise)
        fprintf(fid, '%g %g %g\n',x_spanwise(i), lift_distribution(i), moment_distribution(i));
    end

    fclose(fid);

end