
function err = error_function(c_spline, xmin, xplus, int_points) 

    err_c_condition = (ppval(c_spline, xmin) - ppval(c_spline, xplus))^2;
    err_int_condition = (1-trapz(ppval(c_spline, linspace(xmin, xplus, int_points))))^2; 
    err = err_c_condition + err_int_condition;

end

