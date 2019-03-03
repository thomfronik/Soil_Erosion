
function [err_c_condition, err_int_condition] = error_function(x,c, xmin, xplus, int_points) 

    err_c_condition = abs(interp1(x,c,xmin) - interp1(x,c,xplus));
    err_int_condition = (1-(1/abs(xmin - xplus))*trapz(interp1(x,c, linspace(xmin, xplus, int_points))))^2; 
    err = 10*err_c_condition + err_int_condition;

end

