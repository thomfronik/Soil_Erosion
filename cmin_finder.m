function [cmin, xi_min] = cmin_finder(h_spline, c_spline, hmin)
    
    xi_min = fsolve(@(x) ppval(h_spline,x) - hmin, zeros(size(hmin))); 
    cmin= ppval(c_spline, xi_min); 

end

