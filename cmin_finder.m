function [cmin, xi_min] = cmin_finder(x, h, c, hmin)
    
    xi_min = fsolve(@(y) interp1(x,h, y) - hmin, zeros(size(hmin))); 
    cmin= interp1(x,c,xi_min); 

end

