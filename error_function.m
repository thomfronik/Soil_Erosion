function [e, index] = error_function(cmin, cplus, xi_min, xi_plus)
    err = zeros(size(cmin)); 
    for i = 1:length(xi_plus)
        x_plus = xi_plus(i);
        x_min = xi_min(i);
        c_min = cmin(i);
        c_plus = cplus(i); 
        intt = (1/(x_plus-x_min))*trapz(xi_min(1:i), cmin(1:i),1) + trapz(xi_plus(1:i), cplus(1:i),1); 
        intt 
        err(i) = (c_plus - c_min)^2 + (1- intt).^2;  
    end 
    [e,index] = min(err); 
    end
    
    