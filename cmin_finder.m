function cmin = cmin_finder(h_spline, c_spline, hplus, options)
xi_min = fsolve(@(x) ppval(h_spline,x) - hplus, zeros(1,length(hplus)), options); 
cmin= ppval(c_spline, xi_min); 
