function hmin = rh_condition_height(hplus)
hmin = -0.5*hplus + sqrt(hplus.^4 + 8*hplus)./(2*hplus); 

