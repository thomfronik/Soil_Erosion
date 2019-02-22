function hmin = rh_condition_height(hplus, hr)

    hmin = -0.5*hplus + sqrt(hplus.^4 + 8*hr^3*hplus)./(2*hplus); 

end

