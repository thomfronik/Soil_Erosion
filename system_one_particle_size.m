function y = system_one_particle_size(x,t,hr, eta, delta)

    h = -((x(1).^3-hr^3)^(-1))*(x(1).^3.*((1/eta)*(x(2)-1./(x(1).^(13/3))) - delta) + delta./(x(1).^(1/3)));
    c = 1./(x(1).^(13/3))-x(2); 
    y = [h ; c];

end
