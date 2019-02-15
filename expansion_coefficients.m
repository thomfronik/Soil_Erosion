% Returns expansion coefficients - Vincent Meijer
% NOTE: Can probably be made quicker but currently functional 

function [b0_num, b1_num] = expansion_coefficients(h_r_i, eta_i, delta_i, sigma_1_i, sigma_2_i)

% Initialize a lot of symbols
syms eta delta h_r a0 b0 a1 b1 sigma_1 sigma_2;
% Sigma_1 = sum_i v_i c_i'(xi_r)
% Sigma_2 = sum_i v_i c_i''(xi_r) 

% Expressions for b0, b1 

eqn1 = b0 == ((h_r^(13/3))/3) * a0 + (10/9)*eta*delta*h_r - 13/9;
eqn2 = b1 == ((h_r^(13/3))/3) *(a1 + 2*a0) + 26/27 - (5/27)*eta*delta*h_r; 

% Sub in expression for a0
eqn3 = subs(eqn1, a0, -((eta*h_r^(16/3))/b0) * sigma_1);
% Solve for b0
b0_sol = solve(eqn3, b0); 
% Sub in a0
eqn4 = subs(eqn2, a0, -((eta*h_r^(16/3))/b0) * sigma_1);
% Sub in a1
eqn4 = subs(eqn4, a1, ((-b1*eta*h_r^(16/3))/(2*b0^2)) * sigma_1 - ((eta*h_r^(32/3))/(2*b0^2)) * sigma_2);


% Now sub in numerical values
b0_num = vpa(min(subs(b0_sol, {h_r, eta, delta, sigma_1}, {h_r_i, eta_i, delta_i, sigma_1_i})));
% Sub in b0 
eqn4 = subs(eqn4, b0, b0_num);
% Solve 
b1_sol = solve(eqn4(1), b1);
b1_num = vpa(subs(b1_sol, {h_r, eta, delta, sigma_1,sigma_2}, {h_r_i, eta_i, delta_i, sigma_1_i,sigma_2_i}));




