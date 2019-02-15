% Derive taylor expansions for expression involving theta - Vincent Meijer 

% Create symbolic variables
syms x y;
% Original second term
f = @(x) (1-x^(13/3))/((x^(4/3))*(x^3-1)); 
% Original Third term 
g = @(x) (x^(10/3) - 1)/(x^(1/3)*(x^3 - 1)); 
% Taylor of f 
f_t = taylor(f, x, 1, 'Order',2);
% Taylor of g
g_t = taylor(g, x, 1, 'Order', 2); 

% Make a slight substitution such that the  coefficients are correct
f_t = subs(f_t, x, y+1);
g_t = subs(g_t, x, y+1); 

% Display results 
disp("Coefficients for the second term");
coeffs(f_t)
disp("Coefficients for the third term");
coeffs(g_t)





