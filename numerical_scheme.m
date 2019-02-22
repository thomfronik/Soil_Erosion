% Numerical scheme for solving the coupled system of equations for h and c
% Considers one particle size only 

clear, close all

% Initialise parameters 
eta = 0.01;
xr = 0;

% Initialise constants and scaling factors
g = 9.81; 
n = 0.02;
h0 = 0.0895;
q0 = 0.1097;
v0 = 0.38; 
x0 = q0/v0;
delta = (x0*q0^2*n^2)/(h0^(13/3));

% Initial values
hr = (q0^2/(g*h0^3))^(1/3);
%cr_upper = eta*delta*(1-1/(hr^(13/3)))+1/(hr^(13/3)); 
%cr_lower = 1/(hr^(13/3))-(1/(108*eta*hr^(29/3)))*(10*delta*eta*hr-13)^2; 
cr = eta*delta*(1-hr^(-10/3))+ 1*hr^(-13/3);
eta_lower = (cr-1/(hr^(13/3)))/(delta*(1-1/(hr^(13/3)))); 
eta_upper = fsolve(@(x) 1/(hr^(13/3))-1/(108*x*hr^(29/3))*(10*x*delta*hr-13)^2 - cr, 100);

%%
% Here we insert our analytical solution around hr
% Derivatives of c at xi_r and b0 and b1 from theta expansion

c_prime_r = 1/(hr^(13/3))-cr; 
b0 = (1/18)*(10*eta*delta*hr-13)-(1/2)*sqrt((1/81)*(10*eta*delta*hr-13)^2-(4/3)*eta*(hr^(29/3))*c_prime_r); %maybe minus sign change
c_double_prime_r = -13/(3*hr^(29/3)*eta)*b0 - c_prime_r; 
b1 = ((1-(1/6)*eta*hr^(29/3)*c_prime_r/(b0^2))^(-1))*(-(2/3)*c_prime_r*eta*hr^(29/3)/b0 - (1/6)*(eta*hr^(45/3))/(b0^2)*c_double_prime_r + 26/27 - (5/27)*delta*eta*hr);

% Analytical solutions
h_expansion = @(x) hr + b0/(hr^(13/3)*eta)*(x-xr)+b0*b1/(2*hr^(29/3)*eta^2)*(x-xr).^2; 
c_expansion = @(x) cr + c_prime_r*(x-xr)+(1/2)*c_double_prime_r*(x-xr).^2;  

h_prime = @(h) -(1./(h.^3 - hr.^3)).*(h.^3 .* ((1/eta)*(cr-h.^(-13/3)) - delta) + delta*h.^(-13/3));

% Finding the interval for which we need an analytical expression
h_begin = fsolve(h_prime, 1.001*hr);
h_end = fsolve(h_prime,0.99999*hr); 
x_begin = 1.05*(fsolve(@(x) h_expansion(x) - h_begin, 0));
x_end = 1.05*(fsolve(@(x) h_expansion(x) - h_end, 0));

% Boundaries of analytical solution
M = 10;

epsilon = 10^-2;
xr_plus = x_end;
xr_min = x_begin; 
x_M = linspace(xr_min, xr_plus, M)'; 

h_M = h_expansion(x_M); 
c_M = c_expansion(x_M); 

cr_plus = c_M(end);
cr_min = c_M(1); 
hr_plus = h_M(end); %hr - epsilon;
hr_min = h_M(1); %hr + epsilon;

%% Numerical integration of the system of equations 
% Minimum and maximum xi coordinate respectively 
x_min = xr - 10*eta*x0; 
x_max = xr + 10*eta*x0; 

% Number of points for which we want the solution to be evaluated
N = 500; 

% Intial values for the left and right parts of the solution respectively
y0_L = [hr_min ; cr_min];
y0_R = [hr_plus ; cr_plus];

% Intervals on which we want to know the solutions
x_R = linspace(xr_plus, x_max, N);
x_L = linspace(xr_min, x_min, N);

% Solving the coupled system of equations using ODE45
[x_R, y_R] = ode45(@(t,x) system_one_particle_size(x, t, hr, eta, delta), x_R , y0_R);
[x_L, y_L] = ode45(@(t,x) system_one_particle_size(x, t, hr, eta, delta), x_L, y0_L); 

%% Plotting the piecewise solutions
% entire plot can be replaced by "plot(x,h)" and "plot(x,c)" respectively,
% but this just shows nice colors. 

figure
hold on 
plot(x_L, real(y_L(:,1)))
plot(x_M, real(h_M))
plot(x_R, real(y_R(:,1)))
hold off
xlabel('$\xi$', 'Interpreter', 'latex') 
ylabel('c($\xi$)', 'Interpreter', 'latex') 

figure
hold on 
plot(x_L, real(y_L(:,2)))
plot(x_M, real(c_M))
plot(x_R, real(y_R(:,2)))
hold off
xlabel('$\xi$', 'Interpreter', 'latex') 
ylabel('c($\xi$)', 'Interpreter', 'latex') 

%% Satisfying the Rankine-Hugoniot conditions

% Constructing the full solution
x = [fliplr(x_L')'; x_M(2:end-1); x_R]; 
h = [fliplr(y_L(:,1)')'; h_M(2:end-1); y_R(:,1)]; 
c = [fliplr(y_L(:,2)')'; c_M(2:end-1); y_R(:,2)];

% Constructing splines
h_spline = spline(x,h); 
c_spline = spline(x,c); 

% Defining our h+ and c+ values and then use them to find h- and c-
xplus = x(x > 0); 
hplus = h(h < hr);
cplus = c(h < hr); 
hmin = rh_condition_height(hplus, hr); 
[cmin, xmin] = cmin_finder(h_spline, c_spline, hmin); 


% Finding jump position by finding the minimal error
min_err = inf; 
k = 0; 
int_points = 100; 

for i = 1:length(hplus)

    err = error_function(c_spline, xmin(i), xplus(i), int_points);
    
    if err < min_err
        min_err = err;
        k = i;
    end
    
end

% Jump coordinate
xi_plus = xplus(k);
xi_min = xmin(k); 

%% Plotting the full solution 

x_sol = [ fliplr(xmin(1:k)')' ; xplus(1:k) ];
h_sol = [ fliplr(hmin(1:k)')' ; hplus(1:k) ]; 
c_sol = [ fliplr(cmin(1:k)')' ; cplus(1:k) ]; 

figure
plot(x_sol, h_sol); 
xlabel('$\xi$', 'Interpreter', 'latex') 
ylabel('h($\xi$)', 'Interpreter', 'latex') 

figure
plot(x_sol, c_sol); 
xlabel('$\xi$', 'Interpreter', 'latex') 
ylabel('c($\xi$)', 'Interpreter', 'latex') 





