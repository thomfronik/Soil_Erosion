% Numerical scheme for solving the coupled system of equations for h and c
% Considers one particle size only 

% Initialise constants and scaling factors
g = 9.81; 
n = 0.02;
h0 = 0.0895;
q0 = 0.1097;
v0 = 1; 
x0 = q0/v0;
delta = (x0*q0^2*n^2)/(h0^(13/3));

% Initial values
xr = 0;
hr = (q0^2/(g*h0^3))^(1/3);
cr = 1; 

%%
% Here we insert our analytical solution around hr
% For now its just an ugly perturbation

% Initialise parameters 
eta = 0.1; 

% Derivatives of c at xi_r and b0 and b1 from theta expansion
c_prime_r = 1/(hr^(13/3))-cr; 
b0 = (1/18)*(10*eta*delta*hr-13)-(1/2)*sqrt((1/81)*(10*eta*delta*hr-13)^2-(4/3)*eta*hr^(29/3)*c_prime_r); %maybe minus sign change
c_double_prime_r = -13/(3*hr^(29/3)*eta)*b0 - c_prime_r; 
b1 = ((1-(1/6)*eta*hr^(29/3)*c_prime_r/(b0^2))^(-1))*(-(2/3)*c_prime_r*eta*hr^(29/3)/b0 - (1/6)*(eta*hr^(45/3))/(b0^2)*c_double_prime_r + 26/27 - (5/27)*delta*eta*hr);

h_expansion = @(x) hr + b0/(hr^(13/3)*eta)*(x-xr)+b0*b1/(2*hr^(29/3)*eta^2)*(x-xr).^2; 

% Boundaries of analytical solution
M = 10;
epsilon = 10^-4;
xr_plus = xr + epsilon;
xr_min = xr - epsilon; 
x_M = linspace(xr_min, xr_plus, M); 

h_M = h_expansion(x_M); 

cr_plus = cr;
cr_min = cr; 
hr_plus = h_M(end);%hr - epsilon;
hr_min = h_M(1);%hr + epsilon;

%% Numerical integration of the system of equations 
% Minimum and maximum xi coordinate respectively 
x_min = xr - 5*x0; 
x_max = xr + 5*x0; 

% Number of points for which we want the solution to be evaluated
N = 300; 

% Intial values for the left and right parts of the solution respectively
y0_L = [hr_min ; cr_min];
y0_R = [hr_plus ; cr_plus];

% Intervals on which we want to know the solutions
x_R = linspace(xr_plus, x_max, N);
x_L = linspace(xr_min, x_min, N);

% Solving the coupled system of equations using ODE45
[x_R, y_R] = ode45(@(t,x) system_one_particle_size(x, t, hr, eta, delta), x_R , y0_R);
[x_L, y_L] = ode45(@(t,x) system_one_particle_size(x, t, hr, eta, delta), x_L, y0_L); 

% Plotting the piecewise solutions
figure
hold on
plot(x_L, real(y_L(:,1)))
plot(x_M, real(h_M))
plot(x_R, real(y_R(:,1)))
xlabel('$\xi$', 'Interpreter', 'latex') 
ylabel('h($\xi$)', 'Interpreter', 'latex') 

figure
hold on 
plot(x_L, real(y_L(:,2)))
plot(x_R, real(y_R(:,2)))
xlabel('$\xi$', 'Interpreter', 'latex') 
ylabel('c($\xi$)', 'Interpreter', 'latex') 


