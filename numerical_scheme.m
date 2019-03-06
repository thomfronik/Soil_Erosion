% Numerical scheme for solving the coupled system of equations for h and c
% Considers one particle size only 

clear, close all

% Initialise parameters 
xr = 0;
eta_begin = 1; 
eta_end = 9; 
N_eta = 100; 

% Initialise constants and scaling factors
g = 9.81; 
n = 0.02;
h0 = 0.0895;
q0 = 0.1097;
v0 = 0.4; 
x0 = q0/v0;
delta = (x0*q0^2*n^2)/(h0^(13/3));
S0 = delta*v0*h0/q0; 

eta = linspace(eta_begin, eta_end, N_eta); 

% Loop to find best eta value 
global_min_err = inf;

errors = zeros(size(eta)); 
for i = 1:length(eta)

    % Initial values
    hr = (q0^2/(g*h0^3))^(1/3);
    %cr_upper = eta*delta*(1-1/(hr^(13/3)))+1/(hr^(13/3)); 
    %cr_lower = 1/(hr^(13/3))-(1/(108*eta*hr^(29/3)))*(10*delta*eta*hr-13)^2; 
    cr = eta(i)*delta*(1-hr^(-10/3))+ 1*hr^(-13/3);
    eta_lower = (cr-1/(hr^(13/3)))/(delta*(1-1/(hr^(13/3)))); 
    eta_upper = fsolve(@(x) 1/(hr^(13/3))-1/(108*x*hr^(29/3))*(10*x*delta*hr-13)^2 - cr, 100);


    % % % Analytical solution around hr % % %
    % Derivatives of c at xi_r and b0 and b1 from theta expansion

    c_prime_r = 1/(hr^(13/3))-cr; 
    b0 = (1/18)*(10*eta(i)*delta*hr-13)-(1/2)*sqrt((1/81)*(10*eta(i)*delta*hr-13)^2-(4/3)*eta(i)*(hr^(29/3))*c_prime_r); %maybe minus sign change
    c_double_prime_r = -13/(3*hr^(29/3)*eta(i))*b0 - c_prime_r; 
    b1 = ((1-(1/6)*eta(i)*hr^(29/3)*c_prime_r/(b0^2))^(-1))*(-(2/3)*c_prime_r*eta(i)*hr^(29/3)/b0 - (1/6)*(eta(i)*hr^(45/3))/(b0^2)*c_double_prime_r + 26/27 - (5/27)*delta*eta(i)*hr);

    % Analytical solutions
    h_expansion = @(x) hr + b0/(hr^(13/3)*eta(i))*(x-xr)+b0*b1/(2*hr^(29/3)*eta(i)^2)*(x-xr).^2; 
    c_expansion = @(x) cr + c_prime_r*(x-xr)+(1/2)*c_double_prime_r*(x-xr).^2;  

    h_prime = @(h) -(1./(h.^3 - hr.^3)).*(h.^3 .* ((1/eta(i))*(cr-h.^(-13/3)) - delta) + delta*h.^(-13/3));

    % Finding the interval for which we need an analytical expression

    % WARNING!: This solve crashes for a lot of eta values
    %h_begin = fsolve(h_prime, 1.0001*hr);
    %h_end = fsolve(h_prime,0.99999*hr); 

    % Very ugly way to determine where to start analytical solution
    right_test_h = linspace(hr, 0.5*hr, 1000); 
    left_test_h = linspace(hr, 2*hr, 1000); 
    right_test_h_prime = h_prime(right_test_h); 
    left_test_h_prime = h_prime(left_test_h); 
    right_neg_test_h_prime = right_test_h(right_test_h_prime < 0); 
    left_neg_test_h_prime = left_test_h(left_test_h_prime < 0); 
    h_begin = left_neg_test_h_prime(1); 
    h_end = right_neg_test_h_prime(1);

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
    hr_plus = h_M(end); 
    hr_min = h_M(1); 

    % % % Numerical integration of the system of equations % % %

    % Minimum and maximum xi coordinate respectively 
    x_min = -5; 
    x_max = 5; 

    % Number of points for which we want the solution to be evaluated
    N = 500; 

    % Intial values for the left and right parts of the solution respectively
    y0_L = [hr_min ; cr_min];
    y0_R = [hr_plus ; cr_plus];

    % Intervals on which we want to know the solutions
    x_R = linspace(xr_plus, x_max, N);
    x_L = linspace(xr_min, x_min, N);

    % Solving the coupled system of equations using ODE45
    [x_R, y_R] = ode45(@(t,x) system_one_particle_size(x, t, hr, eta(i), delta), x_R , y0_R);
    [x_L, y_L] = ode45(@(t,x) system_one_particle_size(x, t, hr, eta(i), delta), x_L, y0_L); 

    % % %  Satisfying the Rankine-Hugoniot conditions % % %

    % Constructing the full solution
    x = [ fliplr(x_L')' ; x_M(2:end-1) ; x_R ]; 
    h = [ fliplr(y_L(:,1)')' ; h_M(2:end-1) ; y_R(:,1) ]; 
    c = [ fliplr(y_L(:,2)')' ; c_M(2:end-1) ; y_R(:,2) ];
    
    if ~isreal(c)
        continue
    end
    
    % Constructing splines
    h_spline = spline(x,h); 
    c_spline = spline(x,c); 
    
    % Defining our h+ and c+ values and then use them to find h- and c-
    xplus = x(x > 0); 
    hplus = h(h < hr);
    cplus = c(h < hr); 
    hmin = rh_condition_height(hplus, hr); 
    [cmin, xmin] = cmin_finder(x, h, c, hmin); 

    % Finding jump position by finding the minimal error
    min_err = inf; 
    k = 0; 
    int_points = 100;
    
    % New way of finding the minimal error
    for j = 1:length(hplus)  
        [err_c, errint] = error_function(x,c, xmin(j), xplus(j), int_points);
        
        
        if err_c < 10e-3
            err = (1-(1/(xplus(j)-xmin(j)))*trapz([ fliplr(xmin(1:j)')' ; xplus(1:j) ],[ fliplr(cmin(1:j)')' ; cplus(1:j) ]))^2;
        else
            err = inf;
        end
     
        if err < min_err
            min_err = err;
            k = j;
        end

    end
    
    errors(i) = min_err; 
    if min_err < global_min_err
        
        global_min_err = min_err;
        opt_eta = eta(i);
        opt_c_spline = c_spline; 
        opt_x = x; 
        % Determining the full solution
        x_sol = [ fliplr(xmin(1:k)')' ; xplus(1:k) ];
        h_sol = [ fliplr(hmin(1:k)')' ; hplus(1:k) ]; 
        c_sol = [ fliplr(cmin(1:k)')' ; cplus(1:k) ]; 
        
    end

 end


%% Plotting the piecewise solutions
% entire plot can be replaced by "plot(x,h)" and "plot(x,c)" respectively,
% but this just shows nice colors. 
figure 
hold on
plot(eta, errors);
hold off 
xlabel('$\eta$', 'Interpreter', 'latex'); 
ylabel('error'); 

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

%% Plotting the full solution (1 period)

figure
plot(x_sol, h_sol); 
xlabel('$\xi$', 'Interpreter', 'latex') 
ylabel('h($\xi$)', 'Interpreter', 'latex') 

figure
plot(x_sol, c_sol); 
xlabel('$\xi$', 'Interpreter', 'latex') 
ylabel('c($\xi$)', 'Interpreter', 'latex') 
title('Concentration, one period'); 
%% Plot full solution

% Number of periods for which you want to see the solution
periods = 3; 

% Making the solution periodic
x_full = x_sol; 
h_full = h_sol; 
c_full = c_sol; 

for i = 1 : periods-1
    
    x_full = [x_full; x_sol+i*abs(x_sol(1)-x_sol(end))];
    h_full = [h_full; h_sol];
    c_full = [c_full; c_sol]; 
    
end

% Lenght of full interval 
L = x_full(end)-x_full(1);

% Determining underlying- and redeposited sediment bed shape
z_d = S0*(L-x_full); 
z_m = (1/opt_eta)*(-c_full + max(c_full)); %some integration constant in here 

% Relative heights 
sed_height = z_d + z_m; 
water_height = sed_height + h_full; 

% Plotting
figure
hold on
plot(x_full, z_d, 'k')
plot(x_full, sed_height ,'r')
plot(x_full, water_height, 'b')
x_between = [x_full, fliplr(x_full')'];
water = [sed_height, fliplr(water_height')'];
sediment = [z_d, fliplr(sed_height')']; 
bed = [zeros(size(z_d)), fliplr(z_d')']; 
%fill(x_between, water, [0 0.35 0.7]);
%fill(x_between, sediment, [0.6 0.3 0.15])
%fill(x_between, bed, [0.9 0 0.9]) 
xlabel('$\xi$', 'Interpreter', 'latex') 
ylabel('h($\xi$)', 'Interpreter', 'latex') 

figure
plot(x_full, c_full); 
xlabel('$\xi$', 'Interpreter', 'latex') 
ylabel('c($\xi$)', 'Interpreter', 'latex') 












