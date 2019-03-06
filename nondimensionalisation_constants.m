% Numerical values of our small parameters 

% Data
g = 9.81;
h0 = 0.0895;
z0 = h0; 
q0 = 0.1097;
u0 = q0/h0;
n = 0.02; 
J = 4.79; 
F = 0.075;
mstar = 6; 
S = 0.0149; 
pw = 1000; 
ps = 2000; 
phim = 0.6; 
phib = 0.4; 
omegacr = 0.007; 
vi = [0.038 0.0137 0.0827 0.1369 0.2317]; 
pi = [0.376 0.234 0.2 0.166 0.0237]; 
v0 = sum(pi.*vi); 
x0 = q0/v0;

delta = (x0*q0^2*n^2)/(h0^(13/3));
S0 = delta*v0*h0/q0; 
omega0 = pw*g*S0*q0;
c0 = F*omega0*ps/(h0*v0*g*(ps-pw)); 
t0 = z0*ps*(1-phim)/(v0*c0);
A = h0*g*(ps-pw)/(J*ps); 
eps = x0/(u0*t0); 