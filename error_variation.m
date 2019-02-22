etas = 1:1:5; 
errors = zeros(size(etas));

for i = 1:length(etas)
    eta = etas(i);
    errors(i) = numerical_solver(eta);
end
plot(etas,errors)

