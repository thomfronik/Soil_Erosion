function e = rh_condition_concentration(x0, x_L, x_R, y_L, y_R)

PP_L = spline(x_L, y_L);
PP_R = spline(x_R, y_R);
query_points_L = linspace(x_L(1)+x0,x_L(end), 100);
query_points_R = linspace(x_R(1)+x0,x_R(end), 100);
distance = abs(ppval(PP_L,query_points_L)-ppval(PP_R,query_points_R)); 
e = min(distance); 
disp(e)
end


