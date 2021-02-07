function [ ans ] = initial_value( lambda, eta, data, year,index)

t = index-1;
x = data(index);

ans = x/(exp(lambda*t)-((x*eta)*(exp(lambda*t)-1)));

disp(ans);

end

