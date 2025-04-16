function [ FF, F, dF ] = Fun_Diff( u )
% 输出反应项的函数值及其导数值

theta = 0.8;
theta_c = 1.6;

FF = theta/2 * ( (1+u) .* log( 1 + u ) + ( 1 - u ) .* log( 1 - u ) ) - theta_c/2 * u.^2;
F =  theta/2 * ( log( (1 - u) ./ ( 1 + u ) ) ) + theta_c * u;
dF = - theta ./ ( 1 - u.^2 ) + theta_c;

end

