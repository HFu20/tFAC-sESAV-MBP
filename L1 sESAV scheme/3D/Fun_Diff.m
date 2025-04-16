function [ FF, F, dF ] = Fun_Diff( u )
% 输出反应项的函数值及其导数值

FF = 1/4 * ( u.^2 - 1 ).^2;
F = u - u.^3;
dF = 1 - 3 * u.^2;

end

