function [ F ] = TrapezFun( hx, hy, u )
%% 数值积分——梯形公式
% 输入：hx, hy 网格长度；f 数值积分节点上的函数值

%f(N+1,1) = f(1,1);

Nx = sqrt(length(u));

f = reshape(u, Nx, Nx);
ff = ( f(1:end-1,:) + f(2:end,:) )/2;
fff = ( ff(:,1:end-1) + ff(:,2:end) )/2;

F = hx * hy * sum( sum( fff, 1 ), 2 );

end

