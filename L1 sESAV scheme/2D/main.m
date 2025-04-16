clc,clear
%% 导入数据
data;

t1 = clock;

%% 生成质量矩阵和刚度矩阵
e = ones( Nx, 1 );
Ax = c * epsilon^2/(hx^2) * spdiags([e -2*e e], -1:1, Nx, Nx);
Ax(1,end) = c * epsilon^2/(hx^2);
Ax(end,1) = c * epsilon^2/(hx^2);
Ix = spdiags(e, 0, Nx, Nx);

AA = kron( Ax, Ix ) + kron( Ix, Ax );
I = kron( Ix, Ix );

%% 初始能量
[ G, dG, ddG] = Fun_Diff( u );
Ef(1,1) = TrapezFun( hx, hx, G );
r(1,1) =  Ef(1,1) ;
E_inner(1,1) = -1/2 * TrapezFun( hx, hx, (AA * u).*u ); % 内能
E(1,1) = E_inner(1,1) + Ef(1,1);
E_mod(1,1) = E(1,1);

%% n = 1 时
% 生成系数
for k = 1:1
    b( k, 1 ) = ( ( t( 1+1 ) - t( 1+1 - k ) )^(1-alp) - ( t( 1+1 ) - t( 1 + 2 - k ) )^(1-alp) )/( tau( 1 + 1 - k ) * gamma( 2 - alp ) );
end

% 生成右端项
F_1 = 0; % 时间分数阶导数离散的历史项

% 利用迭代法求t1时刻的预估解
u_hat = u(:,1); % 初始值
tol = 10^(-10); % 终止条件
N_iter = 10^(12);

for kk = 1:N_iter
    % 生成右端项
    [G_hat, dG_hat, ddG_hat] = Fun_Diff( u_hat );
    F_2 = b( 1, 1 ) * u(:,1) + c * ( S * u_hat + dG_hat );

    % 整体的右端项
    F = F_1 + F_2;

    % 迭代公式
    u_update = ( ( b( 1, 1 ) + c * S ) * I - AA )\F;

    % 终止条件判断
%     ccc = u_update - u_hat;
    if abs( u_update - u_hat ) < tol
        break;
    end

    % 更新
    u_hat = u_update;
end

% 利用ESAV方法求t1时刻的解
% 辅助变量
[G_hat, dG_hat, ddG_hat] = Fun_Diff( u_hat );
Ef_hat = TrapezFun( hx, hx, G_hat );
xi(2,1) = exp( r(1,1) )/exp( Ef_hat );
V = FunV( xi(2,1) );

% 非线性项
F = b( 1, 1 ) * u(:,1) + c * V * dG_hat + c * V * S * u_hat;

% 解线性方程组
u(:,2) = ( ( b( 1, 1 ) + c * V * S) * I - AA  )\F;

% 解r——线性方程
r(2,1) = r(1,1) - V * TrapezFun( hx, hx, ( dG_hat - S * ( u(:,2) - u_hat ) ).*( u(:,2) - u(:,1) ) );

% 记录最大的解u
u_max(2,1) = max( abs( u(:,2) ) );

% 能量
[ G, dG, ddG] = Fun_Diff( u(:,2) );

% 离散的能量
Ef(2,1) = TrapezFun( hx, hx, G ); % 非线性势能
E_inner(2,1) = -1/2 * TrapezFun( hx, hx, (AA * u(:,2)).*u(:,2) ); % 内能
E(2,1) = E_inner(2,1) + Ef(2,1);
E_mod(2,1) = E_inner(2,1) + r(2,1);

%% n \geq 2 时
for n = 2:Nt
    % 生成系数
    for k = 1:n
        b( k, 1 ) = ( ( t( n+1 ) - t( n+1 - k ) )^(1-alp) - ( t( n+1 ) - t( n + 2 - k ) )^(1-alp) )/( tau( n + 1 - k ) * gamma( 2 - alp ) );
    end

    % 生成右端项
    F = 0;
    for k = 1:n-1
        F = F + ( b(k,1) - b(k+1,1) ) * u(:,n+1-k);
    end
    F = F + b(n,1) * u(:,1);

    % 预估值 = 外推 + 截断
    u_hat = (tau(n) + tau(n-1) )/(tau(n-1)) *  u(:,n) - tau(n)/(tau(n-1)) *  u(:,n-1); % 外推
    u_upper = 0.9575;   u_lower = -0.9575; % MBP的界
%     u_upper = 1;   u_lower = -1; % MBP的界
    u_hat( u_hat > u_upper ) = u_upper;
    u_hat( u_hat < u_lower ) = u_lower; % 通过截断产生满足MBP的二阶解

    % 辅助变量
    [G_hat, dG_hat, ddG_hat] = Fun_Diff( u_hat );
    Ef_hat = TrapezFun( hx, hx, G_hat );
    xi(n+1,1) = exp( r(n,1) )/exp( Ef_hat );
    V = FunV( xi(n+1,1) );

    % 非线性项
    F = F + c * V * dG_hat + c * V * S * u_hat;

    % 解线性方程组
    u(:,n+1) = ( ( b( 1, 1 ) + c * V * S) * I - AA  )\F;

    % 解r——线性方程
    r(n+1,1) = r(n,1) - V * TrapezFun( hx, hx, ( dG_hat - S * ( u(:,n+1) - u_hat ) ).*( u(:,n+1) - u(:,n) ) );

    % 记录最大的解u
    u_max(n+1,1) = max( abs( u(:,n+1) ) );

    % 能量
    [ G, dG, ddG] = Fun_Diff( u(:,n+1) );

    % 离散的能量
    Ef(n+1,1) = TrapezFun( hx, hx, G ); % 非线性势能
    E_inner(n+1,1) = -1/2 * TrapezFun( hx, hx, (AA * u(:,n+1)).*u(:,n+1) ); % 内能
    E(n+1,1) = E_inner(n+1,1) + Ef(n+1,1);
    E_mod(n+1,1) = E_inner(n+1,1) + r(n+1,1);
end

%% 终止时间
t2 = clock;
UPC_times = etime(t2, t1)

%% 数值解
u_2 = reshape( u(:,end), Nx, Nx );
figure(1);
surf(X,Y,u_2)
shading interp
colormap( 'jet' );
colorbar
view([90, 90]);
axis off;

%% 解的最大值
figure(2);
plot( [0,T]', [1,1]', 'r--', 'LineWidth', 1 )
hold on
plot( t, u_max, 'k-.', 'LineWidth', 1 )

%% 能量耗散律
figure(3);
plot( t, E, 'r-', 'LineWidth', 1 )
hold on
plot( t, E_mod, 'k-.', 'LineWidth', 1 )

% %% 数据存储
% save E_mod_tau0002.mat E_mod;
% save E_tau0002.mat E;
% save u_max_tau0002.mat u_max;
% save u_2_tau0002_T20.mat u_2;
% save t_tau0002.mat t;
