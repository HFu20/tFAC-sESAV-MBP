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
% 生成离散卷积核
b(1,1) = 1/tau(1) * (sig * tau(1))^(1-alp)/gamma(2-alp);

% 生成右端项
F_1 = 0; % 源项

%% 系数矩阵
B_left = ( b(1,1) + S * sig * c ) * I - sig * AA;
B_right = ( b(1,1) - S * ( 1 - sig ) * c ) * I + ( 1 - sig ) * AA;

% 迭代法求解非线性方程
u_star = u(:,1); % 初始值
tol = 10^(-12); % 终止条件
N_iter = 10^(12);

for kk = 1:N_iter
    % 生成右端项
    u_hat = sig * u_star + ( 1 - sig ) * u(:,1);
    [G_hat, dG_hat, ddG_hat] = Fun_Diff( u_hat );
    F_2 = B_right * u(:,1) + c * ( S * u_hat + dG_hat );

    % 整体的右端项
    F = F_1 + F_2;

    % 解方程组
    u_update = B_left\F;

    % 终止条件判断
    if abs( u_update - u_star ) < tol
        break;
    end

    % 更新
    u_star = u_update;
end

% 预估解
u_hat = sig * u_update + (1-sig) * u(:,1);

%% 利用ESAV方法求t1时刻的解
% 辅助变量
[G_hat, dG_hat, ddG_hat] = Fun_Diff( u_hat );
Ef_hat = TrapezFun( hx, hx, G_hat );
xi(2,1) = exp( r(1,1) )/exp( Ef_hat );
V = FunV( xi(2,1) );

% 系数矩阵
B_left = ( b(1,1) + S * sig * c * V ) * I - sig * AA;
B_right = ( b(1,1) - S * ( 1 - sig ) * c * V ) * I + ( 1 - sig ) * AA;

% 右端项
F_2 = B_right * u(:,1) + c * V * ( S * u_hat + dG_hat );
F = F_1 + F_2;

% 解u
u(:,2) = B_left\F;
u_max(2,1) = max( abs( u(:,2) ) );

% 解r——线性方程
r(2,1) = r(1,1) - TrapezFun( hx, hx, ( dG_hat ).*( u(:,2) - u(:,1) ) )...
         + TrapezFun( hx, hx, ( S * ( sig * u(:,2) + (1-sig) * u(:,1) - u_hat ) ).*( u(:,2) - u(:,1) ) );

% 能量
[ G, dG, ddG] = Fun_Diff( u(:,2) );

% 离散的能量
Ef(2,1) = TrapezFun( hx, hx, G ); % 非线性势能
E_inner(2,1) = -1/2 * TrapezFun( hx, hx, (AA * u(:,2)).*u(:,2) ); % 内能
E(2,1) = E_inner(2,1) + Ef(2,1);
E_mod(2,1) = E_inner(2,1) + r(2,1);

%% n \geq 2 时
for n = 2:Nt
    % 生成离散的卷积核
    A = Fun_Kernel(alp, sig, t, tau, rho, n);

    % 生成右端项
    t_hat = sig * t(n+1,1) + ( 1 - sig ) * t(n,1);
    F = 0;
    for k = 1:n-1
        F = F - A(n-k+1,1)*(u(:,k+1)-u(:,k));
    end

    % 预估值 = 外推 + 截断
    u_hat = ( sig * tau(n) + tau(n-1) )/(tau(n-1)) *  u(:,n) - sig * tau(n)/(tau(n-1)) *  u(:,n-1); % 外推
    u_upper = 0.9575;   u_lower = -0.9575; % MBP的界
    u_hat( u_hat > u_upper ) = u_upper;
    u_hat( u_hat < u_lower ) = u_lower; % 通过截断产生满足MBP的二阶解

    % 辅助变量
    [G_hat, dG_hat, ddG_hat] = Fun_Diff( u_hat );
    Ef_hat = TrapezFun( hx, hx, G_hat );
    xi(n+1,1) = exp( r(n,1) )/exp( Ef_hat );
    V = FunV( xi(n+1,1) );

    % 非线性项
    F = F + c * V * dG_hat + c * V * S * u_hat;

    % 系数矩阵
    B_left = ( A(1,1) + S * sig * c * V ) * I - sig * AA;
    B_right = ( A(1,1) - S * ( 1 - sig ) * c * V ) * I + ( 1 - sig ) * AA;

    F = F + B_right * u(:,n);

    % 解线性方程组
    u(:,n+1) = B_left\F;

    % 解r——线性方程
    r(n+1,1) = r(n,1) - V * TrapezFun( hx, hx, ( dG_hat - S * ( sig * u(:,n+1) + (1-sig) * u(:,n) - u_hat ) ).*( u(:,n+1) - u(:,n) ) );

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
figure(2);
surf(X,Y,u_2)
shading interp
colormap( 'jet' );
colorbar
view([90, 90]);
axis off;

%% 解的最大值
figure(3);
plot( [0,T]', [0.9575,0.9575]', 'r--', 'LineWidth', 1 )
hold on
plot( t, u_max, 'k-.', 'LineWidth', 1 )

%% 能量耗散律
figure(4);
plot( t, E, 'r-', 'LineWidth', 1 )
hold on
plot( t, E_mod, 'k-.', 'LineWidth', 1 )

% save u_tau01_VV.mat u_2;
save t_tau003_VV.mat t;
save u_max_tau003_VV.mat u_max;