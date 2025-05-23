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

AA = kron( kron( Ax, Ix ) , Ix ) + kron( kron( Ix, Ax ), Ix ) + kron( kron( Ix, Ix ), Ax );
I = kron( kron( Ix, Ix ), Ix );

%% 初始能量
[ G, dG, ddG] = Fun_Diff( u );
Ef(1,1) = TrapezFun( hx, G );
r(1,1) =  Ef(1,1) ;
E_inner(1,1) = -1/2 * TrapezFun( hx, (AA * u).*u ); % 内能
E(1,1) = E_inner(1,1) + Ef(1,1);
E_mod(1,1) = E(1,1);

%% n = 1 时
% 生成系数
for k = 1:1
    t(2,1)
    b( k, 1 ) = ( ( t( 1+1 ) - t( 1+1 - k ) )^(1-alp) - ( t( 1+1 ) - t( 1 + 2 - k ) )^(1-alp) )/( tau( 1 + 1 - k ) * gamma( 2 - alp ) );
end

% 生成右端项
F_1 = 0; % 源项

% 利用迭代法求t1时刻的预估解
u_hat = u(:,1); % 初始值
tol = 10^(-12); % 终止条件
N_iter = 10^(12);

for kk = 1:N_iter
    % 生成右端项
    [G_hat, dG_hat, ddG_hat] = Fun_Diff( u_hat );
    F_2 = b( 1, 1 ) * u(:,1) + c * ( S * u_hat + dG_hat );

    % 整体的右端项
    F = F_1 + F_2;

    % 解方程组
    u_update = ( ( b( 1, 1 ) + c * S ) * I - AA )\F;

    % 终止条件判断
    if abs( u_update - u_hat ) < tol
        break;
    end

    % 更新
    u_hat = u_update;
end

% 利用ESAV方法求t1时刻的解
% 辅助变量
[G_hat, dG_hat, ddG_hat] = Fun_Diff( u_hat );
Ef_hat = TrapezFun( hx, G_hat );
xi(2,1) = exp( r(1,1) )/exp( Ef_hat );
V = FunV( xi(2,1) );

% 非线性项
F = b( 1, 1 ) * u(:,1) + c * V * dG_hat + c * V * S * u_hat;

% 解线性方程组
u(:,2) = ( ( b( 1, 1 ) + c * V * S) * I - AA  )\F;

% 解r——线性方程
r(2,1) = r(1,1) - V * TrapezFun( hx, ( dG_hat - S * ( u(:,2) - u_hat ) ).*( u(:,2) - u(:,1) ) );

% 记录最大的解u
u_max(2,1) = max( abs( u(:,2) ) );

% 能量
[ G, dG, ddG] = Fun_Diff( u(:,2) );

% 离散的能量
Ef(2,1) = TrapezFun( hx, G ); % 非线性势能
E_inner(2,1) = -1/2 * TrapezFun( hx, (AA * u(:,2)).*u(:,2) ); % 内能
E(2,1) = E_inner(2,1) + Ef(2,1);
E_mod(2,1) = E_inner(2,1) + r(2,1);

%% n \geq 2 时；[0,T0]
for n = 2:N0
    t(n+1,1)

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
    u_upper = 1;   u_lower = -1; % MBP的界
    u_hat( u_hat > u_upper ) = u_upper;
    u_hat( u_hat < u_lower ) = u_lower; % 通过截断产生满足MBP的二阶解

    % 辅助变量
    [G_hat, dG_hat, ddG_hat] = Fun_Diff( u_hat );
    Ef_hat = TrapezFun( hx, G_hat );
    xi(n+1,1) = exp( r(n,1) )/exp( Ef_hat );
    V = FunV( xi(n+1,1) );

    % 非线性项
    F = F + c * V * dG_hat + c * V * S * u_hat;

    % 解线性方程组
    u(:,n+1) = ( ( b( 1, 1 ) + c * V * S) * I - AA  )\F;

    % 解r——线性方程
    r(n+1,1) = r(n,1) - V * TrapezFun( hx, ( dG_hat - S * ( u(:,n+1) - u_hat ) ).*( u(:,n+1) - u(:,n) ) );

    % 记录最大的解u
    u_max(n+1,1) = max( abs( u(:,n+1) ) );

    % 能量
    [ G, dG, ddG] = Fun_Diff( u(:,n+1) );

    % 离散的能量
    Ef(n+1,1) = TrapezFun( hx, G ); % 非线性势能
    E_inner(n+1,1) = -1/2 * TrapezFun( hx, (AA * u(:,n+1)).*u(:,n+1) ); % 内能
    E(n+1,1) = E_inner(n+1,1) + Ef(n+1,1);
    E_mod(n+1,1) = E_inner(n+1,1) + r(n+1,1);
end

%% n > T0 时；[T0,T]
n = n+1;
while n < Nt
    % 更新时间步长
    Par_tau(n,1) = ( ( E(n,1) - E(n-1,1) )/tau(n-1,1) )^2;
    tau_p = tau(n-1,1);
    tau(n,1) = max( tau_min, tau_max/( ( 1 + bet * Par_tau(n,1) )^(1/2) ) );

    % 判断是否终止
    if t(n,1) >= T
        %% 数据存储
        save E_mod_T.mat E_mod;
        save E_T.mat E;
        save u_max_T.mat u_max;
        save u_T.mat v;
        save t_T.mat t;
        save tau_T.mat tau;
        break;
    else
        t(n+1,1) = t(n,1) + tau(n,1);
        if t(n+1,1) > T
            tau(n,1) = T - t(n,1);
            t(n+1,1) = T;
        end
    end

    t(n+1,1)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
    u_upper = 1;   u_lower = -1; % MBP的界
    u_hat( u_hat > u_upper ) = u_upper;
    u_hat( u_hat < u_lower ) = u_lower; % 通过截断产生满足MBP的二阶解

    % 辅助变量
    [G_hat, dG_hat, ddG_hat] = Fun_Diff( u_hat );
    Ef_hat = TrapezFun( hx, G_hat );
    xi(n+1,1) = exp( r(n,1) )/exp( Ef_hat );
    V = FunV( xi(n+1,1) );

    % 非线性项
    F = F + c * V * dG_hat + c * V * S * u_hat;

    % 解线性方程组
    u(:,n+1) = ( ( b( 1, 1 ) + c * V * S) * I - AA  )\F;

    % 解r——线性方程
    r(n+1,1) = r(n,1) - V * TrapezFun( hx, ( dG_hat - S * ( u(:,n+1) - u_hat ) ).*( u(:,n+1) - u(:,n) ) );

    % 记录最大的解u
    u_max(n+1,1) = max( abs( u(:,n+1) ) );

    % 能量
    [ G, dG, ddG] = Fun_Diff( u(:,n+1) );

    % 离散的能量
    Ef(n+1,1) = TrapezFun( hx, G ); % 非线性势能
    E_inner(n+1,1) = -1/2 * TrapezFun( hx, (AA * u(:,n+1)).*u(:,n+1) ); % 内能
    E(n+1,1) = E_inner(n+1,1) + Ef(n+1,1);
    E_mod(n+1,1) = E_inner(n+1,1) + r(n+1,1);

    v = u(:,end);

    if t(n+1) >= 10 && t(n+1) <= 11
        save E_mod_T10.mat E_mod;
        save E_T10.mat E;
        save u_max_T10.mat u_max;
        save u_T10.mat v;
        save t_T10.mat t;
        save tau_T10.mat tau;
    end

    if t(n+1) >= 20 && t(n+1) <= 21
        save E_mod_T20.mat E_mod;
        save E_T20.mat E;
        save u_max_T20.mat u_max;
        save u_T20.mat v;
        save t_T20.mat t;
        save tau_T20.mat tau;
    end

    if t(n+1) >= 30 && t(n+1) <= 31
        save E_mod_T30.mat E_mod;
        save E_T30.mat E;
        save u_max_T30.mat u_max;
        save u_T30.mat v;
        save t_T30.mat t;
        save tau_T30.mat tau;
    end

    if t(n+1) >= 40 && t(n+1) <= 41
        save E_mod_T40.mat E_mod;
        save E_T40.mat E;
        save u_max_T40.mat u_max;
        save u_T40.mat v;
        save t_T40.mat t;
        save tau_T40.mat tau;
    end

    if t(n+1) >= 50 && t(n+1) <= 51
        save E_mod_T50.mat E_mod;
        save E_T50.mat E;
        save u_max_T50.mat u_max;
        save u_T50.mat v;
        save t_T50.mat t;
        save tau_T50.mat tau;
    end

   if t(n+1) >= 60 && t(n+1) <= 61
        save E_mod_T60.mat E_mod;
        save E_T60.mat E;
        save u_max_T60.mat u_max;
        save u_T60.mat v;
        save t_T60.mat t;
        save tau_T60.mat tau;
   end

   if t(n+1) >= 70 && t(n+1) <= 71
        save E_mod_T70.mat E_mod;
        save E_T70.mat E;
        save u_max_T70.mat u_max;
        save u_T70.mat v;
        save t_T70.mat t;
        save tau_T70.mat tau;
   end

    if t(n+1) >= 100 && t(n+1) <= 101
        save E_mod_T100.mat E_mod;
        save E_T100.mat E;
        save u_max_T100.mat u_max;
        save u_T100.mat v;
        save t_T100.mat t;
        save tau_T100.mat tau;
    end

    %%%%%%%%%%%%%%%%%%%%%%
    n = n+1;
end

%% 终止时间
t2 = clock;
UPC_times = etime(t2, t1)

%% 数值解
u_3 = reshape( u(:,end), Nx, Nx, Nx );
figure(2);
isosurface(X,Y,Z,u_3,0);
colormap( 'jet' );
% colorbar
axis([-0.5, 0.5, -0.5, 0.5, -0.5, 0.5])

%% 解的最大值
figure(3);
plot( [0,T]', [1,1]', 'r--', 'LineWidth', 1 )
hold on
plot( t, u_max, 'k-.', 'LineWidth', 1 )

%% 能量耗散律
figure(4);
plot( t, E, 'r-', 'LineWidth', 1 )
hold on
plot( t, E_mod, 'k-.', 'LineWidth', 1 )
