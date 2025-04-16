clc,clear

%% 模型参数
alp = 0.9;
c = 1; % 迁移率
epsilon = 0.01; % 界面宽度

%% 离散点
sig = 1 - alp/2;

%% 求解区域
T = 30;
a = 0;  b = 1;

%% 时间网格剖分参数
Nt = 3 * 10^(1) * T;
gam = 1;
% gam = 2/sigma; % 网格分级参数

%% 时间网格剖分
T0 = min( [ 1/gam, T ]' );
T0 = T;
N0 = ceil( Nt/( T + 1 - 1/gam ) );
N0 = Nt;

% [0,T0]——分级网格
t0 = T0 * ((0:1:N0)'/N0).^(gam);  % 时间分级网格
tau0 = diff( t0 );  % 分级网格步长

% [T0, T]——分级网格
N1 = Nt - N0;
%s_t = rand(N1,1);
s_t = ones(N1,1);
tau1 = (T - T0) * s_t/(sum(s_t));

tau = [tau0; tau1];
t = t0;
for k = 1:N1
    t(end+1,1) = t(end,1) + tau1( k, 1 );
end

for k = 1:Nt-1
    rho(k,1) = tau(k,1)/tau(k+1,1);
end

%% 空间网格剖分
Nx = 128;  % 空间网格剖分次数
hx = (b-a)/Nx;  % mesh size
x = (a+hx:hx:b)';  % mesh grid

%% 稳定化
S = 8.02; % 稳定化常数

%% 问题初值
uu = load('u_intal_08.mat');
u_2 = uu.u_2;
u = reshape( u_2, Nx^2, 1 );
u_max(1,1) = max( abs( u(:,1) ) );

[X, Y] = meshgrid(x, x);
% figure(1);
% surf(X,Y,u_2)
% shading interp
% colormap( 'jet' );
% colorbar
% view([90, 90]);
% axis off;