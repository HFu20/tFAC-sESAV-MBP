clc,clear

%% 模型参数
alp = 0.5;
c = 1; % 迁移率
epsilon = 0.03; % 界面宽度

%% 求解区域
T = 100;
a = -0.5;  b = 0.5;

%% 时间网格剖分参数
gam = (2-alp)/alp; % 网格分级参数

%% 时间网格剖分
T0 = 0.5;
N0 = 20;

% [0,T0]——分级网格
t = T0 * ((0:1:N0)'/N0).^(gam);  % 时间分级网格
tau = diff( t );  % 分级网格步长

% [T0,T]——自适应网格
tau_min = 0.01;
tau_max = 1;
bet = 10^7;
Nt = 10^13;

%% 空间网格剖分
Nx = 80;  % 空间网格剖分次数
hx = (b-a)/Nx;  % mesh size
x = (a+hx:hx:b)';  % mesh grid
y = x;
z = x;

%% 稳定化
S = 2; % 稳定化常数

%% 问题初值
[X, Y, Z] = meshgrid(x, y, z);

u1 = tanh( ( 0.2 - sqrt( ( X - 0.14 ).^2 + Y.^2 + Z.^2 ) ) / epsilon );
u2 = tanh( ( 0.2 - sqrt( ( X + 0.14 ).^2 + Y.^2 + Z.^2 ) ) / epsilon );
u_3 = max( u1, u2 );

save u_T0.mat u_3;

u = reshape( u_3, Nx^3, 1 );
u_max(1,1) = max( abs( u(:,1) ) );

figure(1);
isosurface(X,Y,Z,u_3,0);
colormap( 'jet' );
% colorbar
axis([a, b, a, b, a, b])
% axis off
% axis equal;

