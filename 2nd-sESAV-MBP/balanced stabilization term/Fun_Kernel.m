function [A] = Fun_Kernel(afa, sig, t, tau, rho, n)
%% 生成离散的卷积核

%% 高斯积分权重和系数
w1 = 0.555555555555;  w2 = 0.888888888888;  w3 = 0.555555555555;
s1 = -0.774596669241;  s2 = 0;  s3 = 0.774596669241;

%% 离散的系数
a(1,1) = 1/tau(n) * (sig * tau(n))^(1-afa)/gamma(2-afa);

for k = 1:n-1
    a(n-k+1,1) = 1/(gamma(2-afa)*tau(k)) * ((t(n)+sig*tau(n) - t(k))^(1-afa) - (t(n)+sig*tau(n) - t(k+1))^(1-afa));
    b(n-k,1) = 1/(tau(k)+tau(k+1)) * (w1*(s1*tau(k)/2)*Fun_Omega(t(n)+sig*tau(n)-(t(k)+tau(k)/2+s1*tau(k)/2),1-afa) + w3*(s3*tau(k)/2)*Fun_Omega(t(n)+sig*tau(n)-(t(k)+tau(k)/2+s3*tau(k)/2),1-afa));
%     b(n-k,1) = 1/(tau(k)*(tau(k)+tau(k+1))) * ( ( (t(n)+sig*tau(n)-t(k+1))^(2-afa)/gamma(3-afa) - (t(n)+sig*tau(n)-t(k))^(2-afa)/gamma(3-afa) )...
%         - ( (1/2*tau(k))*(t(n)+sig*tau(n)-t(k+1))^(1-afa)/gamma(2-afa) + (1/2*tau(k))*(t(n)+sig*tau(n)-t(k))^(1-afa)/gamma(2-afa) ) );
end

%% 生成离散的卷积核
A(1,1) = a(1,1) + rho(n-1)*b(1);
A(n,1) = a(n,1) - b(n-1,1);

for k = 2:n-1
    A(n-k+1,1) = a(n-k+1,1) + rho(k-1)*b(n-k+1) - b(n-k);
end

end 

