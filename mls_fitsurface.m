clc;clear;
% 定义域
l = -3; r = 3;
% 生成随机点 n 个，这 n 个点为已知数据
n = 200;
x = (r - l) * rand(1, n) + l;
y = (r - l) * rand(1, n) + l;
% 函数定义
z = 2*(1-x).^2 .* exp(-x.^2 - (y+1).^2) - 10*(x.^4/5 - y.^5) ...
.* exp(-x.^2 - y.^2) - 1/3*exp(-(x+1).^2 - y.^2);
% 点是随机的数据，meshgrid 生成点，也就是这里的划分网格
[xq, yq] = meshgrid(l:0.1:r, l:0.1:r);
% 插值生成对应的网格数据
zq = griddata(x, y, z, xq, yq);
% 绘制网格和点，这里的网格是插值生成的
figure;
subplot(2, 1, 1);
mesh(xq, yq, zq);
hold on
scatter3(x,y,z);
hold off
% 拟合上面的网格
points_sum = size(xq, 1);
simulated_z = zeros(points_sum, points_sum);
% 使用二维线性基[1, x, y] 作为基函数
m = 3;
% 定义格子范围（紧支）需要保证有三个不共线的点，否则A矩阵可能不可逆，这个值很难取以避免奇异矩阵的产生，这里简单的取了很大的值
% 我认为严谨的应该事先判断有几个点，对那些产生奇异的块适当的扩大紧支
smax = (r-l)/5;

% i,j 为拟合点的下标
for i = 1:points_sum
    for j = 1:points_sum
        x_val = xq(1, i); y_val = yq(j, 1);
        A = zeros(m, m);
        B = zeros(m, n);
        % 对每个点加权获得拟合值
        for k = 1:n
            xk = x(k); yk = y(k);
            w = w_func(((x_val-xk)^2 + (y_val - yk)^2)^0.5/smax);
            p = [1; xk; yk];
            A = A + w*(p*p');
            B(:, k) = w * p;
        end
        p_xy = [1, x_val, y_val];
        % 注意这里 x_val 是第 i 列， y_val是第 j 行 注意mesh的对应下标
        simulated_z(j, i) = p_xy * (A\B) * z';
    end
end

subplot(2, 1, 2);
surf(xq, yq, simulated_z);
hold on
scatter3(x,y,z);
hold off

% 这里的s就是半径相当于，也可以设置为格子
function [w] = w_func(s)
    if s <= 1/2
        w = 2/3 - 4*s^2 + 4*s^3;
    elseif s <= 1
        w = 4/3 - 4*s + 4*s^2 - 4/3*s^3;
    else
        w = 0;
    end
end
