clc;
clear;

% 已知点数据
x = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
y = [0, 4, 5, 14, 15, 14.5, 14, 12, 10, 5, 4];
z = [1, 2, 1.5, 3, 2, 2.5, 3, 2.5, 2, 1.5, 1];

n = length(x);
[min_x, max_x] = deal(min(x), max(x));

% 绘制点的数量，也就是估计的点
points_num = 100;
simulated_x = linspace(min_x, max_x, points_num);
simulated_y = zeros(size(simulated_x));
simulated_z = zeros(size(simulated_x));

% 定义格子范围
smax = (max_x - min_x) * 4 / n;

% 使用一维一次基 [1, x] m = 2
% m = 2;
% 使用一维二次基 [1, x, x^2] m = 3
m = 3;

% 拟合每一个拟合点
for j = 1:length(simulated_x)
    x_val = simulated_x(j);
    
    % 预分配 A B，这里的p就是[1, x_val]
    A = zeros(m, m);
    B = zeros(m, n);
    
    % 计算 w 求和
    for i = 1:n
        xi = x(i);
        w = w_func(abs(x_val - xi) / smax);
        % p_i = [1; xi];
        p_i = [1; xi; xi^2];
        A = A + w * (p_i * p_i');
        B(:, i) = w * p_i;
    end
    
    % px = [1; x_val];
    px = [1; x_val; x_val^2];
    simulated_y(j) = px' * (A \ B) * y';
    simulated_z(j) = px' * (A \ B) * z';
end

% 绘制原始数据点和拟合曲线
figure;
scatter3(x, y, z, 'ro'); % 原始数据点
hold on;
plot3(simulated_x, simulated_y, simulated_z, 'b-'); % 拟合曲线
xlabel('X');
ylabel('Y');
zlabel('Z');
title('3D Curve Fitting');
grid on;
hold off;

function [w] = w_func(s)
    if s <= 1/2
        w = 2/3 - 4 * s^2 + 4 * s^3;
    elseif s <= 1
        w = 4/3 - 4 * s + 4 * s^2 - 4/3 * s^3;
    else
        w = 0;
    end
end
