function fittedCurve = mls_fitcurve(x, y)
    % 曲线需要横坐标单调
    % clc;clear;
    % 已知点数据
    % x = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0];
    % y = [0, 4, 5, 14, 15, 14.5, 14, 12, 10, 5, 4];
    n = size(x,2);
    [min_x, max_x] = deal(min(x), max(x));
    rowrank = randperm(size(x, 2));
    x = x(rowrank);
    y = y(rowrank);
    % 绘制点的数量，也就是估计的点
    points_num = 1000;
    simulated_x = linspace(min_x, max_x, points_num);
    simulated_y = zeros(1, points_num);
    % 定义格子范围（紧支）需要保证有三个不共线的点，否则A矩阵可能不可逆，这里取4个点范围，避免共线
    % 这里是比较均匀的，而且点比较多，到曲面就很难选取了
    smax = (max_x-min_x)*4/n;
    % 使用一维二次基 [1, x, x^2] m = 3
    m = 3;
    
    % 拟合每一个拟合点
    for j = 1:points_num
        x_val = simulated_x(j);
        % 预分配 A B, 这里的p就是[1，x_val，x_val^2]
        A = zeros(m, m);
        B = zeros(m, n);
        % 计算 w 求和
        for i = 1:n
            xi = x(i);
            w = w_func(abs(x_val - xi)/smax);
            p_i = [1;xi;xi^2];
            A = A + w * (p_i*p_i');
            B(:, i) = w*p_i;
        end
        px = [1, x_val, x_val^2];
        simulated_y(j) = px * (A\B) * y';
    end
    fittedCurve = simulated_y;
    plot(x, y,'--o', simulated_x, simulated_y);

    function [w] = w_func(s)
        if s <= 1/2
            w = 2/3 - 4*s^2 + 4*s^3;
        elseif s <= 1
            w = 4/3 - 4*s + 4*s^2 - 4/3*s^3;
        else
            w = 0;
        end
    end
end