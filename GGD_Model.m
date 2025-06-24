function output = GGD_Model(input,N,mu,sigma,p)
    % 定义广义高斯分布的概率密度函数
    ggd_pdf = @(x, mu, sigma, p) (p / (2 * sigma * gamma(1 / p))) * exp(-((abs(x - mu) / sigma) .^ p));
    % 数值积分计算广义高斯分布的CDF
    x_range = linspace(-10, 10, N);  % x值范围
    generalized_gaussian_cdf = zeros(1, N);
    for i = 1:N
        generalized_gaussian_cdf(i) = integral(@(x) ggd_pdf(x, mu, sigma, p), -Inf, x_range(i));
    end
    % histogram(generalized_gaussian_cdf);
    % 检查CDF是否为严格单调递增
    [~, unique_idx] = unique(generalized_gaussian_cdf);
    x_range = x_range(unique_idx);
    generalized_gaussian_cdf = generalized_gaussian_cdf(unique_idx);
    u_values = rand(N, 1); % 生成均匀分布的随机数
    w = interp1(generalized_gaussian_cdf, x_range, u_values, 'linear', 'extrap'); % 插值拟合
    histogram(w);
    output = input + w; % 算法所用数据（含GGD噪声）
end