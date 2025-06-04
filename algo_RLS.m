function [theta_est,MSD] = algo_RLS(theta,theta_true,para_num,lambda,input,d)
    N = length(d);
    P = 10^6 * eye(para_num,para_num);
    for k = 2*para_num-1:N
        x_n = input(k:-1:k-512+1); 
        if all(abs(x_n)<=0.5)
            MSD(k,1) = MSD(k-1);
        else
            e = d(k) - x_n' * theta;
            K = P * x_n / (lambda + x_n' * P * x_n); 
            theta = theta+ K * e; 
            P = (P - K * x_n' * P) / lambda;   
            MSD(k,1) = norm(theta-theta_true,2)^2;
       end
    end
    theta_est = theta;
    MSD = 10 * log10(MSD);
end