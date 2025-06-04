function [theta_est,MSD] = algo_RZA_NLMS(theta,theta_true,para_num,rou,esilon,alpha,delta,input,d)
    N = length(d);
    for k = 2*para_num-1:N
        x_n = input(k:-1:k-para_num+1);
        if all(abs(x_n)<=0.5) 
           MSD(k) = MSD(k-1);
        else
           e = d(k) - x_n'*theta;
           theta = theta - rou*sign(theta)./(1+esilon*abs(theta)) + alpha*e*x_n / (x_n'*x_n + delta);
           MSD(k,1) = norm(theta-theta_true)^2;
        end
    end
    theta_est = theta;
    MSD = 10 * log10(MSD);
end