function [CRLB_history] = compute_CRLB(input,para_num,sigma,p)
    compute_coefficient = @(sigma, p) (sigma^2*gamma(1/p) /  (p *(p-1)*gamma(1-1/p)) );
    compute_temp = @(H_old,h_new) (H_old*h_new);
    inv_recursive = @(H_old, h_new, temp) (H_old - temp*temp' / (1+h_new'*temp) );
    compute_CRLB1 = @(coefficient, H_old, h_new, temp) (coefficient*inv_recursive(H_old, h_new, temp)); 
    N = length(input);
    H = zeros(para_num,1);
    for j = 1:para_num-1
        H(:,j) = input(j+para_num-1:-1:j);
    end
     H_old = (H*H')^-1; 
    for k = 2*para_num-1:N
        x_n = input(k:-1:k-para_num+1);
        H(:,k-para_num+1) = x_n;
        coefficient = compute_coefficient(sigma,p);
        temp = compute_temp(H_old,x_n); 
        CRLB_temp = compute_CRLB1(coefficient,H_old,x_n,temp); 
        H_old = inv_recursive(H_old, x_n, temp);
        CRLB_history(k,1) = sum(diag(CRLB_temp));
    end
    CRLB_history = 10*log10(CRLB_history);
    is_diagonally_dominant = all(abs(diag(CRLB_temp)) > sum(abs(CRLB_temp - diag(diag(CRLB_temp))), 2));
end