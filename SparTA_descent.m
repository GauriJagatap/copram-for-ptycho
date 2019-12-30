function [x] = SparTA_descent(x_init,y_abs,Phi_fn,Phit_fn,Nn,s,max_iter,tol1,tol2)

%% initialize parameters
N = sqrt(Nn);

[lRow lCol Nc] = size(y_abs);
error_hist(1,1) = 1;
y_abs2 = y_abs.^2; %quadratic measurements
phi = sqrt(sum(y_abs2,3)/Nc); %signal power

%initialize SPARTA parameters
mu = 9.5;
gamma = 30;

x = x_init;

%% start descent
    It_act = Phi_fn(x);
    It_mag = abs(It_act);
    y_ind = find(It_mag > y_abs/(1+gamma)); %assume gamma=inf
    
    y_est = Phi_fn(x);
    y_norm = y_est(y_ind)./abs(y_est(y_ind));
    y_err = zeros(size(y_abs));
    y_err(y_ind) = y_norm;
    sum_TAF = Phit_fn(y_est - y_abs.*y_err);

    grad_x = (mu/Nc)*sum_TAF;
    arg_TAF = x - grad_x;
    x = truncated_AF(arg_TAF,s);
    
end