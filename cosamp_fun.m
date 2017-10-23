%2d cosamp_fun

function [xhat,xcosamp] = cosamp_fun(yy, Phi_f, PhiT_f, Nn, K, Its)

%---
N = sqrt(Nn);
%yy = yy(:); % 
[lRow lCol Nc]  = size(yy);
aa= zeros(N,N,Its); % stores current sparse estimate
s_cosamp = zeros(N);
kk=1; % current MP iteration
maxiter= 5;
verbose= 0;
tol= 1e-3;

while le(kk,Its),
    rr = yy - Phi_f(s_cosamp); %r = y - A(s) % 720x720 --> 512x512x100
    proxy = PhiT_f(rr); % At(r) % 512x512x100 --> 720x720
    %---Estimate support
    prox_vec = proxy(:);
    [tmp,ww]= sort(abs(prox_vec),'descend'); %refine this part of the code to store sorted mat indices
    tt= union(find(ne(s_cosamp(:),0)),ww(1:(2*K)));
    Kk = length(tt);

    %% Preparation for cg_solve %works on 3K dimensional vectors
    PP_tt = @(z) A_I(Phi_f,z,tt,N); %mod %in: 3K x 1 (720 x 720) ; out: 512x512x100
    
    %%
    PP_transpose_tt = @(z) A_I_transpose(PhiT_f,z,tt); %mod %out: 512x512x100; in: 3K x 1 (720 x 720)  
 
    qq = PP_transpose_tt(yy); %mod
    PPtranspose_PP_tt = @(z) PP_transpose_tt(PP_tt(z)); 
    [w, res, iter] = cgsolve(PPtranspose_PP_tt, qq, tol, maxiter, verbose);

    bb= 0*s_cosamp(:); bb(tt)= w;

    %---Prune
    kk = kk+1;
    [tmp,ww2]= sort(abs(bb),'descend'); 
    s_cosamp=0*bb;
    s_cosamp(ww2(1:K)) = bb(ww2(1:K));
    
    tt = ww2(1:K);
    
    aa(:,:,kk) = reshape(s_cosamp,N,N);  
    
    if kk>1
       if norm(aa(:,:,kk)-aa(:,:,kk-1))<1e-2*norm(aa(:,:,kk))
           break;
       end
    end
end
xcosamp= aa;
xcosamp(:,:,kk:end)=[];
xhat = xcosamp(:,:,end);

opts.debias = 0; %need to change
if opts.debias
    PP_tt = @(z) A_I(Phi_f,z,tt,N);
    PP_transpose_tt = @(z) A_I_transpose(PhiT_f,z,tt);
    qq = PP_transpose_tt(yy);
    PPtranspose_PP_tt = @(z) PP_transpose_tt(PP_tt(z));
    w = cgsolve(PPtranspose_PP_tt,qq, 1e-4, 100, 0);
    xhat = 0*s_cosamp; xhat(tt)= w;
end