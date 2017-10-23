%2-D jsmp_fun.m

function [xhat,xjsmp] = jsmp_fun(yy, Phi_f, PhiT_f, Nn, K, Jh, Jb, Its)

%---
J = Jb*Jh;
N = sqrt(Nn); 

%% zero-pad such that the image size is exactly divisible by block size

[lRow lCol Nc]  = size(yy);
aa= zeros(N,N,Its); % stores current sparse estimate
s_jsmp = zeros(N);
num_blocks = round(Nn/J); %J = block len, K = block sparsity

kk=1; % current MP iteration
maxiter= 5;
verbose= 0;
tol= 1e-3;

while le(kk,Its),
    rr = yy - Phi_f(s_jsmp);
    proxy = PhiT_f(rr);

    %% reshape the blocks
    [proxy, indcss] = reshapeblock(proxy,Nn,Jb);
    
    %%
    proxy = proxy(:); %vectorizing column-wise
    %---Estimate support
   % pick the K largest blocks. Easy.
    proxy_jsmp_block = reshape(proxy, J, num_blocks);
    [trash,blockww] = sort(sum(proxy_jsmp_block.^2,1),'descend');
    newsupp = zeros(J,num_blocks);
    newsupp(:,blockww(1:(2*K))) = 1;
    newsupp = reshape(newsupp, Nn, 1);
    tt1 = find(newsupp==1);
    tt_mod = indcss(tt1);
    tt=union(find(ne(s_jsmp(:),0)),tt_mod);  
%     Kk = length(tt);

    %% Preparation for cg_solve
    PP_tt = @(z) A_I(Phi_f,z,tt,N);
    
    %% 
    PP_transpose_tt = @(z) A_I_transpose(PhiT_f,z,tt);
    qq = PP_transpose_tt(yy);
    PPtranspose_PP_tt = @(z) PP_transpose_tt(PP_tt(z));
    
    %Pseudo-inverse
    [w, res, iter] = cgsolve(PPtranspose_PP_tt, qq, tol, maxiter, verbose);

    bb1= 0*s_jsmp(:); bb1(tt)= w;
    bb1 = reshape(bb1,N,N);
    
    %%
    [bb1, idcss] = reshapeblock(bb1,Nn,Jb);    
    
    %% ---Prune
    kk = kk+1;
    
    bb1_block = reshape(bb1, J, num_blocks);
    [~,blockww] = sort(sum(bb1_block.^2,1),'descend');
    
    newsupp = zeros(J,num_blocks);
    newsupp(:,blockww(1:K)) = 1;
    newsupp = reshape(newsupp, Nn, 1);
    s_jp = bb1(:).*newsupp;
    s_jsmp(idcss) = s_jp;

%     s_jsmp=0*s_jsmp;

    aa(:,:,kk) = reshape(s_jsmp,N,N);
    
    if kk>1
       if norm(aa(:,:,kk)-aa(:,:,kk-1))<5e-3*norm(aa(:,:,kk))
           break;
       end
    end
end
xjsmp= aa;
xjsmp(:,:,kk+1:end)=[];
xhat = xjsmp(:,:,end);