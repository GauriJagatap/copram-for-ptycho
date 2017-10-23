function [proxy, indcss] = reshapeblock(proxy,Nn,Jb)
N = sqrt(Nn);
[rr, cc] = size(proxy);
indcs = reshape((1:Nn),N,N);
proxxy = zeros(Jb*rr,cc/Jb);    
indcss = zeros(Jb*rr,cc/Jb);
for i=1:Jb:cc
    for j=1:rr
        proxxy(Jb*j-Jb+1:Jb*j,ceil(i/Jb)) = proxy(j,i:i+Jb-1);
        indcss(Jb*j-Jb+1:Jb*j,ceil(i/Jb)) = indcs(j,i:i+Jb-1);
    end
end
proxy = proxxy;
indcss = indcss(:);
end