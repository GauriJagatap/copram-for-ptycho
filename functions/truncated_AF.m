function x_out = truncated_AF(x_in,K)
[n1 n2] = size(x_in);
x_in = x_in(:);
    [xx ii] = sort(abs(x_in(:)),'descend');
    x_out = zeros(size(x_in(:)));
    x_out(ii(1:K),1) = x_in(ii(1:K),1);
    x_out = reshape(x_out,n1,n2);
end