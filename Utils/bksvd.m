function [U,S] = bksvd(A, k, iter, bsize)
%-------------------------------------------------------------------------------------
% Randomized Block Krylov Iteration for truncated Singular Value Decomposition
% Computes the left (column) singular vectors and corresponding values
% Described in Musco, Musco, 2015 (http://arxiv.org/abs/1504.05477)
%
% usage : 
%
%  input:
%  * A : matrix to decompose
%  * k : number of singular vectors to compute
%  * iter : number of iterations, default = 8
%  * bsize : block size, which must be >= k, default = k + 10
%
%
%  output:
%  k singular vector/value pairs. 
%  * U : a matrix whose columns are A's approximate top left singular vectors
%  * S : a diagonal matrix whose entries are A's approximate top singular values
%-------------------------------------------------------------------------------------

% Check input arguments and set defaults.
if nargin > 4
    error('bksvd:TooManyInputs','requires at most 4 input arguments');
end
if nargin < 2
    error('bksvd:TooFewInputs','requires at least 2 input arguments');
end
if nargin < 3
    iter = 8;
end
if nargin < 4
    bsize = k + 10;
end
if(k < 1 || iter < 1 || bsize < k)
    error('bksvd:BadInput','one or more inputs outside required range');
end

% Allocate space for Krylov subspace.
K = zeros(size(A,1),bsize*iter);
% Random block initialization.
G = randn(size(A,2),bsize);

% Construct and orthonormalize Krlov Subspace. 
% Orthogonalize at each step using economy size QR decomposition.
[block,R] = qr(A*G,0);
K(:,1:bsize) = block;
for i=2:iter
    [block,R] = qr(A*(A'*block),0);
    K(:,(i-1)*bsize+1:i*bsize) = block;
end
[Q,R] = qr(K,0);

% Rayleigh-Ritz postprocessing
M = Q'*A;
% Economy size dense SVD.
[U,S] = svd(M,0);
U = Q*U(:,1:k);
S = S(1:k,1:k);

end

%-------------------------------------------------------------------------------------
% Copyright (c) 2015 Christopher Musco, Cameron Musco
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.