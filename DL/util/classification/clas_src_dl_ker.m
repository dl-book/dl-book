% Copyright (c) 2017 Bogdan Dumitrescu <bogdan.dumitrescu@acse.pub.ro>
% 
% Permission to use, copy, modify, and/or distribute this software for any
% purpose with or without fee is hereby granted, provided that the above
% copyright notice and this permission notice appear in all copies.
% 
% THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
% WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
% MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
% ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
% WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
% ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
% OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.

function A = clas_src_dl_ker(Y, H, nv, s, kerfun)
% SRC classifier with dictionaries trained for each class, using kernel AK-SVD
% Input:
%   Y        - signal matrix
%   H        - label matrix
%   nv       - number of atoms in each class dictionary (vector)
%   s        - sparsity level
%   kerfun   - handle of kernel function
% Output:
%   A     - concatenated dictionary for kernel sparse representation

% BD 25.12.2017

[m,N] = size(Y);            % m=signal length, N=number of training signals
c = size(H,1);              % number of classes

iters = 50;                 % number of iterations in DL algorithm

% Train a dictionary for each class
A = [];
for i = 1:c
  is = find(H(i,:));        % indices of signals for the current class
  A0 = randn(length(is),nv(i));      % initial dictionary
  A0 = normc(A0);           % cheap, but incorrect normalization; seems to work well
  [Ai, X, errs] = aksvd_ker(Y(:,is), A0, s, kerfun, iters);
  A = [A Ai];
end
