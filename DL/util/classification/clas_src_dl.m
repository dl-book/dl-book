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

function D = clas_src_dl(Y, H, nv, s)
% Design SRC classifier with dictionaries trained for each class (Section 8.3)
% Input:
%   Y     - signal matrix
%   H     - label matrix
%   nv    - number of atoms in each class dictionary (vector)
%   s     - sparsity level
% Output:
%   D     - concatenated dictionary for sparse representation

% BD 24.12.2017

[m,N] = size(Y);            % m=signal length, N=number of training signals
c = size(H,1);              % number of classes
Nv = sum(H,1);              % number of signals in each class

iters = 50;                 % number of iterations in DL algorithm
update = 'aksvd';           % dictionary update method
replatom = 'random';        % atom replacement strategy
params = {};

% Train a dictionary for each class
D = [];
for i = 1:c
  is = find(H(i,:));        % indices of signals for the current class
  D0 = randn(m,nv(i));      % initial dictionary
  D0 = normc(D0);
  [Di, X, errs] = DL(Y(:,is), D0, s, iters, str2func(update), params, 'replatom', replatom);
  D = [D Di];
end
