% Copyright (c) 2017-2018 Bogdan Dumitrescu <bogdan.dumitrescu@acse.pub.ro>
% Copyright (c) 2018 Paul Irofti <paul@irofti.net>
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

function [W, D] = clas_discrim_dl(Y, H, n, s, alpha, init_method)
% Design discriminative DL classifier (Section 8.5.2)
% Input:
%   Y           - signal matrix
%   H           - label matrix
%   n           - number of atoms
%   s           - sparsity level
%   alpha       - trade-off parameter in the objective
%   init_method - 0, random
%                 1, trained dictionary and LS classifier
% Output:
%   W           - classifier matrix
%   D           - dictionary for sparse representation

% BD 23.12.2017

if nargin < 6
  init_method = 0;
end

Ye = [Y; sqrt(alpha)*H];    % extended signal matrix
iters = 50;                 % number of iterations in DL algorithm
iterc = 20;                 % number of iterations for training initial dictionary
update = 'aksvd';           % dictionary update method
replatom = 'random';        % atom replacement strategy
params = {};
m = size(Y,1);              % signal length
c = size(H,1);              % number of classes
me = m+c;

% Discriminative DL
switch init_method
  case 0                    % plain initialization
    D0 = randn(me,n);
    D0 = normc(D0);
  case 1                    % trained dictionaries and LS classifier
    D0 = zeros(me,n);
    %D0c = randn(m,n);
    D0c = Y(:,1:n);
    D0c = normc(D0c);
    D0(1:m,:) = DL(Y, D0c, s, iterc, str2func(update), params, 'replatom', replatom);
    X = omp(Y, D0(1:m,:), s);
    gamma = 1e-3;
    W = H*X' / (X*X' + gamma*eye(n));     % compute initial classifier
    D0(m+1:me,:) = sqrt(alpha) * W;       % and append to dictionary
    D0 = normc(D0);
end

De = DL(Ye, D0, s, iters, str2func(update), params, 'replatom', replatom);

% Extract dictionary and classifier matrix
D = De(1:m,:);
normD = sqrt(sum(D.*D)+eps);
D = normc(D);
W = De(m+1:me,:) / sqrt(alpha);
W = W ./ repmat(normD, c, 1);

