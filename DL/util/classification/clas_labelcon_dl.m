% Copyright (c) 2018 Bogdan Dumitrescu <bogdan.dumitrescu@acse.pub.ro>
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

function [W, D, A] = clas_labelcon_dl(Y, H, Q, n, s, alpha, beta, init_method)
% Design label consistent DL classifier (Section 8.5.3)
% Input:
%   Y           - signal matrix
%   H           - label matrix
%   Q           - label consistency matrix
%   n           - number of atoms
%   s           - sparsity level
%   alpha       -
%   beta        - trade-off parameters in the objective
%   init_method - 0, random
%                 1, trained dictionaries for each class and LS classifier and A
% Output:
%   W           - classifier matrix
%   D           - dictionary for sparse representation
%   A           - label consistency transformation

% BD 22.04.2018

if nargin < 8
  init_method = 0;
end

Ye = [Y; sqrt(alpha)*H; sqrt(beta)*Q];    % extended signal matrix
iters = 50;                 % number of iterations in DL algorithm
iterc = 20;                 % number of iterations for training initial dictionaries
update = 'aksvd';           % dictionary update method
replatom = 'random';        % atom replacement strategy
params = {};
m = size(Y,1);              % signal length
c = size(H,1);              % number of classes
me = m+c+n;

% Label consistent DL
switch init_method
  case 0                    % plain initialization
    D0 = randn(me,n);
    D0 = normc(D0);
  case 1                    % train small dictionaries for each class
    D0 = zeros(me,n);
    jj = 0;
    for i = 1 : c           % train the dictionaries
      jc = find(H(i,:)==1); % indices of signals from class i
      Yc = Y(:, jc);        % those signals
      nc = sum(Q(:,jc(1))); % number of atoms for class i
      %D0c = randn(m,nc);
      D0c = Yc(:,1:nc);
      D0c = normc(D0c);
      Dc = DL(Yc, D0c, s, iterc, str2func(update), params, 'replatom', replatom);
      D0(1:m,jj+1:jj+nc) = Dc;
      jj = jj + nc;
    end
    X = omp(Y, D0(1:m,:), s);
    % compute initial classifier and label consistency transformation
    gamma = 1e-3;
    W = H*X' / (X*X' + gamma*eye(n));
    A = Q*X' / (X*X' + gamma*eye(n));
    D0(m+1:m+c,:) = sqrt(alpha) * W;
    D0(m+c+1:me,:) = sqrt(beta) * A;
    D0 = normc(D0);
end

De = DL(Ye, D0, s, iters, str2func(update), params, 'replatom', replatom);

% Extract dictionary, classifier matrix and label consistency transformation
D = De(1:m,:);
normD = sqrt(sum(D.*D)+eps);
D = normc(D);
W = De(m+1:m+c,:) / sqrt(alpha);
W = W ./ repmat(normD, c, 1);
A = De(m+c+1:me,:) / sqrt(beta);
A = A ./ repmat(normD, n, 1);