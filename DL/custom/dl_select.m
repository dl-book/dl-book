% Copyright (c) 2017-2018 Bogdan Dumitrescu <bogdan.dumitrescu@acse.pub.ro>
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

function [D, X] = dl_select(Y, A, n, s)
% Dictionary selection (from Cehver & Krause, 2011). Algorithm 5.2
% Input:
%   Y     - training signals
%   A     - initial atoms
%   n     - desired number of atoms
%   s     - desired sparsity level
% Output:
%   D     - learned dictionary
%   X     - final representation matrix

% BD 22.05.2017

%[m,N] = size(Y);
[~,na] = size(A);

id = [];    % indices of selected atoms
ia = 1:na;  % indices of available atoms
D = [];     % the dictionary

% greedy selection
for k = 1:n
  res = inf*ones(1,na);
  for i = ia  % find atom minimizing the representation error if added to dictionary
    Dnew = [D A(:,i)];
    X = omp(Y, Dnew, min(k,s));
    res(i) = norm(Y - Dnew*X, 'fro');
  end
  [~,j] = min(res);
  D = [D A(:,j)];
end
X = omp(Y,D,s);  % recompute final representation
