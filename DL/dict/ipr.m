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

function D = ipr(Y, D, X, mutual_coh, Nit)

% Iterative projections and rotations, Algorithm 4.3 (Barchiesi & Plumbley 2013)
% To be run after each iteration of a DL algorithm (which updates D and X)
% Input:
%   Y            - signal matrix
%   D            - current dictionary
%   X            - current representation matrix
%   mutual_coh   - desired mutual coherence
%   Nit          - number of iterations (for projections and rotations)
% Output:
%   D            - updated dictionary

% BD 18.04.2017

if nargin < 5
  Nit = 5;	% number of iterations
end
[m,n] = size(D);

for it = 1 : Nit   % iterated projection and rotation operations

  % project Gram matrix on bounded coherences set
  G = D'*D;  % Gram matrix
  i = find(abs(G) > mutual_coh);
  G(i) = sign(G(i)) * mutual_coh;
  for j = 1:n
    G(j,j) = 1;
  end

  % spectral projection
  [V,S] = eig(G);
  [s,i] = sort(diag(S), 'descend');
  s(find(s<0)) = 0;
  D = diag(sqrt(s(1:m))) * V(:,i(1:m))';

  % rotation on the data set
  [U,~,V] = svd(D*X*Y');
  D = V*U'*D;
  D = normc(D);
end
