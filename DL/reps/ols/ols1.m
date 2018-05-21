% Copyright (c) 2009 Bogdan Dumitrescu <bogdan.dumitrescu@acse.pub.ro>
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


function [indv, x, A, y] = ols1(A, y, M, Mg)
%
% Orthogonal Least Squares (or Optimized Orthogonal Matching Pursuit)
% Input:
%  A      - matrix of basis vectors (not normed)
%  y      - vector to be matched
%  M      - number of nonzero coefficients (basis vectors used for matching)
%  Mg     - the first Mg vectors are selected without matching
% Output:
%  indv   - indices of nonzero coefficients in first M positions, the
%          remaining indices follow
%  x      - least squares solution for the best M positions

% BD 30.12.2009

if nargin < 4
  Mg = 0;
end

[N_Y,N] = size(A);  % N is the number of vectors in basis
indv = 1:N;  % the whole set of indices
vsnorm = sum(A.*A)';  % basis vectors squared norms
spr = A'*y; % scalar products with the residual
for k = 1 : M
  if k > Mg
    %[z,kb] = max(abs(A(k:end,k:end)'*y(k:end) ./ sqrt(vsnorm(k:end))));
    [z,kb] = max(abs(spr(k:end) ./ sqrt(vsnorm(k:end))));
    kb = kb+k-1;
    A(:,[k kb]) = A(:,[kb k]);  % bring best vector in current position
    indv([k kb]) = indv([kb k]);
    vsnorm([k kb]) = vsnorm([kb k]);
    spr([k kb]) = spr([kb k]);
  end
  
  % compute annihilating Householder reflector 
  %sk = norm(A(k:end,k));
  sk = sqrt(vsnorm(k));
  if A(k,k) < 0
    sk = -sk;
  end
  uk = A(k:end,k);
  uk(1) = uk(1) + sk;
  bk = sk * uk(1);
  
  % apply reflector
  A(k,k) = - sk;
  A(k+1:end,k) = 0;
  for j = k+1 : N
    tau = uk' * A(k:end,j) / bk;
    A(k:end,j) = A(k:end,j) - tau*uk;
  end

  % apply reflector to r.h.s
  tau = uk' * y(k:end) / bk;
  y(k:end) = y(k:end) - tau*uk;
  
  % update residual
  %y(k) = 0;
  
  % update norms and scalar product
  vsnorm(k+1:end) = vsnorm(k+1:end) - A(k,k+1:end)'.*A(k,k+1:end)';
  spr(k+1:end) = spr(k+1:end) - A(k,k+1:end)'*y(k);
end
% compute optimal solution
%x = A(1:M,1:M) \ y(1:M);    % it's a triangular system... !!!
opts.UT = true;
x = linsolve(A(1:M,1:M), y(1:M), opts);
