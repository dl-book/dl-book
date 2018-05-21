% Copyright (c) 2018 Bogdan Dumitrescu <bogdan.dumitrescu@acse.pub.ro>
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

function [X, shared] = omp_mask(Y, D, s, shared, varargin)

% Orthogonal Matching Pursuit algorithm for data with missing entries.
% The missing entries are simply ignored (masked). 
% Input:
%   Y     - training signals set
%   Ymask - varargin{1} = signal mask
%   D     - current dictionary
%   s     - sparsity level
%
% Output:
%   X     - sparse representations

% BD 10.04.2018

[~,N] = size(Y);

Ymask = varargin{1};
Y = Y .* Ymask;    % apply mask, to be sure

ompparams = {'checkdict', 'off'};

for i = 1:N    % represent signals one by one, since the available entries are different
  suppy = find(Ymask(:,i));
  Dsmall = D(suppy,:);
  normDsmall = sqrt( sum( Dsmall.*Dsmall ) );
  Dsmall = Dsmall ./ repmat(normDsmall, size(Dsmall,1), 1);
%  X(:,i) = omp(Y(suppy,i), Dsmall, s, shared, varargin);
  X(:,i) = omp_sparse(Dsmall'*Y(suppy,i), Dsmall'*Dsmall, s, ompparams{:});
  X(:,i) = X(:,i) ./ normDsmall';
end
