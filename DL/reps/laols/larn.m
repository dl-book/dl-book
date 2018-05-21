% Copyright (c) 2015 Bogdan Dumitrescu <bogdan.dumitrescu@acse.pub.ro>
% Copyright (c) 2015 Paul Irofti <paul@irofti.net>
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

function index = larn(A, y, k, K, indve)  
%
% Look-Ahead Residual Norm Atom Selection
% Input:
%  A      - matrix of basis vectors (assumed to be normed)
%  y      - vector to be matched
%  k      - current sparsity, this is sufficient because LAOLS
%           puts the atoms in the current support upfront in A
%  K      - target sparsity
%  indve  - index of look ahead candidates
% Output:
%  index  - index of the highest amplitude of the orthogonal projections

% Based on LAOLS implementation by BD 22.01.2015

    if k == 0
        indv = [];
    else
        indv = 1:k;
    end
    L = length(indve);

    nres = zeros(1,L);
    for k = 1 : L  % evaluate best candidates
        ind = [indv indve(k)];   % add candidate to support
        x = A(:,ind) \ y;     % LS solution for this support
        r = y - A(:,ind)*x;
        for j = k+1:K           % append to full support with OMP
            p = abs(A'*r);      % column projections on residual
            p(ind) = 0;
            [p,im] = max(p);
            ind = [ind im];
            x = A(:,ind) \ y;
            r = y - A(:,ind)*x;
        end
        nres(k) = norm(r);
    end
    [p, i] = min(nres); % take column with best residual for full OMP
    index = indve(i);
end