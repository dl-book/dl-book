% Copyright (c) 2016 Paul Irofti <paul@irofti.net>
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

function [D,X,shared] = aksvd_var2s(Y,D,X,iter,replatoms,shared,varargin)
%% Approximate double sparsity K-SVD algorithm with variable dictionary
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   X -- sparse representations
% PARAMETERS:
%   Q -- dense dictionary (default: eye())
%   q -- atom sparsity (default: 4)
%
% OUTPUTS:
%   Q -- updated dictionary
%   D -- updated dictionary
%   X -- updated representations
    Q = shared{1};
    persistent as;

    if iter == 1
        p = inputParser();
        p.KeepUnmatched=true;
        p.addParameter('as', 4, @(x) x <= size(D,1));
        p.parse(varargin{:});
        as = p.Results.as;
    end

    % Update the dense dictionary
    [Q,D,X] = var2s_up(Y,Q,D,X,replatoms);
    shared{1} = Q;
    % Update sparse dictionary and representations
    [D,X] = atom_up(Y,D,X,replatoms,@(Y,D,X,d,x) sparse_up(Y,Q,D,X,x,as));
end

function [q,z] = dense_up(Y,Q,Z,z)
%% Approximate K-SVD dense dictionary atom update
% INPUTS:
%   Y -- training signals set
%   Q -- dense dictionary
%   Z -- D*X factorization
%   z -- representations row using the current atom
%
% OUTPUTS:
%   q -- updated atom
    F = Y - Q*Z;
    q = F*z'/norm(F*z');
end

function [Q,D,X] = var2s_up(Y,Q,D,X,replatoms)
%% Sequential dictionary update loop
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   X -- sparse representations
%   replatoms -- unused dictionary atoms replacement strategy
%   customfunc -- algorithm specific atom update routine
%                 The function has to respect the API: 
%                   [d,x] = func(Y,D,X,d,x)
%
% OUTPUTS:
%   D -- updated dictionary
%   X -- updated representations
    for j = 1:size(Q,2)
        Z = D*X;
        [~, data_indices, z] = find(Z(j,:));
        [~, dindices, ~] = find(D(j,:));
        [~, xindices, ~] = find(X(dindices,:));
        xindices = unique(xindices);
        Q(:,j) = 0;
        if (isempty(data_indices))
            Q(:,j) = new_atom(replatoms,Y,Q,D,j);
            continue;
        end

        % Update dense atom j
        smallY = Y(:,xindices);
        smallZ = Z(:,xindices);
        smallX = X(dindices,xindices);
        F = smallY - Q*smallZ;
        q = F*z'/norm(F*z');
        Q(:,j) = q;

        % Use the existing support, found via OMP, and update the 
        % sparse atom coefficients via least-squares.
        d = smallX'\(F'*q);
        D(j,dindices) = d;

        % Update representations that are using the updated sparse atoms
        for k = dindices
            d = D(:,k);
            D(:,k) = 0;
            [~, xindices, ~] = find(X(k,:));
            smallY = Y(:,xindices);
            smallX = X(:,xindices);
            F = smallY - Q*D*smallX;
            X(k,xindices) = F'*Q*d;
            D(:,k) = d;
        end
    end
end

function [d,x] = sparse_up(Y,Q,D,X,x,q)
%% Approximate K-SVD sparse dictionary atom update
% INPUTS:
%   Y -- training signals set
%   Q -- fixed dictionary
%   D -- sparse dictionary
%   X -- sparse representations
%   x -- sparse representations row using the current atom
%   q -- atom sparsity
%
% OUTPUTS:
%   d -- updated atom
%   x -- updated representations corresponding to the current atom
    F = Y - Q*D*X;
    x = x/norm(x);
    z = F*x';
    d = omp(z,Q,q);
    d = d/norm(Q*d);
    x = F'*Q*d;
end
