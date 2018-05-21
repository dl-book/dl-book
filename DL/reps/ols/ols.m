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

function [X, shared] = ols(Y,D,s,shared,varargin)
%% Orthogonal Least Squares (or Optimized Orthogonal Matching Pursuit)
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   s -- sparsity constraint
%
% OUTPUTS:
%   X -- sparse representations
    X = zeros(size(D,2),size(Y,2));

    for j=1:size(Y,2)
        [indv,x] = ols1(D, Y(:,j), s);
        x_ols = zeros(size(D,2),1);
        x_ols(indv(1:s)) = x;
        X(:,j) = x_ols;
    end
end
