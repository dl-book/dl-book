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

function [X, shared] = laols(Y,D,s,shared,varargin)
%% Look-Ahead Orthogonal Least Squares
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   s -- sparsity constraint
% OPTIONAL:
%   ahead -- number of look-ahead candidates stored in varargin(default: s)
%
% OUTPUTS:
%   X -- sparse representations
    X = zeros(size(D,2),size(Y,2));

    if isempty(varargin)
        ahead = s;
    else
        ahead = varargin{:};
    end

    for j=1:size(Y,2)
        [indv,x] = laols1(D, Y(:,j), s, ahead);
        x_la = zeros(size(D,2),1);
        x_la(indv(1:s)) = x;
        X(:,j) = x_la;
    end
end
