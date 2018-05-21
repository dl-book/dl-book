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

function [D,X,shared] = MOD(Y, D, X, ~, ~, shared, varargin)
%% Method of Optimal Directions (MOD) algorithm
%
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   X -- sparse representations
%
% OUTPUTS:
%   D -- updated dictionary

    % Use sparse pattern from OMP to compute new sparse representations
    X = ompreg(Y,D,X,0);
	
	% I use SVD here to avoid the singularity issues that I encountered 
	% when plainly solving the system
	%[U, S, V] = svd(X',0);
    [U, S, V] = svds(X',min(size(X')));
	for k = 1:size(D,2)
		if S(k,k) < 1e-8
			S(k,k) = 0;
		else
			S(k,k) = 1/S(k,k);
		end
	end
	D = V*S*U'*Y';
	D = D';
    D = normc(D);
end
