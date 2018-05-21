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

function [X,shared] = omp(Y,D,s,shared,varargin)
%% Orthogonal Matching Pursuit algorithm
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   s -- sparsity constraint
%
% OUTPUTS:
%   X -- sparse representations

    if isempty(varargin)
        ompparams = {'checkdict', 'off'};
        do_sparse = true;
    else
        if strcmp(varargin{1}, 'error')
            do_sparse = false;
            error = varargin{2};
            if length(varargin) > 2
                ompparams = {varargin{3:end}};
            else
                ompparams = {'checkdict', 'off'};
            end
        else
            do_sparse = true;
            ompparams = varargin;
        end
    end
    if do_sparse
        X = omp_sparse(D'*Y, D'*D, s, ompparams{:});
    else
        X = omp_err(D'*Y, sum(Y.*Y), D'*D, error, ompparams{:});
    end
end
