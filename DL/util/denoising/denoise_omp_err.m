% Copyright (c) 2017 Paul Irofti <paul@irofti.net>
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

function [Y, D, X] = denoise_omp(Y, D, varargin)
%% Image denoising with sparsity driven OMP
% INPUTS:
%   Y -- 0-mean vectorized patches
%   D -- dictionary
%   s -- sparsity
% OUTPUTS:
%   X -- sparse representations
    p = inputParser;
    p.KeepUnmatched=true;
    p.addParameter('s', 8);
    p.parse(varargin{:});
    s = p.Results.s;

    X = omp(Y,D,s);
    Y = [];
end