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
%   sigma -- noise standard deviation
% OUTPUTS:
%   X -- sparse representations
    p = inputParser;
    p.KeepUnmatched=true;
    p.PartialMatching=false;
    p.addParameter('sigma', 20);
    p.parse(varargin{:});
    sigma = p.Results.sigma;
    
    max_s = size(Y,1)/2;        % maximum density is half the patch
    gain = 1.15;                % default noise gain
    params = {'error', sqrt(size(Y,1))*gain*sigma, 'maxatoms', max_s};

    X = omp(Y,D,max_s,[],params{:});
    Y = [];
end