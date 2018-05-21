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

function [Y, D, X] = denoise_2D(Y, D, varargin)
%% Image denoising with 2D-OMP
% INPUTS:
%   Y -- 0-mean vectorized patches
%   D -- 2-cell containing left handside and right handside dictionaries
%   sigma -- noise standard deviation
% OUTPUTS:
%   Y -- vectorized patches (mean included)
%   X -- sparse representations
    D1 = D{1}; D2 = D{2}; D = [];
    
    p = inputParser;
    p.KeepUnmatched=true;
    p.PartialMatching=false;
    p.addParameter('sigma', 20);
    p.parse(varargin{:});
    sigma = p.Results.sigma;
   
    N = size(Y, 2);
	[p1, n1] = size(D1);
    [n2, p2] = size(D2);
    X = zeros(n1,n2,N);
    %Y2D = vec2patch(Y, p1, p2);
    Y2D = reshape(Y, p1, p2, N);
    Y = zeros(p1,p2,N);
    
    max_s = p1*p2/2;	% maximum density is half the patch
    gain = 1.15;        % default noise gain
    error = sqrt(p1*p2)*gain*sigma;
    %error = gain*sigma;
    
    for k = 1:N
        X(:,:,k) = pair_omp_err(Y2D(:,:,k),D1,D2,max_s,error);
        Y(:,:,k) = D1*X(:,:,k)*D2;
    end
    Y = reshape(Y, p1*p2, N);
end