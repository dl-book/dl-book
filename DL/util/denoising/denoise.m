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

function [Iclean, Yclean, X] = denoise(dfunc, Inoisy, D, varargin)
% INPUTS:
%   dfunc -- denoising function, prototype: X = dfunc(Y, D, sigma)
%   Inoisy -- noisy image
%   D1,D2 -- left handside and right handside dictionaries
% OPTIONAL:
%   dargs -- denoise function arguments
%   p1,p2 -- patch size (defaults to 8)
% OUTPUTS:
%   I -- denoised image
%   Y -- vectorized patches (mean included)
%   X -- sparse representations
    p = inputParser;
    
    p.addOptional('dargs', {});
    p.addOptional('p1', 8);
    p.addOptional('p2', 8);
    
    p.parse(varargin{:});
    res = p.Results;
    p1 = res.p1;
    p2 = res.p2;

    [Ynoisy,Ymean] = denoise_proc_input(Inoisy,p1,p2);
    %{
    % Call original omp denoising method from Elad
    sigma = 20;
    gain = 1.15;
    maxval = 255;
    params.x = Inoisy;
    params.blocksize = [p1 p2];
    params.maxval = maxval;
    params.memusage = 'high';
    params.dict = D;
    params.Edata = sqrt(p1*p2) * sigma * gain;   % target error for omp
    params.codemode = 'error';
    params.sigma = sigma;
    params.noisemode = 'sigma';
    msgdelta = -1;
    Iclean = ompdenoise2(params,msgdelta);
    Yclean = [];
    X = [];
    %}
    [Yclean, D, X] = dfunc(Ynoisy, D, res.dargs{:});
    [Iclean, Yclean] = ...
        denoise_proc_output(Inoisy, Yclean, D, X, Ymean, p1, p2);
    
    % Average with noisy image
    %{
    sigma = 20;      % test_denoise first sigma
    maxval = 255;   % grayscale image
    lambda = maxval/(10*sigma);   % default value from Elad's ompdenoise2
    
    % Compute overlaping pixels in Iclean
    [m,n] = size(Inoisy);
    inds = reshape(1:m*n,[m n]);
    subs = im2col(inds, [res.p1 res.p2], 'sliding');
    Iover = accumarray(subs(:),1);
    Iover = reshape(Iover,m,n);

    %Iclean = Iclean + ...
    %    (lambda*(Inoisy - Iclean))./(Iover + lambda*ones(m,n));

    %}
end