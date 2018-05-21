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

function [I, Y] = denoise_proc_output(Inoisy, Y, D, X, Ymean, p1, p2)
%% Create all possible 0-mean vectorized patches from image I
% INPUTS:
%   I -- image
%   p1,p2 -- patch size
% OUTPUTS:
%   Y -- 0-mean vectorized patches
%   Ymean -- substracted mean
    [m,n] = size(Inoisy);
    % Allow the denoising function to pass Y to us if it knows better
    if isempty(Y)
        Y = D*X;
    end
    Y = Y + ones(p1*p2,1) * Ymean;
    I = avgcol2im(Y, m, n, p1, p2);
end