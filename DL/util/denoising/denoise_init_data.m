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

function [I, Inoisy, Y, Ymean] = denoise_init_data(img, sigma, p1, p2)
%% Add white gaussian noise to image
% INPUTS:
%   img -- image path
%   sigma -- noise standard deviation
%   p1,p2 -- patch size
% OUTPUTS:
%   Inoisy -- noisy image
%   Y -- 0-mean vectorized patches from Inoisy
%   Ymean -- substracted mean
    I = double(imread(img));
    I = I(:, :, 1);
    Inoisy = I + sigma*randn(size(I));
    Y = im2col(Inoisy, [p1 p2], 'sliding');
    Ymean = mean(Y);
    Y = Y - repmat(Ymean,size(Y,1),1);
end