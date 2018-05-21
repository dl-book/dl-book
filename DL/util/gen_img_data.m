% Copyright (c) 2015 Paul Irofti <paul@irofti.net>
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

function [D, Y] = gen_img_data(p, n, m, type, imdir, images)
%% Generates dictionary learning data from images
% p -- patch size when vectorized (i.e. 8x8 patch means p=64)
% n -- number of atoms in the dictionary
% m -- total number of patches
% type -- patch construction (default: 'sliding')
% imdir -- directory path containing the images
% images -- cell array of training figures

% D -- random normalized dictionary
% Y -- processed image patches into column vectors
    if nargin < 4
        type = 'sliding';
    end
    if nargin < 5
        imdir='C:\cygwin\home\bulibuta\wrk\img\';
        images = {'4.1.05.tiff ','4.1.06.tiff','4.2.01.tiff', ...
            '4.2.02.tiff', '4.2.05.tiff','4.2.06.tiff', ...
            '5.1.09.tiff','5.2.08.tiff', '5.2.09.tiff','5.2.10.tiff'};
    end
    
    Y = mkimgsigs(imdir,images,p,m,type);
    if isempty(m) == 0
        Y = Y(:,randperm(size(Y,2),m));
    end

    % Generate dictionary
    D = normc(randn(p,n));
end