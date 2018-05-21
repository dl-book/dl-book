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

function [Y , Irows, Icols, sigs, zcols, Ymean, Ynorm] = ...
    mkimgsigs(imdir,images,p,m, type, center, normalize)
%% Creates random image patches
% imdir -- directory path containing the images
% images -- cell array of training figures
% p -- patch size
% m -- total number of patches
% type -- 'sliding' or 'distinct' patches
% center -- center signals?
% normalize -- normalize signals?

% Y -- processed image patches into column vectors
% Iwidth -- original image width
% sigs -- number of signals including zero columns
% zcols -- indices of the eliminated null columns
% Ymean -- substracted mean from each column

    %[Y, ~, Iwidth, Ymean] = readImages(imdir,images);
    Y = [];
    if nargin < 5
        type = 'sliding';
    end
    if nargin < 6
        center = 1;
    end
    if nargin < 7
        normalize = 1;
    end
    for i=1:length(images)
        I = double(imread([imdir,char(images(i))]));
        I = I(:, :, 1);
        Y = [Y im2col(I, [sqrt(p) sqrt(p)], type)];
    end
    Icols = size(I,2); % XXX: useful only when a single image is used
    Irows = size(I,1);
    Ymean = mean(Y);
    if center
        Y = Y - repmat(Ymean, p, 1);    % 0-mean
    end
    if normalize
        Ynorm = norm(Y);
        Y = Y / norm(Y);                % l2 normalized
    end
    %Y = Y./255;
	sigs = size(Y,2);

	% Find the non-zero columns
	[~,cols] = find(Y);             % there are the non-empty ones but w/ dups
	zcols = setdiff(1:sigs,cols);   % these are the empty ones
	cols = setdiff(1:sigs,zcols);   % these are the non-empty ones w/o dups
	Y = Y(:,cols);
	
	if nargin >= 4 && m > 0
		sigs = size(Y,2);
		patches = randperm(sigs, m);
		Y = Y(1:p, patches);
	end
end