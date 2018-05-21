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

function I = sigs2img(Yp, Ywidth, Ysigs, zcols, Ymean)
%% Restore a processed image that was split by mksigs
% Yp -- processed dataset
% Ywidth -- the width of the original image
% Ysigs -- total number of columns resulted from vectorizing the patches
% zcols -- the number of eliminated zero columns
% Ymean -- the mean substracted from the patches at processing

% I -- the processed image

%% Restore the eliminated zero columns
nzcols = setdiff(1:Ysigs, zcols);
Y = zeros(size(Yp,1), Ysigs);
Y(:,nzcols) = Yp;

%% Denormalize and add back the mean
%Y = round(Y.*255);
%Y = Y + repmat(Ymean,size(Y,1),size(Y,2));

Y = im2uint8(Y - min(Y(:)));

%% Restore the image patch by patch
col = 1;
row = 1;
for i=1:Ysigs
	patch = reshape(Y(:,i),8,8);
	I(row:row + 7, col:col + 7) = patch;
	if col + 8 > Ywidth
		col = 1;
		row = row + 8;
	else
		col = col + 8;
	end
end
%imshow(I, [0 255]);
end