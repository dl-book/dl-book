% Copyright (c) 2016-2018 Paul Irofti <paul@irofti.net>
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

function I = sliding2im(Y, Idim, pos, npatches, patchsz)
%% Rebuild image block from sliding im2col generated matrix Y
% INPUTS:
%   Y -- overlapping patches vectorized as matrix columns
%   Idim -- original image dimension (assumed height=width)
%   pos -- Y column from which to start building the square image block
%   npatches -- number of patches per block dimension
%   patchsz -- patch size (assumed square)
% OUTPUTS:
%   I -- resulting image block
I = [];
for j = 0:npatches-1
    Ij = [];
    for i = 0:npatches-1
        Ysub = Y(:,(pos + j*(Idim-(patchsz-1))*patchsz) + i*patchsz);
        Ii = reshape(Ysub,patchsz,patchsz);
        Ij = [Ij; Ii];
    end
    I = [I Ij];
end
end