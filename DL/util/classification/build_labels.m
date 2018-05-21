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

function H = build_labels(parts, samples)
%% Build node labels for LC-KSVD classification
% INPUTS:
%   parts -- node partitions
%   samples -- samples per node (equals the number of emitters)
%
% OUTPUTS:
%   H -- class labels of input signals
%--------------------------------------------------------------------------
    nodes = max(cell2mat(parts));
    classes = length(parts);
    H = zeros(classes, nodes*samples);
    cmarker = ones(1,samples);
    
    for c = 1:classes
        for n = parts{c}
            H(c, (n-1)*samples + 1:n*samples) = cmarker;
        end
    end
end