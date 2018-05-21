% Copyright (c) 2018 Paul Irofti <paul@irofti.net>
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

function [success, X] = classification(Y, H, D, W, s)
%% Perform classification and return successfulness in percentage
% INPUTS:
%   Y -- test signals set
%   H -- ground truth label matrix
%   D -- dictionary
%   W -- classifier matrix
%   s -- sparsity level
%
% OUTPUTS:
%   success -- succesful classification percentage
%   X       -- sparse representations

X = omp(Y,D,s);          % sparse representation
[~,estimate] = max(W*X); % classification
[~, truth] = max(H);     % ground truth

% Succesful classification percentage
success = sum(estimate == truth)*100/length(truth);