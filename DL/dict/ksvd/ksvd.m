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

function [D,X,shared] = ksvd(Y,D,X,~,replatoms,shared,varargin)
%% K-SVD algorithm
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   X -- sparse representations
%
% OUTPUTS:
%   D -- updated dictionary
%   X -- updated representations
    [D,X] = atom_up(Y,D,X,replatoms,@(Y,D,X,d,x) ksvd_up(Y,D,X));
end

function [d,x] = ksvd_up(Y,D,X)
%% K-SVD atom update
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   X -- sparse representations
%
% OUTPUTS:
%   d -- updated atom
%   x -- updated representations corresponding to the current atom
    [d,s,x] = svds(Y - D*X, 1);
    x = s*x;
end