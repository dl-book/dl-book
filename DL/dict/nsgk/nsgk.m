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

function [D,X,shared] = nsgk(Y,D,X,iter,replatoms,shared,varargin)
%% NSGK algorithm
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   X -- sparse representations
%   iter -- current DL iteration
%
% OUTPUTS:
%   D -- updated dictionary
    [D,X] = nsgk_loop(Y,D,X,iter,replatoms, ...
        @(F,D,d,x) nsgk_up(F,x));
end

function [d,x] = nsgk_up(F,x)
%% NSGK atom update
% INPUTS:
%   F -- approximation error w/o the current atom
%   x -- sparse representations row using the current atom
%
% OUTPUTS:
%   d -- updated atom
    d = F*x';
    d = d/norm(x)^2;
    d = d/norm(d); % comment this line to obtain an algorithm in SGK style
end
