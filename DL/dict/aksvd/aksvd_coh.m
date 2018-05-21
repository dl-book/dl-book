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

function [D,X,shared] = aksvd_coh(Y,D,X,~,replatoms,shared,varargin)
%% AK-SVD algorithm with coherence reduction
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   X -- sparse representations
%
% PARAMETERS:
%   coh -- coherence reduction factor (default: 6)
%
% OUTPUTS:
%   D -- updated dictionary
%   X -- updated representations
    p = inputParser();
    p.KeepUnmatched=true;
    p.addParameter('coh', 6);
    p.parse(varargin{:});
    gamma = p.Results.coh;
    
    [D,X] = atom_up(Y,D,X,replatoms,@(Y,D,X,d,x) coh_up(Y,D,X,d,gamma));
end

function [d,x] = coh_up(Y,D,X,d,gamma)
%% AK-SVD with coherence reduction atom update
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   X -- sparse representations
%   d -- current atom
%   gamma -- coherence reduction factor
%
% OUTPUTS:
%   d -- updated atom
%   x -- updated representations corresponding to the current atom

    % first value of the representation
    x = Y'*d - X'*(D'*d);
    x = x';
    p = size(D,1);

    % power method iteration
    d = Y * x' - D*(X*x') - 2*gamma*D*(D'*d) + 2*sqrt(p)*gamma*d;
    d = d/norm(d);
    x = Y'*d - X'*(D'*d);
    x = x';
end