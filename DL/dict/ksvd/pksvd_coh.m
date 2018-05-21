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

function [D,X,shared] = pksvd_coh(Y,D,X,~,replatoms,shared,varargin)
%% PK-SVD algorithm with coherence reduction
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   X -- sparse representations
%
% PARAMETERS:
%   pa -- number of atoms to update in parallel (default: all)
%   coh -- coherence reduction factor (default: 6)
%
% OUTPUTS:
%   D -- updated dictionary
%   X -- updated representations
    p = inputParser();
    p.KeepUnmatched=true;
    p.addParameter('pa', size(D,2));
    p.addParameter('coh', 6);
    p.parse(varargin{:});
    gamma = p.Results.coh;
    pa = p.Results.pa;
    
    [D,X] = patom_up(Y,D,X,replatoms,pa,@(F,D,d,x) coh_up(F,D,gamma));
end

function [d,x] = coh_up(F,D,gamma)
%% PK-SVD algorithm with coherence reduction atom update
% INPUTS:
%   F -- approximation error w/o the current atom
%   D -- current dictionary
%   gamma -- coherence reduction factor
%
% OUTPUTS:
%   d -- updated atom
%   x -- updated representations corresponding to the current atom
	H = F*F' - 2*gamma*(D*D');
	[d,~] = eigs(H + 2*gamma*sqrt(size(D,1))*eye(size(H,1)), 1);
    x = F'*d;
end