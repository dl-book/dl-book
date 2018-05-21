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

function [D,X,shared] = psgk_coh(Y,D,X,~,replatoms,shared,varargin)
%% P-SGK algorithm with coherence reduction
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   X -- sparse representations
%   iter -- current DL iteration
%
% PARAMETERS:
%   pa -- number of atoms to update in parallel (default: all)
%   coh -- coherence reduction factor (default: 6)
%
% OUTPUTS:
%   D -- updated dictionary
    p = inputParser();
    p.KeepUnmatched=true;
    p.addParameter('pa', size(D,2));
    p.addParameter('coh', 6);
    p.parse(varargin{:});
    gamma = p.Results.coh;
    pa = p.Results.pa;
    
    [D,X] = patom_up(Y,D,X,replatoms,pa,@(F,D,d,x) coh_up(F,D,d,x,gamma));
end

function [d,x] = coh_up(F,D,d,x,gamma)
%% P-SGK algorithm with coherence reduction atom update
% INPUTS:
%   F -- approximation error w/o the current atom
%   D -- current dictionary
%   d -- current atom
%   x -- sparse representations row using the current atom
%   gamma -- coherence reduction factor
%
% OUTPUTS:
%   d -- updated atom
    % Representation with the old atom
    x0 = F'*d;
    % Update atom
    d = F*x0 - 2*gamma*D*(D'*d) + 2*sqrt(size(D,1))*gamma*d;
    d = d/norm(x0)^2;
    d = d/norm(d); % comment this line to obtain the original SGK algorithm
end
