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

function [D,X,shared] = sgk_reg(Y,D,X,iter,replatoms,shared,varargin)
%% SGK algorithm with regularization
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   X -- sparse representations
%   iter -- current DL iteration
%
% PARAMETERS:
%   reg -- regularization factor (default: 0.01)
%   vanish -- regularization vanishing factor (default: 0.95)
%
% OUTPUTS:
%   D -- updated dictionary

    persistent mu;
    persistent vanish;
    
    if iter == 1
        p = inputParser();
        p.KeepUnmatched=true;
        p.addParameter('reg', 0.01);
        p.addParameter('vanish', 0.95);
        p.parse(varargin{:});
        mu = p.Results.reg;
        vanish = p.Results.vanish;
    else
        mu = mu * vanish;
    end

    % Use sparse pattern from OMP to compute new sparse representations
    X = ompreg(Y,D,X,mu);
    
    [D,X] = atom_up(Y,D,X,replatoms,@(Y,D,X,d,x) reg_up(Y,D,X,x));
end

function [d,x] = reg_up(Y,D,X,x)
%% SGK with regularization atom update
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   X -- sparse representations
%   x -- sparse representations row using the current atom
%
% OUTPUTS:
%   d -- updated atom
    d = (Y*x' - D*(X*x'))/(x*x');
    d = d/norm(d); % comment this line to obtain the original SGK algorithm
end
