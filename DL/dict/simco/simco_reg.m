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

function [D, X, shared] = simco_reg(Y,D,X,iter,~,shared,varargin)
%% Regularized SimCO
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
%   X -- updated representations

%  Based on software writen by:
%  Wenwu Wang <w.wang@surrey.ac.uk>
%  Department of Electronic Engineering 
%  University of Surrey
    persistent mu;
    persistent vanish;
    persistent IPara;
    persistent regstop;
    
    if iter == 1
        p = inputParser();
        p.KeepUnmatched=true;
        p.addParameter('reg', 0.01);
        p.addParameter('vanish', 0.95);
        p.addParameter('itN', 1);
        p.addParameter('regstop', Inf);
        p.parse(varargin{:});
        mu = p.Results.reg;
        vanish = p.Results.vanish;
        regstop = p.Results.regstop;
            
        IPara.mu = mu;
        IPara.I = 1:size(D,2);
        IPara.DebugFlag = 0;
        IPara.itN = p.Results.itN;
        IPara.gmin = 1e-5; % the minimum value of gradient
        IPara.Lmin = 1e-6; %t4-t1 should be larger than Lmin
        IPara.t4 = 1e-2; %the initial value of t4
        IPara.rNmax = 3; %the number of iterative refinement in Part B in DictLineSearch03.m
    else
        if regstop > iter
            mu = mu * vanish;
        else
            mu = 0;
        end
        IPara.mu = mu;
    end

    [D,X,~] = DictUpdate03(Y,D,X,IPara);
end

