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

function [D,X,shared] = nsgk_loop(Y,D,X,iter,replatoms,customfunc)
%% NSGK dictionary update loop
%
% This is an NSGK adaptation of atom_up. It was needed for keeping track
% of past iteration dictionary and representations.
%
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   X -- sparse representations
%   replatoms -- unused dictionary atoms replacement strategy
%   customfunc -- algorithm specific atom update routine
%                 The function has to respect the API: 
%                   [d,x] = func(Y,D,X,d,x)
%
% OUTPUTS:
%   D -- updated dictionary
%   X -- updated representations
    persistent D0;
    persistent X0;
    
	if iter == 1
        X0 = X;
        D0 = D;
	end
    Z = Y + D0*X0 - D0*X;
    for j = 1:size(D,2)  
        [~, data_indices, ~] = find(X0(j,:));
        if (isempty(data_indices))
            D(:,j) = new_atom(replatoms,Y,D,X0,j);
            continue;
        end
        d = D(:,j);
        x = X0(j,data_indices);
        %x = X0(j,:);
        
        F = Z(:,data_indices) - D*X0(:,data_indices) + d*x;
        %F = Z - D*X0 + d*x;
        
        [D(:,j),X0(j,data_indices)] = customfunc(F,D,d,x);
    end
	X0 = X;
    D0 = D;
end