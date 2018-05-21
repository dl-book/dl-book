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

function [D,X] = patom_up(Y,D,X,replatoms,pa,customfunc)
%% Parallel dictionary update loop
% INPUTS:
%   Y -- training signals set
%   D -- current dictionary
%   X -- sparse representations
%   pa -- number of atoms to update in parallel
%   replatoms -- unused dictionary atoms replacement strategy
%   customfunc -- algorithm specific atom update routine
%                 The function has to respect the API: 
%                   [d,x] = func(Y,D,X,d,x)
%
% OUTPUTS:
%   D -- updated dictionary
%   X -- updated representations
    for j = 1:size(D,2)/pa
        Ei = Y - D*X;
        for p = (j-1)*pa + 1:j*pa
            [~, data_indices, x] = find(X(p,:));
            if (isempty(data_indices))
                D(:,p) = new_atom(replatoms,Y,D,X,p);
                continue;
            end
            
            d = D(:,p);         
            F = Ei(:,data_indices) + d*x;

            [D(:,p),X(p,data_indices)] = customfunc(F,D,d,x);
        end
    end
end