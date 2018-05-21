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

function atom = new_atom(replatoms,Y,D,X,j)
%% Replace unused atom j
% INPUTS:
%   replatoms -- unused dictionary atoms replacement strategy
%       'zero': return the zero column
%       'random': return a random generated atom
%       'no': perform no replacement
%       'worst': replace with the worst represented signal
%
%   Y -- training signals set
%   D -- current dictionary
%   X -- sparse representations
%   j -- atom's column index in the dictionary D
%
% OUTPUTS:
%   atom -- replacement atom

    switch replatoms
        case 'zero'
            atom = 0;
        case 'random'
            atom = randn(size(D,1), 1);
            atom = atom / norm(atom);
        case 'no'
            atom = D(:,j);  % do nothing
        case 'worst'
            E = Y - D*X;
            errnorms = sqrt(sum(E.^2,1));
            [~,idx] = max(errnorms);
            atom = Y(:,idx)/norm(Y(:,idx));
    end
end