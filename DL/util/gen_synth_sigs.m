% Copyright (c) 2017 Paul Irofti <paul@irofti.net>
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

function Yawgn = gen_synth_sigs(D0,p,n,m,s,SNRdB,thresh)

    if nargin < 7
        thresh = 0;
    end
    % Traning signals generated as a linear combination of s atoms
    % at random locations with i.i.d. coefficients
    Y0 = zeros(p,m);
    for i = 1:m
        pos = randperm(n,s);
        coeff = randn(s,1);
        if nargin >= 7
            coeff = coeff + thresh * sign(coeff);
        end
        for s0 = 1:s           
            Y0(:,i) = Y0(:,i) + coeff(s0)*D0(:,pos(s0));
        end
    end
    
    % White Gaussian noise is added to the resulting signals.
    if SNRdB == Inf
        Yawgn = Y0;
    else
        Yawgn = awgn(Y0,SNRdB,'measured','db');
    end
end