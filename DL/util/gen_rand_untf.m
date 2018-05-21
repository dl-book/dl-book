% Copyright (c) 2017-2018 Bogdan Dumitrescu <bogdan.dumitrescu@acse.pub.ro>
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

function D = gen_rand_untf(m,n)

% Generate random unit norm tight frame of size mxn as described in Algorithm 4.1

% BD 10.4.2017

% generate random tight frame
D = randn(m,n);
[Q,~] = qr(D',0);
D = Q' * sqrt(n/m);  % rows have correct norm

% force atoms to unit norm
atom_norms = sum(D.*D);
for i = 1 : n-1
  if abs(atom_norms(i) - 1) > 1e-12  % do nothing if it happens that the norm is 1
    s1 = sign(atom_norms(i)-1);
    j = i+1;
    while sign(atom_norms(j)-1) == s1  % find atom with norm on the other side of 1
      j = j+1;
    end
    % compute tangent of rotation angle
    an1 = atom_norms(i);
    an2 = atom_norms(j);
    cp = D(:,i)'*D(:,j);
    t = (cp + sign(cp)*sqrt(cp*cp - (an1-1)*(an2-1))) / (an2-1);
    % compute rotation
    c = 1 / sqrt(1+t*t);
    s = c*t;
    % new atoms and updated norm
    D(:,[i j]) = D(:,[i j]) * [c s; -s c];
    atom_norms(j) = an1+an2-1;
  end
end
