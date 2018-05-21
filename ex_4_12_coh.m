% Copyright (c) 2017 Bogdan Dumitrescu <bogdan.dumitrescu@acse.pub.ro>
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

% Coherence distribution for a UNTF
% Example 4.12. Generates Figure 4.3.

% BD 8.04.2017

m = 64;
n = 3*64;
csize = 20;

% generate UNTF
D = gen_rand_untf(m,n);

% compute scalar products and draw
G = abs(D'*D); % the scalar products
G = tril(G,-1); % G is symmetric and the diagonal is 1
v = sort(abs(G(:)));
v = v(n*(n+1)/2+1:end); % cut the zeros
plot(v)
hold on
plot(1:length(v), sqrt((n-m)/m/(n-1))*ones(1,length(v)), 'r')
hold off
grid
xlabel('\# product', 'interpreter', 'latex', 'FontSize', csize );
ylabel('Atom scalar products', 'interpreter', 'latex', 'FontSize', csize );

fprintf('Percentage of atom products larger than twice the Welch bound: %.4f\n',...
        100*sum(v > 2*sqrt((n-m)/m/(n-1))) / length(v))

%print -depsc fig_allcoh.eps
