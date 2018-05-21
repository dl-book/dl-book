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

% Coherence distribution for dictionaries trained with AK-SVDc (Alg.4.2)
% Example 4.13. Generates Figure 4.4.

% BD 13.04.2017

%%-------------------------------------------------------------------------
m = 64;         % problem dimension
n = 192;        % number of atoms in the dictionary
N = 1024;       % number of training signals
s = 4;          % sparsity constraint
gamma_v = [0.1 1 10];    % coherence trade-off parameter
iters = 100;    % number of DL iterations
snr = 20;       % signal-to-noise ratio in dB

update = 'aksvd_coh';    % dictionary update method
replatom = 'random';     % atom replacement strategy

%%-------------------------------------------------------------------------

% Generate data
[D, ~, Y] = gen_synth_data(m,n,N,s,snr);

% Dictionary Learning
Dcoh = {};
errs = {};
criteria = {};
for i = 1 : length(gamma_v)
  params = {'coh', gamma_v(i)};
  erropts = {'coherr', 'coh', gamma_v(i)};
  
  [Dcoh{i}, X, errs{i}, criteria{i}] = ...
      DL(Y, D, s, iters, str2func(update), params, ...
         'replatom', replatom, 'erropts', erropts);
end

% Plot inner atom products in absolute value
linewid = 2;
colors = {'g','y','k'};
csize = 20;
figure(1)
for i = 1 : length(gamma_v)
  G = abs(Dcoh{i}'*Dcoh{i}); % the scalar products for learned dictionary
  G = tril(G,-1); % G is symmetric and the diagonal is 1
  v = sort(abs(G(:)));
  v = v(n*(n+1)/2+1:end); % cut the zeros
  fprintf('gamma = %5.2f\n', gamma_v(i))
  fprintf('  mutual coherence: %.4f\n', max(v))
  fprintf('  RMS coherence:    %.4f\n', norm(v)/sqrt(length(v)) )
  plot(v, colors{i}, 'LineWidth', linewid)
  hold on
end
plot(1:length(v), sqrt((n-m)/m/(n-1))*ones(1,length(v)), 'r')
hold off
grid
xlabel('\# product', 'interpreter', 'latex', 'FontSize', csize );
ylabel('Atom scalar products', 'interpreter', 'latex', 'FontSize', csize );
h=legend('$\gamma=0.1$', '$\gamma=1$', '$\gamma=10$');  % lazy !!!
set(h, 'interpreter', 'latex', 'FontSize', csize, 'Location', 'NorthWest' );

%print -depsc fig_coh_gamma.eps

% Plot errors
figure(2)
for i = 1 : length(gamma_v)
  plot(errs{i}, colors{i}, 'LineWidth', linewid);
  hold on
end
hold off
grid
xlabel('iteration', 'interpreter', 'latex', 'FontSize', csize );
ylabel('RMSE', 'interpreter', 'latex', 'FontSize', csize );
h=legend('$\gamma=0.1$', '$\gamma=1$', '$\gamma=10$');  % lazy !!!
set(h, 'interpreter', 'latex', 'FontSize', csize );

%print -depsc fig_err_gamma.eps


