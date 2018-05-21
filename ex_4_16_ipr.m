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

% Example 4.16: Plot atom inner products for DL bounding coherence with IPR
% Generates Figure 4.5
%-------------------------------------------------------------------------
m = 64;         % problem dimension
n = 192;        % number of atoms in the dictionary
N = 1024;       % number of training signals
s = 4;          % sparsity constraint
iters = 100;    % DL iterations
snr = 20;       % signal-to-noise ratio in dB

update = 'aksvd';    % dictionary update method
replatom = 'random';     % atom replacement strategy

imgdata = 0;                % test with images?
type = 'sliding';           % allow overlapping patches
imdir = [pwd '\img\'];      % image directory
images = {'lena.png'};   % images
%%-------------------------------------------------------------------------
%addpath(genpath('DL'))

% Generate data
if imgdata == 1
    [D, Y] = gen_img_data(m, n, N, type, imdir, images);
else
    [D, ~, Y] = gen_synth_data(m,n,N,s,snr);
end

% Dictionary Learning
params = {};
Dcoh = {};
errs = {};

% IPR with AK-SVD
[Dcoh{1}, X, errs{1}, criteria] = ...
      DL(Y, D, s, iters, str2func(update), params, ...
         'replatom', replatom, 'postopts', 'ipr', 'postoptsargs', {'mutual_coh_ipr', 0.2, 'nit_ipr', 5});

% IPR with AK-SVDc, gamma = 1;
update = 'aksvd_coh';
erropts = {'coherr', 'coh', 1};

[Dcoh{2}, X, errs{2}, criteria] = ...
      DL(Y, D, s, iters, str2func(update), params, ...
         'replatom', replatom, 'erropts', erropts, 'postopts', 'ipr', 'postoptsargs', {'mutual_coh_ipr', 0.2, 'nit_ipr', 5});

% Plot inner atom products in absolute value
linewid = 2;
colors = {'b','y'};
csize = 20;
figure(1)
for i = 1:2
  G = abs(Dcoh{i}'*Dcoh{i}); % the scalar products for learned dictionary
  G = tril(G,-1); % G is symmetric and the diagonal is 1
  v = sort(abs(G(:)));
  v = v(n*(n+1)/2+1:end); % cut the zeros
  fprintf('mutual coherence %.4f\n', max(v))  % mutual coherence
  fprintf('RMS coherence %.4f\n', norm(v)/sqrt(length(v)))  % RMS coherence
  plot(v, colors{i}, 'LineWidth', linewid)
  hold on
end
plot(1:length(v), sqrt((n-m)/m/(n-1))*ones(1,length(v)), 'r')
hold off
grid
xlabel('\# product', 'interpreter', 'latex', 'FontSize', csize );
ylabel('Atom scalar products', 'interpreter', 'latex', 'FontSize', csize );
h = legend('IPR AK-SVD', 'IPR AK-SVDc $\gamma=1$');
set(h, 'interpreter', 'latex', 'FontSize', csize, 'Location', 'NorthWest' );

%print -depsc fig_coh_ipr.eps

% Plot errors
figure(2)
for i = 1 : 2
  plot(errs{i}, colors{i}, 'LineWidth', linewid);
  hold on
end
hold off
grid
xlabel('iteration', 'interpreter', 'latex', 'FontSize', csize );
ylabel('RMSE', 'interpreter', 'latex', 'FontSize', csize );
h = legend('IPR AK-SVD', 'IPR AK-SVDc $\gamma=1$');
set(h, 'interpreter', 'latex', 'FontSize', csize );

%print -depsc fig_err_ipr.eps


