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

% Example 5.9: compare dictionary selection with AK-SVD

% BD 22.05.2017

m = 16;         % problem dimension
n = 40;         % number of atoms in the dictionary
na = 200;       % number of candidate atoms
N = 1000;       % number of training signals
s = 4;          % sparsity constraint
iters = 100;    % DL iterations for AK-SVD
snr = 20;       % signal-to-noise ratio in dB

imgdata = 0;                % test with images?
type = 'sliding';           % allow superposing patches
imdir = [pwd '\img\'];      % image directory
images = {'lena.png', 'barbara.png', 'peppers.png'};   % images
%%-------------------------------------------------------------------------
% Generate data
if imgdata == 1
  [~, Y] = gen_img_data(m, n, N, type, imdir, images);
else
  [~, ~, Y] = gen_synth_data(m,n,N,s,snr);
end

% generate candidate atoms
if 0   % random
  A = normc(randn(m,na));
else   % from signals
  A = normc(Y(:,randperm(N,na)));
end
X = omp(Y,A,s);
[U,~,V] = svd(A*X*Y');
A = V*U'*A;

% dictionary selection
[D, X] = dl_select(Y, A, n, s);

fprintf('Dictionary selection RMSE: %f\n', norm(Y-D*X,'fro')/sqrt(m*N) )
fprintf('Dictionary selection coh.: %f\n', max(max(D'*D-eye(n))) )

% AK-SVD
update = 'aksvd';        % dictionary update method
replatom = 'random';     % atom replacement strategy
params = {};
D0 = normc(randn(m,n));

[Dak, Xak, errs] = DL(Y, D0, s, iters, str2func(update), params, 'replatom', replatom);
fprintf('AK-SVD RMSE: %f\n', norm(Y-Dak*Xak,'fro')/sqrt(m*N) )
fprintf('AK-SVD coh.: %f\n', max(max(Dak'*Dak-eye(n))) )

% AK-SVD with IPR
[Dcoh, Xcoh, errs] = DL(Y, D0, s, iters, str2func(update), params, ...
         'replatom', replatom, 'postopts', 'ipr', 'postoptsargs', {'mutual_coh_ipr', 0.4, 'nit_ipr', 5});
fprintf('AK-SVD+IPR RMSE: %f\n', norm(Y-Dcoh*Xcoh,'fro')/sqrt(m*N))
fprintf('AK-SVD+IPR coh.: %f\n', max(max(Dcoh'*Dcoh-eye(n))) )

