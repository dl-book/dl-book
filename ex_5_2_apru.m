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

% APrU (version OS) vs AK-SVD
% Generates values similar to those from Table 5.1,
% but the values of s and lambda must be changed manually

%-------------------------------------------------------------------------
m = 64;            % problem dimension
n = 128;
N = 8192;          % number of training signals
s = 10;            % sparsity level for AK-SVD
lambda = 0.00024;  % trade-off parameter for 1-norm penalty
iters = 50;        % number of DL iterations

update = 'aksvd';    % dictionary update method
replatom = 'random'; % atom replacement strategy

type = 'sliding';        % patches can be superposed
imdir = [pwd '\img\'];   % image directory
images = {'lena.png', 'barbara.png', 'peppers.png'};   % images
%%-------------------------------------------------------------------------

% Generate data (ONLY ONCE if multiple values of lam are used!!!)
[D0, Y] = gen_img_data(m, n, N, type, imdir, images);

% AK-SVD
params = {};
[D1, X1, errs1] = DL(Y, D0, s, iters, str2func(update), params, 'replatom', replatom);
fprintf('AK-SVD RMSE: %f\n', min(errs1))

% APrU
[D2, X2, errs2] = apru(Y, D0, lambda, iters, s);
fprintf('APrU RMSE:   %f\n', min(errs2))
fprintf('APrU average sparsity level %.2f\n', sum(sum(X2~=0))/N)

