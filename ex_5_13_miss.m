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

% Example 5.13. AK-SVD with missing data

% BD 9.12.2017

%%-------------------------------------------------------------------------
m = 64;         % problem dimension
n = 128;        % number of atoms in the dictionary
N = 8000;       % 4000 and 8000 number of training signals
s = 6;          % sparsity constraint
iters = 50;     % DL iterations
missing_data_ratio = 0.3;   % how many data are missing, on average

type = 'distinct';          % distinct patches
imdir = [pwd '\img\'];      % image directory
images = {'barbara.png', 'boat.png', 'house.png', 'lena.png', 'peppers.png'};
%%-------------------------------------------------------------------------

% Generate data
[D0, Y] = gen_img_data(m, n, N, type, imdir, images);
Ymask = rand(m, N);
Ymask = (Ymask > missing_data_ratio);
fprintf('Ratio of available data: %.4f\n', sum(sum(Ymask))/m/N)

%%-------------------------------------------------------------------
% AK-SVD on full data

[Df, Xf, errsf] = DL(Y, D0, s, iters, str2func('aksvd'), {}, 'replatom', 'random');

%%-------------------------------------------------------------------
% AK-SVD on incomplete data

params = {'replatoms', 'random'};

%erropts = {'err_missing_data', 'data_mask', Ymask};

%sp_opts = {omp_error_goal, ones(n,1), 0};

[D, X, errmask, errs] = ...
      DL(Y, D0, s, iters, str2func('aksvd_mask'), params, ...
      'spfunc', str2func('omp_mask'), 'data_mask', Ymask, 'outopts', 'last');

% error on full data with dictionary designed for incomplete data
Xx = omp(Y, D, s);
norm(Y-D*Xx,'fro') / sqrt(m*N)
    
%----------------------------------------------------
% Graphs
linewid = 2;
csize = 20;

%---------------------------
figure(1)
plot(errmask, 'b', 'LineWidth', linewid)
hold on
plot(errs, 'r', 'LineWidth', linewid)
plot(errsf, 'k', 'LineWidth', linewid)
hold off
grid
%axis([0 50 3e-4 4e-4])
xlabel('Iterations', 'interpreter', 'latex', 'FontSize', csize );
ylabel('RMSE', 'interpreter', 'latex', 'FontSize', csize );
h=legend('incomplete data only', 'incomplete data extended', 'full data');
set(h, 'interpreter', 'latex', 'FontSize', csize );

%print -depsc fig_aksvd_miss_50.eps
