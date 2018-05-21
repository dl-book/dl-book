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

% Example 6.1/ Error reached for several dictionary sizes
% Generates Figure 6.1 (subject to different random noise)

%-------------------------------------------------------------------------
m = 64;         % problem dimension
nv = [128 192 256 320 384 448 512 576 640 704 768 832 896 960 1024 1152 1280 1408 1536 1664 1792 1920];        % number of atoms in the dictionary
N = 8192;       % number of training signals
sv = [4 6 8 10 12];   % sparsity constraint
iters = 50;     % DL iterations
snr = 20;       % signal-to-noise ratio in dB

update = 'aksvd';      % dictionary update method
replatom = 'random';   % atom replacement strategy

type = 'sliding';      % allow overlapping patches
imdir = [pwd '\img\']; % image directory
images = {'lena.png', 'barbara.png', 'peppers.png'};   % images
%-------------------------------------------------------------------------

% Generate data
[D0, Y] = gen_img_data(m, max(nv), N, type, imdir, images);

% Dictionary Learning
params = {};
errmin = zeros(length(nv), length(sv));
for i = 1 : length(nv)
  for j = 1 : length(sv)
    n = nv(i)
    s = sv(j)
    [D, X, errs] = DL(Y, D0(:,1:n), s, iters, str2func(update), params, 'replatom', replatom);
    errmin(i,j) = min(errs);
  end
end

% Plot error evolution
linewid = 2;
csize = 20;
figure(1)
plot(nv, errmin)
grid
xlabel('$n$', 'interpreter', 'latex', 'FontSize', csize );
ylabel('RMSE', 'interpreter', 'latex', 'FontSize', csize );
h = legend('$s=4$', '$s=6$', '$s=8$', '$s=10$', '$s=12$');
set(h, 'interpreter', 'latex', 'FontSize', csize );

%print -depsc fig_err_size.eps
