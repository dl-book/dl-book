% Copyright (c) 2017 Bogdan Dumitrescu <bogdan.dumitrescu@acse.pub.ro>
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

%% Error reached for several initializations
%clear; clc; close all; fclose all; format compact;
%%-------------------------------------------------------------------------
m = 64;         % problem dimension
n = 128;        % number of atoms in the dictionary
N = 4000;       % number of training signals
s = 6;          % sparsity constraint
iters = 50;     % DL iterations
Ninit = 5;      % number of different initializations

update = 'aksvd';    % dictionary update method
replatom = 'random'; % atom replacement strategy

type = 'distinct';   % 'distinct' or 'sliding' patches
imdir = 'img\';      % image directory
images = {'barbara.png', 'boat.png', 'house.png', 'lena.png', 'peppers.png'};
%%-------------------------------------------------------------------------

% Generate data
[D0, Y] = gen_img_data(m, n, N, type, imdir, images);

% Dictionary Learning
params = {};
for i = 1 : Ninit
  D0 = normc(randn(m,n));   % random vectors
  [D, X, errsr{i}] = DL(Y, D0, s, iters, str2func(update), params, 'replatom', replatom);
end

for i = 1 : Ninit
  D0 = normc(Y(:, randperm(N,n)));  % random signals
  [D, X, errss{i}] = DL(Y, D0, s, iters, str2func(update), params, 'replatom', replatom);
end

% Plot error evolution
linewid = 2;
% colors = {'b','y'};
csize = 20;

%---------------------------
figure(1)
for i = 1 : Ninit
  plot(errsr{i})
  hold on
end
hold off
grid
axis([0 50 4e-4 5e-4])
%axis([0 50 3e-4 4e-4])
xlabel('Iterations', 'interpreter', 'latex', 'FontSize', csize );
ylabel('RMSE', 'interpreter', 'latex', 'FontSize', csize );

%---------------------------
figure(2)
for i = 1 : Ninit
  plot(errss{i})
  hold on
end
hold off
grid
axis([0 50 4e-4 5e-4])
%axis([0 50 3e-4 4e-4])
xlabel('Iterations', 'interpreter', 'latex', 'FontSize', csize );
ylabel('RMSE', 'interpreter', 'latex', 'FontSize', csize );

