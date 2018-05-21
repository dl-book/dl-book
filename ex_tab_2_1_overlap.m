% Copyright (c) 2017, 2018 Paul Irofti <paul@irofti.net>
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

%% Basic Algorithms: DL denoising overlapping patches example
%% run_denoising.m needs to be executed before this
clc; close all; fclose all; format compact;
%%-------------------------------------------------------------------------
p = 8;                  % patch size
s = 6;                  % sparsity
N = 40000;              % total number of patches
n = 256;                % dictionary size
iters = 50;             % DL iterations
replatom = 'worst';    % replace unused atoms
stds = [5 10 20 30 50]; % noise standard deviation
ots = '20171017180427'; % data timestamp, copy from the filename
%%-------------------------------------------------------------------------
iimg = 4;               % 'lena'
istd = 3;               % '20dB'
imethod = 4;            % 'aksvd'
%%-------------------------------------------------------------------------
updates = {'MOD', 'sgk', 'ksvd', 'aksvd', 'nsgk', 'paksvd', 'pnsgk'};

% Figures output
figdir = 'fig\';
figprefix = 'tab_2_1_overlap';
%%-------------------------------------------------------------------------
datadir = 'data\';
dataprefix = 'denoise';

imdir = 'img\';
img_test = {'barbara.png', 'boat.png', 'house.png', 'lena.png', 'peppers.png'};
%%-------------------------------------------------------------------------
img = img_test{iimg};
sigma = stds(istd);             % noise standard deviation

%% Load
fname = [datadir dataprefix '-' img '-std' num2str(sigma) ...
        '-n' num2str(n) '-' ots '.mat'];
load(fname, 'Inoisy', 'Iall', 'Yall'); % Load DL-based results
I = double(imread([imdir img]));
I = I(:, :, 1);

%% Build patches
Y = Yall{imethod};
Iclean = Iall{imethod};
nblocks = 9;
npatches = 8;
Idim = 512;
pos = round(size(Y,2)/2) - 25*(Idim - (p-1)) - 30;

ipos = 265-30;
r = ipos:ipos+(npatches-1)*p;
I3 = {I(r,r), Inoisy(r,r), Iclean(r,r)};
I3title = {'Original', 'Noisy', 'Cleaned'};

for i = 1:length(I3)
    imshow(I3{i},[0 255]);
    mkpdf(figdir,[figprefix '_' I3title{i}]);
end

set(0,'defaultTextInterpreter','latex');
%subplot(sqrt(nblocks),sqrt(nblocks)+1,1)
for i = 0:sqrt(nblocks)-1
    %subplot(sqrt(nblocks),sqrt(nblocks)+1,i*(sqrt(nblocks)+1) + 1);
    %imshow(I3{i+1},[0 255]);
    %mkpdf(figdir,[figprefix '_' num2str(i) num2str(j) '.pdf']);
    for j = 0:sqrt(nblocks)-1
        %subplot(sqrt(nblocks),sqrt(nblocks)+1,i*(sqrt(nblocks)+1) + j + 1 + 1)
        pos = pos + j*(Idim - (p - 1)) + i;
        Iij = sliding2im(Y, Idim, pos, npatches, p);
        imshow(Iij,[0 255])
        mkpdf(figdir,[figprefix '_' num2str(i) num2str(j)]);
    end
end