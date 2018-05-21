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

%% Basic Algorithms: DL denoising with training on noisy image
%% run_denoising.m needs to be executed before this
%% Assumes BM3D was installed and enabled during run_denoising!
clear; clc; close all; fclose all; format compact;
%%-------------------------------------------------------------------------
p = 8;                  % patch size
s = 6;                  % sparsity
N = 40000;              % total number of patches
n = 256;                % dictionary size
iters = 50;             % DL iterations
replatom = 'worst';     % replace unused atoms
stds = [5 10 20 30 50]; % noise standard deviation
ots = '20171017180427'; % data timestamp, copy from the filename
%%-------------------------------------------------------------------------
updates = {'MOD', 'sgk', 'ksvd', 'aksvd', 'nsgk', 'paksvd', 'pnsgk'};
% Install BM3D and decomment bellow
%dlcompare = [4 8];        % AK-SVD and BM3D
dlcompare = [4 7];        % AK-SVD and PNSGK

% Figures output
figdir = 'fig\';   %racheta
figprefix = 'fig_3_10_denoising';

color = ['y', 'm', 'c', 'r', 'g', 'b', 'k', 'k'];
%%-------------------------------------------------------------------------
datadir = 'data\';
dataprefix = 'denoise';

imdir = 'img\'; %racheta
img_test = {'barbara.png', 'boat.png', 'house.png', 'lena.png', 'peppers.png'};
%%-------------------------------------------------------------------------
ups = length(updates);

f = figure(1);
subplot(2,5,1);
for iimg = 1:length(img_test)
    %% Load
    img = img_test{iimg};
    [~,imgname,~]=fileparts(img_test{iimg});
    allstdpsnr = zeros(length(stds), ups+1);
    allstdssim = zeros(length(stds), ups+1);
    for i = 1:length(stds)
        sigma = stds(i);
        fname = [datadir dataprefix '-' img '-std' num2str(sigma) ...
        '-n' num2str(n) '-' ots '.mat'];
        fprintf('%s sigma=%2d:\n', img, sigma);

        load(fname, 'psnrall', 'ssimall'); % Load DL-based results
        %load(fname, 'psnr_bm3d', 'ssim_bm3d'); % Others
        allstdpsnr(i,1:end-1) = cell2mat(psnrall);
        %allstdpsnr(i,end) = psnr_bm3d;
        allstdssim(i,1:end-1) = cell2mat(ssimall);
        %allstdssim(i,end) = ssim_bm3d;
    end
    
    %% Plot
    %for i = 1:ups+1
    sp=subplot(2,5,iimg);
    hold on;
    for i = dlcompare
        plot(stds(:), allstdpsnr(:,i), color(i), 'Linewidth', 1);
        grid on;
        box on;
        xticks(10:10:50)
        xlim([5,50])
        ylim([25,40.001])
        text(0.5,0.95,imgname,...
                'units','normalized','verticalalignment','top', ...
                'horizontalalignment','center','Interpreter','LaTex')
    end
    hold off;
    
    sp=subplot(2,5,iimg+5);
    hold on;
    for i = dlcompare
        plot(stds(:), allstdssim(:,i), color(i), 'Linewidth', 1);
        grid on;
        box on;
        xticks([5 10:10:50])
        set(sp,'XTickLabels',{'5', '', '20', '', '40', '50'});
        xlim([5,50])
        ylim([0.7,1])
        text(0.5,0.95,imgname,...
                'units','normalized','verticalalignment','top', ...
                'horizontalalignment','center','Interpreter','LaTex')
    end
    hold off;
end

set(gca,'TickLabelInterpreter','latex');
lgd = legend({'AK-SVD', 'BM3D'});