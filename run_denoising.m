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
function run_denoising()
clear; clc; close all; fclose all; format compact;
%%-------------------------------------------------------------------------
p = 8;                  % patch size
s = 6;                  % sparsity
N = 40000;              % total number of patches
n = 256;                % dictionary size
iters = 50;             % DL iterations
replatom = 'worst';    % replace unused atoms
stds = [5 10 20 30 50]; % noise standard deviation
%%-------------------------------------------------------------------------
updates = {'MOD', 'sgk', 'ksvd', 'aksvd', 'nsgk', 'paksvd', 'pnsgk'};
methods = [
  % Name	Function        Dictionary index
  {'modomp ', @denoise_omp,	1};
  {'sgkomp ', @denoise_omp,	2};
  {'ksvdomp', @denoise_omp,	3};
  {'akomp', @denoise_omp,	4};
  {'nsgkomp', @denoise_omp,	5};
  {'pakomp', @denoise_omp,	6};
  {'pnsomp', @denoise_omp,	7};
];
%%-------------------------------------------------------------------------
datadir = 'data\';   %racheta
dataprefix = 'denoise';

imdir = 'img\'; %racheta
img_test = {'barbara.png', 'boat.png', 'house.png', 'lena.png', 'peppers.png'};
%%-------------------------------------------------------------------------
addpath(genpath('DL'));
%addpath('BM3D')
ts = datestr(now, 'yyyymmddHHMMss');
%%-------------------------------------------------------------------------
% INITIALIZATION
%%-------------------------------------------------------------------------
ups = length(updates);
Dall = cell(ups,1);
Xall = cell(ups,1);
D0 = odctdict(p^2,n);
%%-------------------------------------------------------------------------
% DENOISING
%%-------------------------------------------------------------------------
f = {'name', 'func', 'dict'};
m = cell2struct(methods, f, 2);
%-------------------------------------------------------------------------
results = sprintf('%s\n', dataprefix);
Yall = {};
Xall = {};
Iall = {};
psnrall = {};
ssimall = {};
%-------------------------------------------------------------------------
for iimg = 1:length(img_test)
    img = img_test{iimg};
    for sigma = stds
        do_img(img, sigma);        
        clc; disp(results);
    end
end
%-------------------------------------------------------------------------
function do_img(img, sigma)
    fname = [datadir dataprefix '-' img '-std' num2str(sigma) ...
        '-n' num2str(n) '-' ts '.mat'];
    results = sprintf('%s%s sigma=%2d:\n', results, img, sigma);
    clc;disp(results);

    %% Initial Data
    [I, Inoisy, Ynoisy, Ynmean] = ...
        denoise_init_data([imdir,char(img)], sigma, p, p);    
    save(fname, 'Inoisy', 'Ynoisy', 'Ynmean');
    
    function do_denoise(name, dfunc, D, i)
        [Iall{i}, Yall{i}, Xall{i}] = ...
            denoise(dfunc, Inoisy, D, {'sigma', sigma, 's', s}, p, p);
        psnrall{i} = psnr(Iall{i},I,255);
        ssimall{i} = ssim(Iall{i}, I, 'DynamicRange', 255);
        results = sprintf('%s %s psnr=%f ssim=%f\n', ...
            results, name, psnrall{i}, ssimall{i});
        clc; disp(results);
    end

    Dall = cell(ups,1);
    Xtrainall = cell(ups,1);
    errsall = zeros(ups,iters);
    for i = 1:ups
        fprintf('%s', updates{i}(1));
        [Dall{i}, Xtrainall{i}, errsall(i,:)] = ...
            DL(Ynoisy(:,randperm(size(Ynoisy,2), N)), D0, s, iters, ...
            str2func(updates{i}), 'replatom', replatom);
    end
    for i = 1:size(methods,1)
        do_denoise(m(i).name, m(i).func, Dall{m(i).dict}, i);
    end
    
    %% Save DL-based results
    save(fname, 'Yall', 'Xall', 'Iall', ...
        'Dall', 'Xtrainall', 'errsall', ...
        'psnrall', 'ssimall', '-append');

    %% BM3D: enable after installing
    %{
    [~, Ibm3d] = BM3D(1, Inoisy, sigma);
    Ibm3d = Ibm3d*255; % undo BM3D scaling
    psnr_bm3d = psnr(Ibm3d,I,255);
    ssim_bm3d = ssim(Ibm3d, I, 'DynamicRange', 255);
    results = sprintf('%s bm3d    psnr=%f ssim=%f\n', ...
        results, psnr_bm3d, ssim_bm3d);
    save(fname, 'Ibm3d', 'psnr_bm3d', 'ssim_bm3d', '-append');
    clear Ibm3d pnsr_bm3d ssim_bm3d
    clc; disp(results);
    %}
end
end