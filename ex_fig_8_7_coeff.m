% Copyright (c) 2018 Paul Irofti <paul@irofti.net>
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

%% Test coeff spread results with discriminative DL and LC-KSVD
clear; close all; fclose all; format compact;
%%-------------------------------------------------------------------------
emitter=[2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32];     % train overflow
sensor_nodes = [8 9 10 12 24 31]; % sensors placed in the water network

n = 256;         % dictionary size
s = 4;           % sparsity prior
alpha = 4;       % classification penalty
beta = 16;       % label consistent penalty
init_method = 1; % train small dictionaries for each class at init

% Hanoi partitions built to improve classification. 
parts = {[1 2 3 4 16 17 18], [5 6 7 8 13 14 15], [9 10 11 12],  ...
   [19 20 21 22], [23 24 25 26 27 28 29 30 31]};
%--------------------------------------------------------------------------

%% Get residues generated via EPANET emulation
load('residues.mat', 'R', 'H_train', 'R_test');

%% Save junction labels
H_train_junc = H_train;
%% Build labels for subgraphs    
samples = length(emitter);
H_train_parts = build_labels(parts, samples);

p = 1;
subplot(2,2,1);
f = figure(1);

for class_parts = [0 1] % Group a few nodes in one class to reduce errors
    disp(['class_parts=' num2str(class_parts)]);

    %% Build labels for subgraphs    
    if class_parts
        H_train = H_train_parts;
    else
        H_train = H_train_junc;
    end

    %% Classification: Learning
    Y_train = double(R(sensor_nodes,:));
    
    % Compute labels
    c = size(H_train,1);        % number of classes
    nc = floor(n/c);            % evenly divide atoms per number of classes
    nr = nc*c;                  % total number of atoms (nr <= n)
    Q_train = zeros(nr, size(Y_train,2));
    jj = 0;
    for i = 1 : c                 % allocate atoms/labels for each signal
      jc = find(H_train(i,:)==1); % indices of signals from class i
      Q_train(jj+1:jj+nc,jc) = 1;
      jj = jj + nc;
    end
   
    % Perform DL
    [W1, D1] = clas_discrim_dl(Y_train, H_train, ...
        nr, s, alpha, init_method);
    [W2, D2, ~] = clas_labelcon_dl(Y_train, H_train, Q_train, ...
        nr, s, alpha, beta, init_method);

    if class_parts == 1
        c = length(parts);
        title_txt = ['Partitions (c=' num2str(c) ')'];
    else
        c = size(Y_train,2)/length(emitter);
        title_txt = ['Junctions (c=' num2str(c) ')'];
    end

    %% Classification: Test data
    Y_test = double(R_test(sensor_nodes,:));
    X1 = omp(Y_test, D1, s);
    X2 = omp(Y_test, D2, s);

    % Contiguously order signals from the same class
    X1ord = []; X2ord = [];
    for i = 1 : c
      jc = find(H_train(i,:)==1); % indices of signals from class i
      X1c = X1(:, jc);            % those signals
      X2c = X2(:, jc);            % those signals
      X1ord = [X1ord X1c];
      X2ord = [X2ord X2c];
    end
    X1 = X1ord;
    X2 = X2ord;

    %% Plot data
    sp=subplot(2,2,p);
    hold on;
    [ii,jj]=find(X1);
    scatter(jj,ii,1,'b');
    scatter([],[],[],'r');    %dummy for legend
    step = size(X1,2)/c;
    xlim([1,500])
    ylim([1,256])
    set(sp,'XTick',[1 100:100:400 496]);
    set(sp,'YTick',[1 64:64:256]);
    box on;
    title(title_txt);
    hold off;
    p=p+1;

    sp=subplot(2,2,p);
    hold on;
    scatter([],[],[],'b');    %dummy for legend
    [ii,jj]=find(X2);
    scatter(jj,ii,1,'r');
    xlim([1,496])
    ylim([1,256])
    set(sp,'XTick',[1 100:100:400 496]);
    set(sp,'YTick',[1 64:64:256]);
    box on;
    title(title_txt);
    hold off;
    p=p+1;

end %% clas_parts

set(gca, 'TickLabelInterpreter', 'latex');
lgd = legend({'discriminative DL', 'LC-DL'});
set(lgd,'Orientation','horizontal');
set(lgd,'Position',[0.50 0.03 0 0]);