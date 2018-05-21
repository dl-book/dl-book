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

%% Test classification results with discriminative DL and LC-KSVD
clear; close all; fclose all; format compact;
%%-------------------------------------------------------------------------
emitter=[2 4 6 8 10 12 14 16 18 20 22 24 26 28 30 32];     % train overflow
emitter_test=[3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33];% test overflow
nsensors = 2:10;
sensors = {                 % nodes with sensors
[4  12]
[12  15  29]
[9  12  15  30]
[8   9  12  29  31]
[8   9  10  12  24  31]
[7   9  12  24  25  28  31]
[7  12  24  25  26  28  29  31]
[7  12  24  25  26  28  29  30  31]
[14  15  23  24  25  26  28  29  30  31]
};
n = 256;         % dictionary size
s = 4;           % sparsity prior
alpha = 4;       % classification penalty
beta = 16;       % label consistent penalty
init_method = 1; % train small dictionaries for each class at init

% Hanoi partitions built to improve classification. 
parts = {[1 2 3 4 16 17 18], [5 6 7 8 13 14 15], [9 10 11 12],  ...
   [19 20 21 22], [23 24 25 26 27 28 29 30 31]};
%--------------------------------------------------------------------------
load('residues.mat', 'R', 'H_train', 'R_test', 'H_test');
% class_parts, sensors, class_methods
results = zeros(2,length(nsensors),2);

%% Save junction labels
H_train_junc = H_train;
H_test_junc = H_test;
%% Build labels for subgraphs    
samples = length(emitter);
H_train_parts = build_labels(parts, samples);

samples = length(emitter_test);
H_test_parts = build_labels(parts, samples);

for isensors = 1:length(nsensors)
    sensor_nodes = sensors{isensors};
    disp(['Sensor Nodes: ' num2str(sort(sensor_nodes))]);

    for class_parts = [0 1]
        disp(['class_parts=' num2str(class_parts)]);

        %% Build labels for subgraphs    
        if class_parts
            H_train = H_train_parts;
            H_test = H_test_parts;
        else
            H_train = H_train_junc;
            H_test = H_test_junc;
        end

        %% Classification: Learning
        Y_train = double(R(sensor_nodes,:));

        % Compute labels
        c = size(H_train,1);        % number of classes
        nc = floor(n/c);            % evenly divide atoms per classes
        nr = nc*c;                  % total number of atoms (nr <= n)
        Q_train = zeros(nr, size(Y_train,2));
        jj = 0;
        for i = 1 : c                 % allocate atoms for each signal
          jc = find(H_train(i,:)==1); % indices of signals from class i
          Q_train(jj+1:jj+nc,jc) = 1;
          jj = jj + nc;
        end

        % Perform DL
        [W1, D1] = clas_discrim_dl(Y_train, H_train, ...
            nr, s, alpha, init_method);
        [W2, D2, ~] = clas_labelcon_dl(Y_train, H_train, Q_train, ...
            nr, s, alpha, beta, init_method);

        %% Classification: Test data
        Y_test = double(R_test(sensor_nodes,:));

        accuracy1 = classification(Y_test, H_test, D1, W1, s);
        accuracy2 = classification(Y_test, H_test, D2, W2, s);

        results(class_parts + 1, isensors, 1) = accuracy1;
        results(class_parts + 1, isensors, 2) = accuracy2;
        fprintf('discriminative DL %.03f and LC-KSVD %.03f \n', ...
            accuracy1, accuracy2);
    end %% clas_parts
end %% ss

%% Plot
figure(1);
plot(nsensors,results(1, :,1));
hold on;
plot(nsensors,results(1, :,2));
legend({'discriminative DL', 'LC-DL'});
title(['Junctions (c=' num2str(size(H_train_junc,1)) ')']);
hold off;

figure(2);
plot(nsensors,results(2, :,1));
hold on;
plot(nsensors,results(2, :,2));
legend({'discriminative DL', 'LC-DL'});
title(['Partitions (c=' num2str(size(H_train_parts,1)) ')']);
hold off;