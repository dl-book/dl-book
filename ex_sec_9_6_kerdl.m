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

% Examples described in Section 9.6: classification with kernel
% methods for a simple two-class problem.
%
% DL algorithms: aksvd_ker - kernel AK-SVD (Algorithm 9.2)

% BD 25.12.2017

type_clas = 1;  % classifier used for the actual classification
                % 0 - discriminative DL, like in Section 9.6.3
                % 1 - SRC with kernel DL, like in Section 9.6.1-2 (used for figure 9.2 in the book!)

% general parameters
m = 2;
c = 2;        % number of classes
Nc = 100;     % number of training signals per class
N = c*Nc;
n = 8;        % number of atoms
nv = [8; 8];  % number of atoms for SRC
s = 1;        % sparsity level
nstd = 0.1;   % noise power (for training signals)
nstdt = 0.25; % noise power for test signals
gen_new_signals = 1; % set to 0 to run again on same data

% kernel parameters
ker_code = 0;
switch ker_code
  case 0,     % RBF
    sigma = 1;
    mykernel = @(x,y) exp( - (x-y)'*(x-y)/2/sigma/sigma);
  case 1,     % polynomial
    a = 1;
    d = 2;
    mykernel = @(x,y) (x'*y + a)^d;
end

% create classes
if gen_new_signals
  Y = zeros(m,N);
  H = zeros(c,N);
  for i = 1:Nc
    Y(:,i) = nstd*randn(m,1);
    v = randn(m,1);
    v = v/norm(v);  % norm 1, random direction
    Y(:,Nc+i) = v*(1 + nstd*randn);
    H(1,i) = 1;
    H(2,Nc+i) = 1;
  end
end

% train classifier
switch type_clas
  case 0,       % discriminative DL
    alpha = 0.1;
    iternum = 50;
    A = randn(N,n);           % initial dictionary
    A = normc(A);             % incorrect normalization !!!
    W = zeros(2,n);
    [A, W, X] = aksvd_ker_discr(Y, A, H, W, alpha, s, mykernel, iternum);
  case 1,       % SRC with kernel DL
    A = clas_src_dl_ker(Y, H, nv, s, mykernel);
    n = sum(nv);
    A1 = A(:, 1:nv(1));
    A2 = A(:, nv(1)+1:end);
end

% test result on training signals
K = zeros(N,N);
for i = 1:N         % compute kernel matrix from training signals
  for j = 1:i
    K(i,j) = mykernel(Y(:,i), Y(:,j));
    K(j,i) = K(i,j);
  end
end
K1 = K(1:Nc,1:Nc);
K2 = K(Nc+1:N,Nc+1:N);

% test classifier
Ntc = 2000;       % number of test signals for each class
Nt = c*Ntc;
Yt = zeros(m,Nt);
for i = 1:Ntc  % generate test signals
  Yt(:,i) = nstdt*randn(m,1);
  v = randn(m,1);
  v = v/norm(v);  % norm 1, random direction
  Yt(:,Ntc+i) = v* (1+nstdt*randn);
end

c1 = [];    % vectors with index signals for each class
c2 = [];
switch type_clas
  case 0,       % kernel discriminative DL
    for i = 1 : Nt  % test signals one by one, for memory reasons
      kz = zeros(N,1);
      for j = 1 : N  % kernelized scalar products with training signals
        kz(j) = mykernel(Yt(:,i), Y(:,j));
      end
      x = omp_ker(K, kz, A, 1);     % representation of test signal
      h = W*x;                      % label vector
      if h(1) > h(2)
        c1 = [c1 i];
      else
        c2 = [c2 i];
      end
    end
  case 1,       % SRC with kernel DL
    for i = 1 : Nt
      z = Yt(:,i);
      k1 = zeros(1,Nc);
      k2 = zeros(1,Nc);
      for j = 1 : Nc  % since s=1, we only need to find which atom is closest from the signal
        k1(j) = mykernel(z, Y(:,j));
        k2(j) = mykernel(z, Y(:,Nc+j));
      end
      if max(abs(k1*A1)) > max(abs(k2*A2))
        c1 = [c1 i];
      else
        c2 = [c2 i];
      end
    end
end

% plot training points
figure(1)
plot(Y(1,1:Nc), Y(2,1:Nc), 'or')
hold on
plot(Y(1,Nc+1:N), Y(2,Nc+1:N), 'og')
hold off
grid
%print -depsc fig_train_sig.eps

% plot test points
figure(2)
plot(Yt(1,c1), Yt(2,c1), '.r')
hold on
plot(Yt(1,c2), Yt(2,c2), '.g')
hold off
grid
axis([-1.5 1.5 -1.5 1.5])
%print -depsc fig_kerdl_test_8.eps
