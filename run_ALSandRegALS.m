% script    run_AlsandRegALS.m

% clear all
% close all

% Author: Carmeliza Navasca, April 2010 
% This code compares method ALS with PARAFAC 
% Modified by Pentti Paatero for also testing Rot_ALS, May 14, 2010 

%%%%%% choose dimensions %%%%%
I=10;  %I=10;%I=4;
J=9;   %J=10;%J=5;
K=8;   %K=10;%K=10;
R=7;
size_vec=[I J K R];
svdinit=1;   % =1: Initialize using singular vectors
             % =0: Initialize using random vectors

Counts=[];

%load myrandvalues

for taskcount=1:50;

%%%%%% Create data %%%%%%%%
data_type=5;   % non-negative factors
% %data_type=4;
%[A_in,B_in,C_in,X_in]=create_data(data_type,size_vec);

tc=taskcount;
A_in=5*myrand40x10(1:I,        1:R,tc);  % originally  =5*myrand...
B_in=5*myrand40x10(I+1:I+J,    1:R,tc);  % ..
C_in=5*myrand40x10(I+J+1:I+J+K,1:R,tc);  % ..

X_in=tensor_create(A_in,B_in,C_in);


% At this step, you can add noise to X_in 
sigma=10.0;     % Possible values for sigma: e.g. 0.1 to 10.0
% Noise = sigma*randn(I,J,K);

Noise=sigma*myrandn12x11x10(1:I,1:J,1:K,tc);
  
X_in = X_in + Noise;
normN=norm_fro(Noise);

% if you add noise, then use param_vec(3)=1

%%%%% Run ALS %%%%%%%%%%

% Initialize unknown matrices
%A_init= randn(I,R);
%B_init= randn(J,R);
%C_init= randn(K,R);

if (svdinit),   % Initialize unknown matrices by SVD components of X
  % Unfolding the data tensor X into 3 matrices X1 X2 and X3
  X1=tens2mat(X_in,1);
  X2=tens2mat(X_in,2);
  X3=tens2mat(X_in,3);
  [SVDA,S,V]=svd(X1,'econ');
  [SVDB,S,V]=svd(X2,'econ');
  [SVDC,S,V]=svd(X3,'econ');
  A_init= SVDA(:,1:R);
  B_init= SVDB(:,1:R);
  C_init= SVDC(:,1:R);
else;
  A_init=myrand40x10(1:I,        1:R,tc+1);
  B_init=myrand40x10(I+1:I+J,    1:R,tc+1);
  C_init=myrand40x10(I+J+1:I+J+K,1:R,tc+1);
end;

% Choose parameters
% param_vec=[40 1e-5 2];
%param_vec=[10000 1e-4 2];
param_vec=[1000 1e-6 1];
iMax=1000;

tic
[niter1,time1,A_est1,B_est1,C_est1,Fit_vec1,cost1] =  ...
     Method_ALS(X_in,A_init,B_init,C_init,size_vec,param_vec);     % parafac ALS
toc 
tic
[niter,time,A_est,B_est,C_est,Fit_vec,cost]        =  ...
     Method_ALSreg2(X_in,A_init,B_init,C_init,size_vec,param_vec); % parafac ALS + Tikhonov reg
toc


s2=sigma^2;
Counts=  [Counts; tc, niter1, niter, niter2, reject2, round(cost1/s2), round(cost/s2),round(cost2/s2)];

end

disp('sigma, svdinit')
disp([sigma svdinit])

disp(' Itercounts and final Cost for ALS, ALS_Reg, and Rot_ALS')
Counts


%%%%% Exploit Results %%%%%%%%

% Control that the estimated matrices are equal to the input matrices
% up to scaling and permutation

%disp(' Scaling and permutation on A:')
%pinv(A_in)*A_est

%disp(' Scaling and permutation on B:')
%pinv(B_in)*B_est

%disp(' Regul.ALS: Scaling and permutation on C:')
%pinv(C_in)*C_est

%disp(' Scaling and permutation on A:')
%pinv(A_in)*A_est2

%disp(' Scaling and permutation on B:')
%pinv(B_in)*B_est2

%disp(' Rot.ALS: Scaling and permutation on C:')
%pinv(C_in)*C_est2

% disp(['Convergence has been reached in (seconds) '])
% time

% See evolution of the cost function in last fit
figure
semilogy(1:niter,Fit_vec,'*-b');
hold on
semilogy(1:niter1,Fit_vec1,'+-r');

hold off
grid on
legend('Reg ALS', 'ALS')

disp(' Final Cost functions of ALS and ALS_Reg')

[cost1, cost]

