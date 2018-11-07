function [niter,time,A_est,B_est,C_est,Fit_vec,cost]=Method_ALSreg2(X,A_init,B_init,C_init,size_vec,param_vec)

% Author: Carmeliza Navasca  
% This is the ALS with a Tikhonov regularization term, Feb 2008.

%  K    _______             
%      /      /|         a1/             a2/                       aR/
%     /______/ |          /________   +   /________   +   ---   +   /________
%    |  J   |  |          |   b1          |   b2                    |   bR
%    |      | /    -->  c1|             c2|                       cR|
% I  |______|/            |               |                         |
%        X 

% The method used here is an Alternating Least Squares Algorithm (ALS)
% The aim is to calculate the following unknown matrices that have been used to build X,
% by using alternating least squares conditional updates: 
% A=[a1 ... aR] of size (I,R)
% B=[b1 ... bR] of size (J,R)
% C=[c1 ... cR] of size (K,R)

%%%%%%%%%% Inputs %%%%%%%%%%%
%  -- X: tensor to be decomposed, of size (I,J,K)
%  -- B_init, C_init: initialisation of unknown matrices 
%  (only B_init and C_init are needed if you start with the update of A) 
%  -- size_vec = [I J K R]
%  -- param_vec: user parameters with the 3 following scalars:
%     param_vec(1): time allowed to reach convergence (e.g 40)
%     param_vec(2): convergence criterion, i.e. tolerance to decide whether
%     convergence has been reached (e.g. 1e-5)
%     param_vec(3): convergence method (1 or 2)
%           --> if 1: we compare X_est to X_est of the previous iteration,
%                   to evaluate the relative amount of correction brought by the
%                   current iteration
%           --> if 2: we compare X_est to the data tensor X itself 
%                   good criterion when an exact solution exists, 
%                   i.e. there is no additive noise so that the model can exactly be fitted

%%%%%%%%%%% Outputs %%%%%%%%%%
% -- niter: number of iterations
% -- time: duration of the calculation
% -- A_est,B_est,C_est: estimation of unknown matrices
% -- Fit_vec: evolution of the error


% Take the dimensions
I=size_vec(1);
J=size_vec(2);
K=size_vec(3);
R=size_vec(4);

% Initialise the estimated matrices and estimated model
B_est=B_init;
C_est=C_init;
A_est=A_init;
X_est=tensor_create(A_est,B_est,C_est);

% Take the parameters
time_max=param_vec(1);      % when running time exceeds time_max, 
                            % we stop the algorithm because we probably reached a local minimu   
conv_crit=param_vec(2);     % convergence criterion for global algorithm
conv_method=param_vec(3);   % flag to choose convergence criterion for global algorithm
                                
                                
% Initialise variables
niter=0;
time=0;
FitNew=norm_fro(X-X_est)^2;

% Unfolding the tensor of data X into 3 matrices X1 X2 and X3
X1=tens2mat(X,1);
X2=tens2mat(X,2);
X3=tens2mat(X,3);


%%%%%%  MAIN  LOOP %%%%%%%%% 
    
    tic
    Fit_vec=[];
    alpha=10^(0);

  while (FitNew>conv_crit) && (time<time_max)
   
   FitNew;   
   niter=niter+1;                       % number of iterations
   X_est_old=X_est;                     % store old estimate of the model
    
   % Alternating updates
       
       %size(X1)
       newX1=[X1 alpha*A_est];
       
       D1=[kat_rao(B_est,C_est)' alpha*eye(R) ];
       A_est=newX1/(D1);   % update A  
       
       newX2=[X2 alpha*B_est];
       D2=[kat_rao(C_est,A_est)' alpha*eye(R)];
       B_est=newX2/(D2); % update B
       
       newX3=[X3 alpha*C_est];
       D3=[kat_rao(A_est,B_est)' alpha*eye(R)];
       C_est=newX3/(D3); % update C
   
   % Normalisation step: normalize in the 3 modes equally
    normA=sqrt(sum(A_est.*(A_est').'));
    normB=sqrt(sum(B_est.*(B_est').'));
    normC=sqrt(sum(C_est.*(C_est').'));
    prod=normA.*normB.*normC;

    A_est=A_est*diag(1./normA)*diag(prod.^(1/3));
    B_est=B_est*diag(1./normB)*diag(prod.^(1/3));
    C_est=C_est*diag(1./normC)*diag(prod.^(1/3));

    % Build estimated model from normalized updated matrices
    X_est = tensor_create(A_est,B_est,C_est);
   
    
    % Define Stop Criterion
       if conv_method==1         % we compare X_est to X_est_old
  	     FitNew=norm_fro(X_est_old-X_est)^2;            
      elseif conv_method==2     % we compare X_est to X
         FitNew=norm_fro(X-X_est)^2;
      end

      Fit_vec=[Fit_vec FitNew];     % To see the evolution of the error
      time=toc;    
      if time>100000, time=0;
      end;
      alpha=.95*alpha; %alpha=0.85*alpha; %good for small scale 
  end
  cost = norm_fro(X-X_est)^2

