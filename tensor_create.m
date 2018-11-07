%*************************************************************
% create tensor from factor matrices: A, B, C
%*************************************************************
function X=tensor_create(A,B,C)

% A is I*R
% B is J*R
% C is K*R
% D is 1*R

I=size(A,1);
J=size(B,1);
K=size(C,1);
R=size(A,2);

X=zeros(I,J,K);


%*************************
%We build the tensor X
%*************************

    for k=1:K
        X(:,:,k)=A*diag(C(k,:))*B.';
    end
   

