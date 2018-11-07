function C = kat_rao(A,B)
% Calculate the Kathri Rao product of 2 matrices 
[A1,A2]=size(A);
[B1,B2]=size(B);
if A2~=B2,
    error('Kathri Rao product: # of columns in A and B are different')
else
    C=zeros(A1*B1,A2);
    for j=1:A2,
        C(:,j)=kron(A(:,j),B(:,j));     % column wise kronecker product
    end
end
