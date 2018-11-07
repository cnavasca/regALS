function X_mat = tens2mat(X,mode)

% This function gives the 3 equivalent matrices of a third order tensor 
% by unfolding operation (3 ways of slicing)
% Warning: The equivalent matrices are horizontal here, i.e. X1 (I*JK) X2(J*KI) X3(K*IJ) 
%   with mode=1, the result is such that K varies faster than J along the 
%   rows of X-mat. This may be right but it is essential to know what you get!
%   with mode=2, I varies faster than K
%   with mode=3, J varies faster than I

[I J K]=size(X); 

if mode==1
%Left-right slicing;
%we build a I*JK matrix by concatenation
X1=[];
for j=1:J
    X1=[X1  reshape(X(:,j,:),I,K)];  
end
X_mat=X1;


elseif mode==2
%Front-back slicing
%we build a J*KI matrix 
X2=[];
for k=1:K
    X2=[X2  X(:,:,k).'];     %No need to use reshape X(:,:,k) is already a matrix
end                        
X_mat=X2;


elseif mode==3
%Top-Bottom slicing
% We build a K*IJ matrix
X3=[];
for i=1:I
    X3=[X3 reshape(X(i,:,:),J,K).'];
end
X_mat=X3;

end
