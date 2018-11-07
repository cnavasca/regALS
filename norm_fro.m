function [frob_norm]=norm_fro(X)
% Calculate the frobenius norm of a tensor or a matrix 
frob_norm=sqrt(sum(sum(sum(abs(X).^2))));
