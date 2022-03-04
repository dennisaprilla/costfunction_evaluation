%======================================================================
% function to compute P, i.e, E-step of CPD registration
%======================================================================
function [ P1, Pt1, Px, negativeLogLikelihood ] = computeEStep( X, Y, sigma2, w )
% The elements m, n of the matrix P are defined as:
%
%                      exp(C1*(Xn-Ym)^2)
% P(m,n) = ----------------------------------------,
%            sum-over-k(exp(C1*(Xn-Yk)^2)) + C2
%
% where C1 = -1/(2*sigma2), M is the number of points in Y, N is the number
% of points in X and D is the dimensionality of points (ex D=3 in case of 3
% dimensional points),
% and   C2 = ( (M*w)/(N*(1-w)) )*((2*pi*sigma2)^(D/2))
%
% The outputs are: P1 = P*1, Pt1 = P'*1, Px = P*X,
%  where 1 is a column vector of all ones.


D = size(X, 2);
M = size(Y, 1);
N = size(X, 1);

% Compute C2 as given in above equation
c   = power(2*pi*sigma2, D/2) *( w /(1-w) *M/N );

% Compute elements of P matrix
pMatrix = zeros(M, N, 'like', X);
for col = 1 : D
    pMatrix = pMatrix + (X(:,col)' - Y(:,col)).^2;
end
pMatrix     = exp(pMatrix*(-1/(2*sigma2)));

% Compute Pt1, P1, Px
pMatrixColSum  = sum(pMatrix);
Pt1         = (pMatrixColSum./(pMatrixColSum +c))';
pMatrix     = pMatrix./(pMatrixColSum +c);
P1          = sum(pMatrix, 2);

Px  = zeros(M, D, 'like', X);
for col = 1 : D
    Px(:, col)  = pMatrix*X(:, col);
end
negativeLogLikelihood  = -sum(log(pMatrixColSum +c)) + D*N*log(sigma2)/2;

end