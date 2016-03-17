%% Blind deconvolution
% Notation is based on the paper "Blind Deconvolution using Convex 
% Programming", by Ali Ahmed, Benjamin Recht, and Justin Romberg, from 
% July 22, 2013

%%
close all
clear

%% Generate image
% How does a usefull C_full look like???

creat_C_full; %creates a matrix C_full, also returns size of the image (s1,s2).
m_full = zeros(size(C_full,2),1);
m_full([1 4 5]) = [1 1 1]; % This may be any reasonable vector.
[ x, C ] = generate_image( C_full, m_full );
N = size(C,2);

L = length(x);
mat = @(x)reshape(x,s1,s2);
%% Blurr the image.
[ y, B, w, h ] = blurr_image( x, mat );
K = size(B,2);

%% Set up matrices in Fourier domain
B_hat = fft(full(B)); %fft doesn't take sparse input
C_hat = fft(C);
y_hat = fft(y);

%% Define linear operator A
A = zeros(L,K*N);
for i=1:size(C_hat,2)
    Del = diag(sqrt(L)*C_hat(:,i));
    A(:,(i-1)*K+1:i*K) = Del * B_hat;
end

%% Use Boyds cvx solver.

cvx_begin
    variable X(K*N) 
    minimise( norm_nuc(reshape(X,K,N)) )
    subject to
    A*X(:) == y_hat
cvx_end
