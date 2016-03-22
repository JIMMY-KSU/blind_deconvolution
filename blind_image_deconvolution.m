%% Blind deconvolution
% Notation is based on the paper "Blind Deconvolution using Convex 
% Programming", by Ali Ahmed, Benjamin Recht, and Justin Romberg, from 
% July 22, 2013

%%
close all
clear

%% Generate image
% How does a usefull C_full look like???

% creat_C_full; %creates a matrix C_full, also returns size of the image (s1,s2).
% m_full = zeros(size(C_full,2),1);
% m_full([1 4 5]) = [1 1 1]; % This may be any reasonable vector.
% [ x, C ] = generate_image( C_full, m_full );
% N = size(C,2);

%x = rgb2gray(im2double(imread('images/new/shapes.png')));
x = zeros(20,20);
x(8:12, 8:12) = ones(5,5);
[s1, s2] = size(x);
x = x(:);


L = length(x);
mat = @(x)reshape(x,s1,s2); %s1,s2 are the dimensions of the image
%% Blurr the image.
[ y, B, w_gt, h_gt ] = blurr_image( x, mat );
K = size(B,2);

%% Get subspace information of blurred image or original image
[ C, S, N, m_gt ] = get_subspace( y, L, mat ); %change y to x for original image

%% Set up matrices in Fourier domain
B_hat = fft(full(B)); %fft doesn't take sparse input
C_hat = fft(full(C));
y_hat = fft(y);

%% Define linear operator A
A = zeros(L,K*N);
for i=1:size(C_hat,2)
    Del = diag(sqrt(L)*C_hat(:,i));
    A(:,(i-1)*K+1:i*K) = Del * B_hat;
end

%% Use Boyds cvx solver.

cvx_begin
    variable X(K,N) 
    minimise( norm_nuc(X) )
    subject to
    A*X(:) == y_hat
cvx_end
