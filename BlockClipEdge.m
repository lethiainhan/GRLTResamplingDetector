function [B,Cl,X,BCX]=BlockClipEdge(Z,Zapp,Znoise,N)
%
%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%
%
% This function converts an noisy image channel to blocks, then applies the 
% clipping and edge removing to obtain index of homogeneous blocks
%
%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%
%
% Z: matrix represeting noisy color channel of a noisy image (in double format)
% Zapp: approximate image structure after denoising by wavelet-based filter 
% Znoise: noise extracted using wavelet-based filter given from NoiseExtractChannel
% N: quantization levels (N = 2^8 = 256)
%
%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%
%
% B: 8x8 blocks of noisy image channel
% Cl: index matrix for pixels retained after clipping
% X:  index matrix for pixels impacted by edge
% BCX: index matrix for 8x8 blocks of noisy image channel after clipping and edge detection
%
%%%%%%%%%%%%%%% MAIN CODE %%%%%%%%%%%%%%%

% Convert noisy image channel to blocks (homogeneous blocks)
B = im2vec(Z,8);
B = mean(B); 
B = repmat(B,64,1);
B = im2vec(Z,8) - B;
B = B ~= 0; 
B = sum(B); B = B';
B = reshape(B,size(Z,1)/8,size(Z,2)/8);
B = kron(B,ones(8)); 

% find index matrix for pixels retained after clipping
Cl = Z <= 0.1*N | Z >= 0.9*N;
Cl = im2vec(Cl,8);
Cl = sum(Cl); 
Cl = Cl';
Cl = reshape(Cl,size(Z,1)/8,size(Z,2)/8);
Cl = kron(Cl,ones(8)); 

% index matrix for pixels impacted by edge
signoise = 1.4826*mad(Znoise(:),1);
x = bdct(Zapp);
x = im2vec(x,8);
x = x(2:64,:);
x = 1.4826*mad(x,1);
X = x > signoise; 
y = x(X==0);
if numel(y) > 1
    Y = sort(y); 
    seuil = Y(floor(0.9*numel(y)));
    X = x > seuil;
end
X = X';
X = reshape(X,size(Zapp,1)/8,size(Zapp,2)/8);
X = kron(X,ones(8));

% define index matrix for 8x8 blocks after clipping and edge detection
BCX = (X == 0 & Cl == 0 & B~=0);

end