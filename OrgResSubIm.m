function [Zorg,Zres]=OrgResSubIm(ImName,p,q,KerName,ColorNum)
%
%%%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%%
%
% This function generates p sub-images (one channel among red green blue)
% following different interpolation kernels: linear, cubic convolution
%
%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%
%
% ImName: name of original TIFF image in string format
% p: magnfication rate
% q: minification rate
% KerName: name of used kernel interpolation in string format: 'linear', 'v5cubic'
% ColorNum: red = 1, green = 2, blue = 3
%
%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%
% Zorg: [1,p] cell array for sub-images of original TIFF image
% Zres: [1,p] cell array for sub-images of resampled TIFF image
%
%%%%%%%%%%%%%%% MAIN CODE %%%%%%%%%%%%%%%

% ==== Start the re-sampling process ==== %

% read original TIFF image, take its red channel and get size
Im=imread(ImName); I=double(squeeze(Im(:,:,ColorNum))); [r,c]=size(I);

% magnification with rate p
M=zeros(r,c*p); M(:,1:p:end)=I;

% Interpolation with kernel defined by KerName
U=M;
for j=1:r
    U(j,:)=interp1(1:p:(c*p),I(j,:),1:(c*p),KerName);
end
U=U(:,~all(isnan(U)));

% minification with rate q
D=U(:,1:q:end);  

% ==== End the re-sampling process ==== %

% ==== Start the sub-images extraction ==== %

% pre-allocation
Zorg=cell(1,p); Zres=Zorg;

% define the i-th sub-image
for i=1:p
    Zorg{i}=I(:,i:p:end);
    Zres{i}=D(:,i:p:end);
end

% ==== End the sub-images extraction ==== %

end