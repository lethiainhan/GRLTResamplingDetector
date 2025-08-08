function [Set,SampMu,SampVar,NumPix] = SetExpVarNumSeg(Zc,Zappc,Znoisec,N,NumMin,stp)
%
%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%
%
% This function provides 
% + the set of indices, denoised values, noise values, noisy value of pixels 
%   in homogeneous segments, 
% + the estimated mean and variance of pixels in homogeneous segments
% + the mean of noisy value of pixels in homogeneous segments
% + the number of idd pixels in homogeneous segments 
%
%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%
%
% Zc: matrix of center values represeting noisy color channel of a noisy
%     image (in double format) after 
% Zappc: matrix of center values of approximate image structure after clipping
%        and edge removing
% Znoisec: matrix of center values of extracted noise after clipping and edge 
%          removing
% N: quantization levels (N = 2^8 = 256)
% NumMin: minimum pixels number in homogeneous segments accepted for mean 
%         and variance computation 
% stp: discretization step: stp = 1, 0.5, 0.25, 0.1 ...
%
%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%
%
% Set: set of indices, denoised values, noise values, noisy values of pixels 
%      in homogeneous segements:
%      + 1st column of Set presents segment index
%      + 2nd column of Set presents denoised values
%      + 3rd column of Set presents noise values
%      + 4rd column of Set presents noisy value
% SampMu: colum vector of samples mean of iid pixels in homogeneous segments 
% SampVar: colum vector of samples mean of iid pixels in homogeneous segments 
% NumPix: number of iid pixels in homogeneous segments
%
%%%%%%%%%%%%%%% MAIN CODE %%%%%%%%%%%%%%%

%%%% Transform center values data into column
ColZc = Zc(:); ColZappc = Zappc(:); ColZnoisec = Znoisec(:);

%%%% Set length of homogeneous segments to Del = 1
u = linspace(0,N-1,fix(N/stp)); Del = u(2) - u(1);

%%%% Compute the of indices, denoised values, noise values, noisy values of pixels

% preallocation for Set and preset
Set = zeros(numel(ColZc),4); idx = 1;

% compute elements of Set when the number of pixels in segment > NumMin 
for i = 1:length(u)
    % find indices of pixels in homogeneous segment i
    tmp = find((ColZappc>=u(i)-Del/2) & (ColZappc<u(i)+Del/2));
    % just consider the segment whose number of pixels is greater than NumMin
    if size(tmp,1) >= NumMin 
        Set(idx:idx+size(tmp,1)-1,:) = [ i*ones(size(tmp)) ColZappc(tmp) ColZnoisec(tmp)  ColZc(tmp) ];
        idx = idx + size(tmp,1);
    end
end

% remove untraiteted homogeneous segments (whose the number of pixels < NumMin)
Set = Set(Set(:,1) > 0,:); clear tmp;

%%%% Compute the samples mean and samples variance in each homogeneous segment

% preallocation and preset
SampMu = zeros(length(u),1);
SampVar = zeros(length(u),1);
NumPix = zeros(length(u),1);

% coefficient for "3 sigma" phenomenon of normal distribution
ThreeCoef = 3;

for i = 1:length(u)
    % find index of homogeneous segment
    idx = find(Set(:,1) == i); 
    
    % take denoised, noise, noisy value of pixels in homogeneous segements
    tmp = Set(idx,2:4);
    
    % aberrant values treatment for homogeneous segements 
    if size(tmp,1) >= 2 
        
        % mean of denoised values in a segment (for Gaussian dist, median = mean)
        m1 = median(tmp(:,1));
        
         % standard deviation of denoised values in a segment
        s1 = 1.4826*mad(abs(tmp(:,1)));
        
        % standard deviation of noise values in a segment
        s2 = 1.4826*mad(tmp(:,2),1);
        
        % find index of pixels that do not verify "3 sigma" condition and set to 0
        index1 = find( abs(tmp(:,1) - m1) > ThreeCoef*s1 | abs(tmp(:,2)) > ThreeCoef*s2 );
        Set(idx(index1),:) = 0;
        
        % some specific treatments
        while numel(index1) > 0 && s2 ~= 0 
            index2 = find(abs(tmp(:,1) - m1) <= ThreeCoef*s1 & abs(tmp(:,2)) <= ThreeCoef*s2 );
            tmp = tmp(index2,:);
            idx = idx(index2);     
            m1 = median(tmp(:,1));
            s1 = 1.4826*mad(abs(tmp(:,1)));
            s2 = 1.4826*mad(tmp(:,2),1);       
            index1 = find( abs(tmp(:,1) - m1) > ThreeCoef*s1 |  abs(tmp(:,2)) >= ThreeCoef*s2 );
            Set(idx(index1),:) = 0;
        end
        
        % just consider the segment whose number of pixels is greater than NumMin
        if size(tmp,1) >= NumMin
            NumPix(i) = size(tmp,1);           
            SampMu(i) = mean(tmp(:,1));
            SampVar(i) = var(tmp(:,2));           
        else
            Set(Set(:,1) == i,:) = 0;
        end
        
    end
   
end

% Final results
Set = Set(Set(:,1) > 0,:);
index = find(SampMu > 0 );
SampMu = SampMu(index);
SampVar = SampVar(index);
NumPix = NumPix(index);

end