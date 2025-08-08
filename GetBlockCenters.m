function [Zc,Zappc,Znoisec,Bc,Clc,Xc,BCXc]=GetBlockCenters(Z,Zapp,Znoise,B,Cl,X,BCX,m,n,offsetmn,bsize)
%
%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%
%
% This function get the center values of Z, Zapp, Znoise, B, Cl, X, BCX. 
% These values is defined by the position [m,n] of blocks with size bsize
% and by the offsetmn <= bsize/2.
%
% if m=4, n=4, offsetmn=1 and bsize = 8, then the center values are at the 
% positions [4,4], [4,5], [5,4] and [5,5] of 8x8 blocks
%
% It is recommnded to get more center values than the ones at positions [4,4]
% in order to have more idd pixels data
%
%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%
%
% Z: matrix represeting noisy color channel of a noisy image (in double format)
% Zapp: approximate image structure after denoising by wavelet-based filter 
% Znoise: noise extracted using wavelet-based filter given from NoiseExtractChannel
% B: 8x8 blocks of noisy image channel
% Cl: index matrix for pixels retained after clipping
% X: index matrix for pixels impacted by edge
% BCX: index matrix for 8x8 blocks of noisy image channel after clipping and edge detection
% m,n,offsetmn: parameters defining center positions of blocks with size bsize 
% bsize: size of consider block
%
%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%
%
% Zc: matrix of center values associated with Z
% Zappc: matrix of center values associated with Zappc
% Znoisec: matrix of center values associated with Znoisec
% Bc: matrix of center values associated with B
% Clc: matrix of center values associated with Cl
% Xc: matrix of center values associated with X
% BCXc: matrix of center values associated with BCX
%
%%%%%%%%%%%%%%% MAIN CODE %%%%%%%%%%%%%%%

% the following code is designed for m=4, n=4, offsetmn=1 and bsize = 8.

switch offsetmn
    case 0 
        Zc = Z(m:bsize:end,n:bsize:end); 
        Zappc = Zapp(m:bsize:end,n:bsize:end); 
        Znoisec = Znoise(m:bsize:end,n:bsize:end); 
        Bc = B(m:bsize:end,n:bsize:end); 
        Clc = Cl(m:bsize:end,n:bsize:end); 
        Xc = X(m:bsize:end,n:bsize:end); 
        BCXc = BCX(m:bsize:end,n:bsize:end); 
        
    case 4
        Zc = Z(1:end,1:end);
        Zappc = Zapp(1:end,1:end);
        Znoisec = Znoise(1:end,1:end);   
        Bc = B(1:end,1:end);
        Clc = Cl(1:end,1:end);
        Xc = X(1:end,1:end);
        BCXc = BCX(1:end,1:end);
        
    otherwise % only 1 2 3
        
        % preallocation
        Zc=zeros(size(Z,1)/bsize*(2*offsetmn),size(Z,2)/bsize*(2*offsetmn));
        Zappc=zeros(size(Z,1)/bsize*(2*offsetmn),size(Z,2)/bsize*(2*offsetmn));
        Znoisec=zeros(size(Z,1)/bsize*(2*offsetmn),size(Z,2)/bsize*(2*offsetmn));
        Bc=zeros(size(Z,1)/bsize*(2*offsetmn),size(Z,2)/bsize*(2*offsetmn));
        Clc=zeros(size(Z,1)/bsize*(2*offsetmn),size(Z,2)/bsize*(2*offsetmn));
        Xc=zeros(size(Z,1)/bsize*(2*offsetmn),size(Z,2)/bsize*(2*offsetmn));
        BCXc=zeros(size(Z,1)/bsize*(2*offsetmn),size(Z,2)/bsize*(2*offsetmn));
        
        % get values
        for ci=1:(2*offsetmn)
            for ri=1:(2*offsetmn)
                Zc(ri:(2*offsetmn):end,ci:(2*offsetmn):end)=Z(m-offsetmn+ri:bsize:end,n-offsetmn+ci:bsize:end);
                Zappc(ri:(2*offsetmn):end,ci:(2*offsetmn):end)=Zapp(m-offsetmn+ri:bsize:end,n-offsetmn+ci:bsize:end);
                Znoisec(ri:(2*offsetmn):end,ci:(2*offsetmn):end)=Znoise(m-offsetmn+ri:bsize:end,n-offsetmn+ci:bsize:end);
                Bc(ri:(2*offsetmn):end,ci:(2*offsetmn):end)=B(m-offsetmn+ri:bsize:end,n-offsetmn+ci:bsize:end);
                Clc(ri:(2*offsetmn):end,ci:(2*offsetmn):end)=Cl(m-offsetmn+ri:bsize:end,n-offsetmn+ci:bsize:end);
                Xc(ri:(2*offsetmn):end,ci:(2*offsetmn):end)=X(m-offsetmn+ri:bsize:end,n-offsetmn+ci:bsize:end);
                BCXc(ri:(2*offsetmn):end,ci:(2*offsetmn):end)=BCX(m-offsetmn+ri:bsize:end,n-offsetmn+ci:bsize:end);
            end
        end
        
end

end