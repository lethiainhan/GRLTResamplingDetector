function Znoise = NoiseExtractChannel(Z,sig)
%
%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%
%
% This function extracts the residue (noise) from an image channel using 
% wavelet-based filter
% Mihcak, M.K., Kozintsev, I., Ramchandran, K., 1999. Spatially adaptive
% statistical modeling of wavelet image coefficients and its application
% to denoising, in: 1999 IEEE International Conference on Acoustics,
% Speech and Signal Processing. (ICASSP), IEEE. pp. 3253–3256.
%
%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%
%
% Z: matrix represeting noisy color channel of a noisy image (in double format)
% sig: estimated standard deviation of noise given by 4*StdEst2D(Z) 
%
%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%
%
% Znoise: residue (noise) extracted using wavelet-based filter
%
%%%%%%%%%%%%%%% MAIN CODE %%%%%%%%%%%%%%%

% define the number of wavelet decomposition levels (L = 3 or 4)
L = 4;

% define the scaling function of an orthogonal wavelet filter
qmf = MakeONFilter('Daubechies',8);

% get size of image channel
[M,N] = size(Z); m = 2^L;

% use padding with mirrored image content 
minpad=2;    % minimum number of padded rows and columns as well
nr = ceil((M+minpad)/m)*m;  nc = ceil((N+minpad)/m)*m;  % dimensions of the padded image (always pad 8 pixels or more)
pr = ceil((nr-M)/2);      % number of padded rows on the top
prd= floor((nr-M)/2);     % number of padded rows at the bottom
pc = ceil((nc-N)/2);      % number of padded columns on the left
pcr= floor((nc-N)/2);     % number of padded columns on the right
Z = [Z(pr:-1:1,pc:-1:1),     Z(pr:-1:1,:),     Z(pr:-1:1,N:-1:N-pcr+1);
    Z(:,pc:-1:1),           Z,                Z(:,N:-1:N-pcr+1);
    Z(M:-1:M-prd+1,pc:-1:1),Z(M:-1:M-prd+1,:),Z(M:-1:M-prd+1,N:-1:N-pcr+1)];
    
% pre-compute noise variance
NoiseVar = sig^2;   

% wavelet transform
[wave_trans,L] = mdwt(Z,qmf,L);

% extract the noise from the wavelet coefficients 

for i=1:L
    % indicies of the block of coefficients
    Hhigh = (nc/2+1):nc; Hlow = 1:(nc/2);
    Vhigh = (nr/2+1):nr; Vlow = 1:(nr/2);
   
    % Horizontal noise extraction
    wave_trans(Vlow,Hhigh) = WaveNoise(wave_trans(Vlow,Hhigh),NoiseVar);
         
    % Vertical noise extraction
    wave_trans(Vhigh,Hlow) =  WaveNoise(wave_trans(Vhigh,Hlow),NoiseVar);
      
    % Diagonal noise extraction
    wave_trans(Vhigh,Hhigh) = WaveNoise(wave_trans(Vhigh,Hhigh),NoiseVar);
	
    % re-define the dimension
    nc = nc/2; nr = nr/2;
end

% Last, coarest level noise extraction
wave_trans(1:nr,1:nc) = 0;

% Inverse wavelet transform
[image_noise,~] = midwt(wave_trans,qmf,L);

% Crop to the original size
Znoise = image_noise(pr+1:pr+M,pc+1:pc+N);

end