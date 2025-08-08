function sig=StdEst2D(Z)
%
%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%
%
% This function estimates the noise standard deviation (AWGN model) using  
% Immerkaer's method. (see Alessandro Foi website)
%
%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%
%
% Z: matrix represeting a color channel of a noisy image (in double format)
%
%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%
%
% sig: estimated standard deviation of noise
%
%%%%%%%%%%%%%%% MAIN CODE %%%%%%%%%%%%%%%

LAPL = [1 -2 1;-2 4 -2;1 -2 1];
LAPL = LAPL*sqrt(pi/2/sum(LAPL(:).^2));
YY = conv2(Z,LAPL,'valid');
sig = mean(abs(YY(:)));

end
    
