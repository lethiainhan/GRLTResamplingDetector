function VarSig = EstVarSig(SampMu,SampVar,NumPix,VecParOrg)
%
%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%
%
% This function estimate the sum of squared interpolation coefficients
%
%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%
%
% SampMu: colum vector of samples mean of iid pixels in homogeneous segments 
% SampVar: colum vector of samples mean of iid pixels in homogeneous segments 
% NumPix: number of iid pixels in homogeneous segments
% VecParOrg: vector of parameters estimated from an original sub-image
%
%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%
%
% VarSig: the estimate of sum of squared interpolation coefficients
%
%%%%%%%%%%%%%%% MAIN CODE %%%%%%%%%%%%%%%

% define f function describing the relationship Var = f(Mu;a,b,gam)
ffun = @(mu) VecParOrg(1)/(VecParOrg(3)^2)*(mu.^(2-VecParOrg(3))) ...
           + VecParOrg(2)/(VecParOrg(3)^2).*(mu.^(2-2*VecParOrg(3)));

% compute the estimate of varsig
VarSig = sum( (NumPix-1).*SampVar./ffun(SampMu) )/ sum(NumPix-1);

end




