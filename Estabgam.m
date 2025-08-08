function VecPar = Estabgam(SampMu,SampVar,NumPix)
%
%%%%%%%%%%%%% FUNCTION %%%%%%%%%%%%%%
%
% This function estimate the parameters a, b, gamma of the generalized
% signal-dependent noise model form the sample means and sample variances
% of homogeneous segments
%
%%%%%%%%%%%%%%% INPUT %%%%%%%%%%%%%%%
%
% SampMu: colum vector of samples mean of iid pixels in homogeneous segments 
% SampVar: colum vector of samples mean of iid pixels in homogeneous segments 
% NumPix: number of iid pixels in homogeneous segments
%
%%%%%%%%%%%%%%% OUTPUT %%%%%%%%%%%%%%%
%
% VecPar: vector of parameters: 
%         + 1st position of VecPar presents a
%         + 2nd position of VecPar presents b
%         + 3rd position of VecPar presents gam
%
%%%%%%%%%%%%%%% MAIN CODE %%%%%%%%%%%%%%%

% Find initial estimates of a and b by setting gam = 1. In such configuration,
% SampVar = a*SampMu + b, so we can use polynomial curve fitting (least-squares 
% method) to find initial estimates of a and b. 
ab = polyfit(SampMu,SampVar,1);

% initial value of VecPar
VecParIni = [ab,1];

% compute the *minus* log-likelihood function
LLfun = @(p) MinusLogLikeFunc(SampMu,SampVar,p,NumPix); 

% setoptimization options structure
opts = optimset('MaxFunEvals',10000,'MaxIter',10000,'TolX',1e-8,'TolFun',1e-8);

% search estimates of parameters using Nelder-Mead simplex method
VecPar = fminsearch(LLfun,VecParIni,opts);

end

function LL = MinusLogLikeFunc(Mu,Var,Par,Num)
% This function compute the *minus* log-likelihood function, i.e., - Eq. (89) of thesis

% define eps to avoid 0 that make log = -infty
eps = 0.00000000000001;

% define f function describing the relationship Var = f(Mu;a,b,gam)
ffun = @(mu) max(Par(1)./(Par(3)^2).*(mu.^(2-Par(3))) + Par(2)./(Par(3)^2).*(mu.^(2-2*Par(3))),eps ) ;

% define the *minus* log-likelihood function (by removing constant parts wrt Par)
LL = sum( (Num-1).*(log(ffun(Mu)) + Var./ffun(Mu)) );

end





