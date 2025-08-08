% GLRT based detector

clear

%%%%% define parameters
N = 2^8; m = 4; n = 4; offsetmn = 4; bsize = 8; NumMin = 100; stp = 1;
p = 4; q = 3; KerName = 'v5cubic'; ColorNum = 1; alp0 = 0.5; 

%%%%% define original image and generate resampling image related to a color channel 
ImName = 'NikonD200.tiff'; 
[Zorg,Zres]=OrgResSubIm(ImName,p,q,KerName,ColorNum);

%%%%% traitement for each kind of image
%%% prealocation
LLR = cell(2,p); VarSig = LLR; the = LLR; bet = LLR; deci = LLR;
for i = 1:2
    %%% Get original image
    if i == 1
        Zim = Zorg;
    %%% Get resampled image
    else
        Zim = Zres;
    end
        
    %%%% traitement for each sub-image
    for j = 1:p
        
        %%% get noisy sub-image and resize to be multiple of 8x8 block
        Zsub = Zim{j};
        [r,c] = size(Zsub); 
        Zsub = Zsub(1:8*floor(r/8),1:8*floor(c/8));
        
        %%% clear unnecessary variables
        clear r c;
        
        %%% get noise and denoised sub-image related to the noisy sub-image
        sig=4*StdEst2D(Zsub);
        Zsubnoise = NoiseExtractChannel(Zsub,sig);
        Zsubapp = Zsub - Zsubnoise;
        
        %%% clear unnecessary variables
        clear sig;
        
        %%% get center values of noisy sub-image, denoised sub-image and noise        
        [B,Cl,X,BCX] = BlockClipEdge(Zsub,Zsubapp,Zsubnoise,N);
        [Zsubc,Zsubappc,Zsubnoisec,~,~,~,BCXc] = GetBlockCenters...
                    (Zsub,Zsubapp,Zsubnoise,B,Cl,X,BCX,m,n,offsetmn,bsize);                                
        Zsubc = Zsubc.*BCXc; Zsubappc = Zsubappc.*BCXc; Zsubnoisec = Zsubnoisec.*BCXc;
        
        %%% clear unnecessary variables
        clear B Cl X BCX Zsub Zsubapp Zsubnoise;
        
        %%% get data associated with homogeneous segments
        [Set,Mu,Var,NumPix] = SetExpVarNumSeg(Zsubc,Zsubappc,Zsubnoisec,N,NumMin,stp);
        Set = Set(Set(:,1)~=1,:); 
             
        %%% estimate parameters of camera model from original sub-images and plot
        if j == 1
            % backup sample mean and sample variance
            MuOrg = Mu; VarOrg = Var;
            % estimate camera parameters
            ParOrg = Estabgam(Mu,Var,NumPix);
            % clear unnecessary variables
            clear Set SampMu SampVar NumPix Zsubc Zsubappc Zsubnoisec;
            % define fitted function
            ffun = @(mu) ParOrg(1)/(ParOrg(3)^2)*(mu.^(2-ParOrg(3))) + ...
                         ParOrg(2)/(ParOrg(3)^2)*(mu.^(2-2*ParOrg(3)));
        end
                
        %%% Decision from interpolated sub-images
        if j >= 2
            
            % estimate the sum of squared interpolation coefficients
            VarSig{i,j} = EstVarSig(Mu,Var,NumPix,ParOrg);
            
            % calibration if VarSig{i,j} >= 1
            if VarSig{i,j} >= 1
                VarSig{i,j} = 0.99999;
            end
            
            % repeat copies of elements of Mu wrt NumPix
            RepMu = repelem(Mu,NumPix);  
            
            % compute the log-likelihood ratio
            LLR{i,j} = sum( 0.5*(1-1/VarSig{i,j})*((Set(:,3)).^2)./ffun(RepMu)...
                          - 0.5*log(VarSig{i,j}) );
            
            % define shape and scale parameters under H0 and H1
            ShaPar = 0.5 * sum(NumPix); 
            ScaParH0 = 1/VarSig{i,j} - 1; ScaParH1 = 1 - VarSig{i,j};
            
            % compute A = 0.5 * log(1/VarSig{i,j}) * sum(NumPix) 
            A = 0.5 * log(1/VarSig{i,j}) * sum(NumPix);
            
            % compute decision threshold
            the{i,j} = A - gaminv(alp0,ShaPar,ScaParH0);
            
            % compute power func
            bet{i,j} = gamcdf(A - the{i,j}, ShaPar, ScaParH1);
            
            % decision for H1: 1 and H0: 0
            deci{i,j} = LLR{i,j} > the{i,j};
            
            % clear unnecessary variables
            clear Set SampMu SampVar Zsubc Zsubappc Zsubnoisec;
            
        end        
        
    end
    
end

%%%%%  Print results to the console

clc
fprintf(' When \x3B1%d = %f, then \n',0,alp0);
fprintf(' --------------------------------------\n');
fprintf('  + for the considered authentic image,\n');
fprintf('    - the values of LLR for sub-images 2, 3 and 4 are %f, %f and %f \n',LLR{3},LLR{5},LLR{7});
fprintf('    - the decision thesholds for sub-images 2, 3 and 4 are %d, %d and %d \n',the{3},the{5},the{7});
fprintf('    - the decisions for sub-images 2, 3 and 4 are H%d, H%d and H%d \n',deci{3},deci{5},deci{7});
fprintf(' --------------------------------------\n');
fprintf('  + for the considered resampled image,\n');
fprintf('    - the values of LLR for sub-images 2, 3 and 4 are %d, %d and %d \n',LLR{4},LLR{6},LLR{8});
fprintf('    - the decision thesholds for sub-images 2, 3 and 4 are %d, %d and %d \n',the{4},the{6},the{8});
fprintf('    - the decisions for sub-images 2, 3 and 4 are H%d, H%d and H%d \n',deci{4},deci{6},deci{8});
fprintf(' --------------------------------------\n');