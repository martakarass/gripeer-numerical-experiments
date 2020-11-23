function  out =  griPEER(y, X, Z, Q, varargin)

%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
% This function solves the griPEER optimization problem of the form:
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% argmin_{Beta, b}   { -2loglik(X*Beta + Z*b) + 
%           tau||Beta(2:end)||^2 + lambdaQ*b'*Q*b + lambdaR||b||^2 +  }   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% where:
% loglik is the loglikelihood function of member of one-parameter 
% exponential family of distributions
%-------------------------------------------
% y is n-dimensional vector of observations
%-------------------------------------------
% X is n by m design matrix of confounding variables (with ones in first
% colum corresponding to the intercept)
%-------------------------------------------
% Beta is m-dimensional vector of interest
%-------------------------------------------
% Z is n by p design matrix of penalized variables
%-------------------------------------------
% b is p-dimensional vector of interest
%-------------------------------------------
% Q is p by b symmetric, positive-semidefiniete matrix, which imposes the
% penalty on coefficients in vector b. If Q is not of that form but has
% only non-negative entries and zeros on diagonal, then Q is labeled as
% adjacency matrix and normalized Laplacian is used instead.
%-------------------------------------------
% tau is an optional, user-specified parameter providing additional penalty
% on the confounding variables (set as zero by default)
%-------------------------------------------
% lambda1 and lambda2 are positive regularization parameters which are 
% automaticly adjusted by taking the advantage of the connection of 
% considered problem with generlized linear mixed model formulation.
%-------------------------------------------

%%%%%%%%%%%%%%%%%%         OPTIONAL INPUTS           %%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------
% 'nboot' --- number of bootstrap samples produced to calcuate confidence
% interval. By default 'nboot' = 500;
%-------------------------------------------
% 'alpha' --- statistical significance level.By default 'alpha' = 0.05
%-------------------------------------------
% 'UseParallel' --- logical value indicating if parallel computatiation
% should be used when bootstrap is employed. By default 'UseParallel' = true
%-------------------------------------------
% 'ciType' --- default set as 'asymptotic', which generates only
% analitically derived confidence interval. To get also the bootstrap
% confidence interval set this as 'both' (this may signifficantly increase
% the time of computing).
%-------------------------------------------
% 'bootType' --- type of bootstrap confidence interval. Check the 
% documentation of 'bootci' function for more details 
%-------------------------------------------
% 'tau' --- the regulatization parameter used in the penalty applied to
% Beta(2:end). It is set as 0 by default but can be changed to any positive
% number.
%-------------------------------------------
% 'family' --- 'binomial' or 'poisson' family can be chosen. The first
% option is default one.
%-------------------------------------------
% 'lambdaR' --- the user may select if this parameter should be automatically
% selected (the default option) or if it should be fixed at some given level.
% If the second option is preferred (and 'lambdaR' is set as fixed number),
% then only 'lambdaQ' is adjusted automatically.
%-------------------------------------------
% 'lambsGrid' --- this should be a 2-by-k matrix, where k is a number of
% considered starting points. The function will check every column as a
% starting point and then will select the best option.
%-------------------------------------------
% 'maxIter' --- the maximal number of iteration for lambdas selection. Set
% as 30 by default.
%-------------------------------------------
% 'stopCrit' --- lambdas convergence tollerance, the default value is 1e-6.
%-------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%         OUTPUTS         %%%%%%%%%%%%%%%%%%%%%%%%
% 'beta'   --- estimate of Beta
% 'b'      --- estimate of b
% 'CI'     --- confidence interval for the [Beta;b] derived analitically 
% 'lambQ'  --- optimal value of lambdaQ
% 'lambR'  --- optimal value of lambdaR
% 'CIboot' --- confidence interval for the [Beta;b] based on bootstrap
%              (only if 'ciType' is set as 'both')
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------
%-------------------------------------------------------------------------

% The griPEER toolbox is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% The griPEER toolbox is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-------------------------------------------
%         Authors:    Damian Brzyski
%         Date:       July 27, 2017
%         email:      damian.brzyski@pwr.edu.pl
%-------------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Checks
if exist('glm_logistic', 'file')==0
cd    error(['The toolbox ',char(39), 'penalized',char(39),' ', 'should be installed before running griPEER function. Go to https://www.jstatsoft.org/article/view/v072i06, download the toolbox, change Current Folder to the directory containing the toolbox and type ',char(39), 'install_penalized',char(39),'.'])
end
if any(abs(mean(Z,1)>1e-10))
   error('Columns in design matrix Z should be centered to 0 means') 
end
%-----------------
if size(X,2) ==0
    error('It should be at least one column in design matrix X')
end
%-----------------
if any(abs(mean(X(:,2:end),1) > 1e-10))
    error('Besides first column (corresponding to the intercept), all columns in design matrix X should be centered to 0 means') 
end
%-----------------
if max(max(abs(bsxfun( @minus, sqrt(sum([X(:,2:end), Z].^2,1)), sqrt(sum([X(:,2:end), Z].^2,1))'))))>1e-10
    error('Columns in X (beyond the first column) and columns in Z should have the same Euclidean norms')
end
%-----------------
if unique(X(:,1))~=1
    error('First column of design matrix X should represent the intercept and consist of ones')
end

%% Warning turn off
warn_id   = 'stats:classreg:regr:lmeutils:StandardGeneralizedLinearMixedModel:BadFinalPLSolution';
warning('off',warn_id)

%% Additional parameters (AP)
AP               = inputParser;
defaultTau       = 0;
defaultNboot     = 500;
defaultAlpha     = 0.05;
defaultParallel  = false;
defaultType      = 'per';
defaultCItype    = 'asymptotic';
defaultFamily    = 'binomial';
defaultLambdaR   = 'auto';
defaultInLambs   = [0.0001; 0.00001];
defaultGridd     = [0.0001; 0.00001];
defaultMaxIter   = 30;
defaultStopCrit  = 1e-6;

%------------------------------------
addRequired(AP, 'y');
addRequired(AP, 'X');
addRequired(AP, 'Z');
addRequired(AP, 'Q');
addOptional(AP, 'tau', defaultTau, @(x)(x>=0) );
addOptional(AP, 'nboot', defaultNboot, @isnumeric);
addOptional(AP, 'alpha', defaultAlpha, @isnumeric);
addOptional(AP, 'UseParallel', defaultParallel, @islogical);
addOptional(AP, 'bootType', defaultType, @(x) any(validatestring(x,{'norm', 'per', 'cper', 'bca', 'stud'} ) ) );
addOptional(AP, 'ciType', defaultCItype, @(x) any(validatestring(x,{'both', 'asymptotic'} ) ) );
addOptional(AP, 'family', defaultFamily, @(x) any(validatestring(x,{'binomial', 'poisson'} ) ) );
addOptional(AP, 'lambdaR', defaultLambdaR, @(x) (x>=0) );
addOptional(AP, 'initLambs', defaultInLambs, @(x) all(x>0) );
addOptional(AP, 'lambsGrid', defaultGridd, @(x) all(all(x>0)) );
addOptional(AP, 'maxIter', defaultMaxIter, @(x) (x>0) );
addOptional(AP, 'stopCrit', defaultStopCrit, @(x) (x>0) );

%-------------------------------------
 parse(AP, y, X, Z, Q, varargin{:}) %-
%-------------------------------------

tau         = AP.Results.tau;
nboot       = AP.Results.nboot;
alpha       = AP.Results.alpha; 
type        = AP.Results.bootType;
ciType      = AP.Results.ciType;
family      = AP.Results.family;
lambdaR     = AP.Results.lambdaR;
intitLambs  = AP.Results.initLambs;
lambsGrid   = AP.Results.lambsGrid;
maxIter     = AP.Results.maxIter;
stopCrit    = AP.Results.stopCrit;
options     = statset('UseParallel', AP.Results.UseParallel);  

%% More checks
%=============================================
if and(strcmp(lambdaR, 'auto'), ~isequal(size(intitLambs),[2,1])   )
    error('intitLambs should be a vector of length 2, since lambdaQ and lambdaR are chosen to be automatically adjusted.')
end
%-----------------
if and(  and(~strcmp(lambdaR, 'auto'), length(intitLambs)>1),    ~isequal(defaultInLambs, intitLambs)  )
    warning('only first entry of intitLambs was used,  since lambdaR was provided.')
end
%-----------------
if and(strcmp(lambdaR, 'auto'), size(lambsGrid,1)<2)
    error('It should be two rows in lambsGrid matrix, if lambdaQ and lambdaR are automatically adjusted.')
end
%-----------------
if strcmp(lambdaR, 'auto')
    gridd    = lambsGrid(1:2,:);
else
    gridd    = unique(lambsGrid(1,:));
end
%-----------------
if and(and(size(lambsGrid,1)>1, ~strcmp(lambdaR, 'auto')),  ~isequal(defaultGridd, lambsGrid) )        
    warning('Only one row of lambsGrid was taken, since lambdaR was provided.')
end
%-----------------
if and(size(lambsGrid,1)>2, strcmp(lambdaR, 'auto'))               
    warning('Only two first rows of lambsGrid were taken (first gives a grid for lambdaQ and the second oone gives a grid for lambdaR)')
end

%% Objects
Xr            =  X(:, 2:end);
[n,p]         =  size(Z);
m             =  size(X,2);                               
colsNorm      =  mean(sqrt(sum([X(:,2:end), Z].^2,1))); 
max_lam_iter  =  maxIter;

%------------ Options for optimization software ------------------
opts                          = optimoptions('fsolve','Display','off');
opts.Algorithm                = 'trust-region';
opts.FunctionTolerance        = 1e-10;
opts.OptimalityTolerance      = 1e-10;
opts.SpecifyObjectiveGradient = true;

%% Checks for the situation with lambdaR defined by the user as 0 
%-----------------------------------------------
if and(lambdaR==0, cond(Q) > 1e8)
   error('When lambdaR is fixed as zero, conditional number of Q is assumed to be not larger than 1e8. Consider the modiffication of Q by adding the identity matrix with small multiplicative constant.')
end
%-----------------
if m == 1 && tau > 0
    error('tau>0 is not used, when matrix X contains only one column of ones corresponding to intercept')
end

%% Family of distrubution
if strcmp(family, 'binomial')
    psiPrim                =  @(x) exp(x)./(1 + exp(x));
    psiBis                 =  @(x) exp(x)./( (1+exp(x)).^2);
    glm_chosen_family      =  @(y, X) glm_logistic(y, X, 'nointercept');
else
    psiPrim                =  @(x) exp(x);
    psiBis                 =  @(x) exp(x);
    glm_chosen_family      =  @(y, X) glm_poisson(y, X, 'nointercept');
end

%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ______________________________________________________
%|                                                      |
%|---------  Condition: tau == 0 versus tau >0  --------|
%|______________________________________________________|
%

lamQ_prev  = intitLambs(1);
if strcmp(lambdaR, 'auto')
    lamR_prev  = intitLambs(2);
else
    lamR_prev  = lambdaR;
end
stopp      = 0;
countt     = 1;
Lambdas    = zeros(2, max_lam_iter);
bs         = zeros(p, max_lam_iter);
Gradients  = zeros(2, max_lam_iter);
Hessians   = zeros(4, max_lam_iter);
Ispos      = zeros(1, max_lam_iter);

%--------- CONDITION: tau == 0
if tau == 0   
    %-------------------------------------------------  
    % WHILE: START
    while  stopp==0
    Qt_prev       =  lamQ_prev*Q + lamR_prev*eye(p);
    L             =  chol(Qt_prev);
    ZL            =  Z* (L^(-1));  
    model         =  glm_chosen_family(y, [X, ZL]);
    fit           =  penalized(model, @p_ridge, 'penaltywt', [zeros(m,1); ones(p,1)], 'lambda', 1/n, 'standardize', false);   
    Beta          =  fit.beta(1:m); 
    b             =  L^(-1)* (fit.beta((m+1):end));
    theta         =  X*Beta + Z*b;
    W_k           =  diag( psiBis(theta).^(-1) ); 
    W_k_inv       =  diag( psiBis(theta) );
    y_k           =  W_k*(y - psiPrim(theta) ) + theta;
    P             =  eye(n) - X*(X'*W_k_inv*X)^(-1)*X'*W_k_inv;
    PZ            =  P*Z;   
    yt_k          =  P*y_k;   
    Omega_k       =  PZ'* W_k_inv *PZ;
    q_k           =  PZ'* W_k_inv *yt_k;
    
    % finding the update of lambda
    %-------------------------------------------------
    %_________________________________________________
    fun           =  @(x)log(  det( (abs(x(1))*Q + abs(x(2))*eye(p) + Omega_k)*(abs(x(1))*Q + abs(x(2))*eye(p))^(-1) )  ) - q_k'*( abs(x(1))*Q + abs(x(2))*eye(p) + Omega_k )^(-1)*q_k;
    val_prop      =  zeros(1,size(gridd, 2));
    lambda_prop   =  zeros(2,size(gridd, 2));
    
    % looking at the border
    gridd2        =  unique(gridd(1,:));
    val_prop2     =  zeros(1,size(gridd2, 2));
    lambda_prop2  =  zeros(2,size(gridd2, 2));
    
    % looking over the grid of starting points
    if strcmp(lambdaR, 'auto')
        if options.UseParallel
            parfor gg = 1:size(gridd, 2)            
                lambda_prop(:, gg)  =  abs(fsolve(@(x) fbnd_o(Q, Omega_k, q_k, p, x), gridd(:,gg), opts));
                val_prop(gg)        =  fun(lambda_prop(:, gg));
            end
        else
            for gg = 1:size(gridd, 2)            
                lambda_prop(:, gg)  =  abs(fsolve(@(x) fbnd_o(Q, Omega_k, q_k, p, x), gridd(:,gg), opts));
                val_prop(gg)        =  fun(lambda_prop(:, gg));
            end
        end
        if options.UseParallel
            parfor gg = 1:length(gridd2)            
                lambda_prop2(:, gg)  =  [0; abs(fsolve(@(x) fbnd_zero_lambQ(Omega_k, q_k, p, x), gridd2(:,gg), opts))];
                val_prop2(gg)        =  fun(lambda_prop2(:, gg));
            end
        else
            for gg = 1:length(gridd2)            
                lambda_prop2(:, gg)  =  [0; abs(fsolve(@(x) fbnd_zero_lambQ(Omega_k, q_k, p, x), gridd2(:,gg), opts))];
                val_prop2(gg)        =  fun(lambda_prop2(:, gg));
            end
        end
        lambda_prop3       =  [lambda_prop, lambda_prop2];
        val_prop3          =  [val_prop, val_prop2];
    else
        if options.UseParallel
            parfor gg = 1:size(gridd, 2)
                lambda_prop(:, gg)  =  [abs(fsolve(@(x) fbnd_o_fix_lambR(Q, Omega_k, q_k, p, lamR_prev, x), gridd(:,gg), opts)); lamR_prev];
                val_prop(gg)        =  fun(lambda_prop(:, gg));
            end
        else
            for gg = 1:size(gridd, 2)
                lambda_prop(:, gg)  =  [abs(fsolve(@(x) fbnd_o_fix_lambR(Q, Omega_k, q_k, p, lamR_prev, x), gridd(:,gg), opts)); lamR_prev];
                val_prop(gg)        =  fun(lambda_prop(:, gg));
            end
        end
        lambda_prop3       =  lambda_prop;
        val_prop3          =  val_prop;
    end
    [~, minIdx]        =  min(val_prop3);
    lambda_k           =  lambda_prop3(:, minIdx);
    Lambdas(:, countt) =  lambda_k;
    estim              =  [Beta; b];
    bs(:, countt)      =  estim((m+1):end);
    %------------------
    if strcmp(lambdaR, 'auto')
        [F, H] = fbnd_o(Q, Omega_k, q_k, p, lambda_k);
    else
        [F, H] = fbnd_o_fix_lambR(Q, Omega_k, q_k, p, lamR_prev, lambda_k(1));
    end
    Gradients(:, countt) =  F;
    Hessians(:, countt)  =  H(:);
    %------------------
    if and( H(1,1) >= 0,  det(H)>=0 )
        Ispos(countt) = 1;
    end    
    %-------------------------------------------------
    
    % checking stop conditions
    lambda_prev = [lamQ_prev; lamR_prev];
    if (norm(lambda_k - lambda_prev)/norm(lambda_prev)) < stopCrit || countt >= max_lam_iter
        stopp = 1;
    else
        countt = countt + 1;
    end
    %-------------------------------------------------
    lamQ_prev     =  lambda_k(1);
    lamR_prev     =  lambda_k(2);
    end
    
    % final estimate
    lamQ_final    =  lamQ_prev;
    lamR_final    =  lamR_prev;
    Qt_final      =  lamQ_final*Q + lamR_final*eye(p);
    L_final       =  chol(Qt_final);
    L_final_inv   =  L_final^(-1);
    ZL            =  Z*L_final_inv;
    model         =  glm_chosen_family(y, [X, ZL]);
    fit           =  penalized(model, @p_ridge, 'penaltywt', [zeros(m,1); ones(p,1)], 'lambda', 1/n, 'standardize', false);   
    Estimate      =  [fit.beta( 1:m ); L_final_inv * (fit.beta( (m+1):end) )];
    Beta_est      =  Estimate( 1:m ); 
    b_est         =  Estimate( (m+1):end );
    
    % confidence interval
    theta_final   =  [X,Z]*Estimate;
    Psi           =  diag(psiBis(theta_final));
    XZ            =  [X, Z];
    Qt_final_ext  =  blkdiag(zeros(m,m), Qt_final);
    estimVar      =  (XZ'*Psi*XZ + Qt_final_ext)^(-1)*XZ'*Psi*XZ*(XZ'*Psi*XZ + Qt_final_ext)^(-1);
    deltaa        =  norminv(1-alpha/2) * sqrt(diag(estimVar));
    estimBand     =  [ [Beta_est; b_est] - deltaa, [Beta_est; b_est] + deltaa ];
    
    % bootstrap confidence interval (optionally)
    if strcmp(ciType, 'both')
        DATASET       = [y, X, Z];
        if strcmp(family, 'binomial')
            boot          = @(Dataset) boot_griPEER_binomial(L_final_inv, m, colsNorm, Dataset);
        else
            boot          = @(Dataset) boot_griPEER_poisson(L_final_inv, m, colsNorm, Dataset);
        end
        CI_boot       = ( bootci(nboot,{boot, DATASET}, 'Options', options, 'type', type, 'alpha', alpha) )'; 
    end
  
%--------- CONDITION: tau > 0
else
    %-------------------------------------------------       
    % while loop
    while  stopp==0
    Qt_prev       =  [[tau*eye(m-1), zeros(m-1,p)]; [zeros(p, m-1), lamQ_prev*Q + lamR_prev*eye(p)]];
    L             =  chol(Qt_prev);
    XZL           =  [Xr, Z]*( L^(-1) );
    model         =  glm_chosen_family(y, [ones(n,1), XZL]);
    fit           =  penalized(model, @p_ridge, 'penaltywt', [zeros(1,1); ones(p+m-1,1)], 'lambda', 1/n, 'standardize', false);
    int           =  fit.beta(1); 
    Beta_b        =  L^(-1)* (fit.beta(2:end));
    theta         =  int + [Xr, Z]*Beta_b;
    W_k           =  ( diag(psiBis(theta)) )^(-1);  
    W_k_inv       =  W_k^(-1);
    y_k           =  W_k*(y - psiPrim(theta)) + theta;
    P             =  eye(n) - ones(n,1) * (ones(n,1)'*W_k_inv*ones(n,1))^(-1) * ones(n,1)'*W_k_inv;
    PZ            =  P*Z;   
    yt_k          =  P*y_k;
    Omega_k       =  PZ'*(W_k + Xr*Xr'/tau)^(-1)*PZ;
    q_k           =  PZ'*(W_k + Xr*Xr'/tau)^(-1)*yt_k;
    
    % finding the update of lambda
    %-------------------------------------------------
    %_________________________________________________

    fun          =  @(x)log(  det( (abs(x(1))*Q + abs(x(2))*eye(p) + Omega_k)*(abs(x(1))*Q + abs(x(2))*eye(p))^(-1) )  ) - q_k'*( abs(x(1))*Q + abs(x(2))*eye(p) + Omega_k )^(-1)*q_k;
    lambda_prop  =  zeros(2,size(gridd, 2));
    val_prop     =  zeros(1,size(gridd, 2));
    
    % looking at the border
    gridd2        =  unique(gridd(1,:));
    val_prop2     =  zeros(1,size(gridd2, 2));
    lambda_prop2  =  zeros(2,size(gridd2, 2));
    
    % looking over the grid of starting points
    if strcmp(lambdaR, 'auto')
        if options.UseParallel
            parfor gg = 1:size(gridd, 2)
                lambda_prop(:, gg)  =  abs(fsolve(@(x) fbnd_o(Q, Omega_k, q_k, p, x), gridd(:,gg), opts));
                val_prop(gg)        =  fun(lambda_prop(:, gg));
            end
        else
           for gg = 1:size(gridd, 2)
                lambda_prop(:, gg)  =  abs(fsolve(@(x) fbnd_o(Q, Omega_k, q_k, p, x), gridd(:,gg), opts));
                val_prop(gg)        =  fun(lambda_prop(:, gg));
            end 
        end
        if options.UseParallel
            parfor gg = 1:length(gridd2)            
                lambda_prop2(:, gg)  =  [0; abs(fsolve(@(x) fbnd_zero_lambQ(Omega_k, q_k, p, x), gridd2(:,gg), opts))];
                val_prop2(gg)        =  fun(lambda_prop2(:, gg));
            end
        else
            for gg = 1:length(gridd2)            
                lambda_prop2(:, gg)  =  [0; abs(fsolve(@(x) fbnd_zero_lambQ(Omega_k, q_k, p, x), gridd2(:,gg), opts))];
                val_prop2(gg)        =  fun(lambda_prop2(:, gg));
            end
        end
        lambda_prop3       =  [lambda_prop, lambda_prop2];
        val_prop3          =  [val_prop, val_prop2];
    else
        if options.UseParallel
            parfor gg = 1:size(gridd, 2)
                lambda_prop(:, gg)  =  [abs(fsolve(@(x) fbnd_o_fix_lambR(Q, Omega_k, q_k, p, lamR_prev, x), gridd(:,gg), opts)); lamR_prev];
                val_prop(gg)        =  fun(lambda_prop(:, gg));
            end
        else
            for gg = 1:size(gridd, 2)
                lambda_prop(:, gg)  =  [abs(fsolve(@(x) fbnd_o_fix_lambR(Q, Omega_k, q_k, p, lamR_prev, x), gridd(:,gg), opts)); lamR_prev];
                val_prop(gg)        =  fun(lambda_prop(:, gg));
            end
        end
        lambda_prop3       =  lambda_prop;
        val_prop3          =  val_prop;    
    end
    [~, minIdx]        =  min(val_prop3);
    lambda_k           =  lambda_prop3(:, minIdx);
    Lambdas(:, countt) =  lambda_k;
    estim              =  [int; Beta_b];
    bs(:, countt)      =  estim((m+1):end);
        %------------------
    if strcmp(lambdaR, 'auto')
        [F, H] = fbnd_o(Q, Omega_k, q_k, p, lambda_k);
    else
        [F, H] = fbnd_o_fix_lambR(Q, Omega_k, q_k, p, lamR_prev, lambda_k(1));
    end
    Gradients(:, countt) =  F;
    Hessians(:, countt)  =  H(:);
    %------------------
    if and( H(1,1) >= 0,  det(H)>=0 )
        Ispos(countt) = 1;
    end   
    %------------------
    %_________________________________________________
    
    %-------------------------------------------------
    
    % checking stop conditions
    lambda_prev = [lamQ_prev; lamR_prev];
    if (norm(lambda_k - lambda_prev)/norm(lambda_prev)) < stopCrit || countt >= max_lam_iter
        stopp = 1;
    else
        countt = countt + 1;
    end
    %-------------------------------------------------
    lamQ_prev     =  lambda_k(1);
    lamR_prev     =  lambda_k(2);
    end
    
    % final estimate
    lamQ_final    =  lamQ_prev;
    lamR_final    =  lamR_prev;
    Qt_final      =  [  [ tau*eye(m-1), zeros(m-1,p) ]; [ zeros(p, m-1), lamQ_final*Q + lamR_final*eye(p) ]  ];
    L_final       =  chol(Qt_final);
    L_final_inv   =  L_final^(-1);
    XZL           =  [Xr, Z] * L_final_inv;
    model         =  glm_chosen_family(y, [ones(n,1), XZL]);
    fit           =  penalized(model, @p_ridge, 'penaltywt', [zeros(1,1); ones(p+m-1,1)], 'lambda', 1/n, 'standardize', false);
    int_est       =  fit.beta(1); 
    Beta_ni_b_est =  L_final_inv * (fit.beta(2:end));
    Beta_est      =  [int_est; Beta_ni_b_est( 1:(m-1) )];
    b_est         =  Beta_ni_b_est( m:end );
    
    % confidence interval
    theta_final   =  [X,Z]*[int_est; Beta_ni_b_est];
    Psi           =  diag( psiBis(theta_final) );
    XZ            =  [X, Z];
    Qt_final_ext  =  blkdiag(0, Qt_final);
    estimVar      =  (XZ'*Psi*XZ + Qt_final_ext)^(-1)*XZ'*Psi*XZ*(XZ'*Psi*XZ + Qt_final_ext)^(-1);
    deltaa        =  norminv(1-alpha/2) * sqrt(diag(estimVar));
    estimBand     =  [ [Beta_est; b_est] - deltaa, [Beta_est; b_est] + deltaa ];
    
    % bootstrap confidence interval (optionally)
    if strcmp(ciType, 'both')
        DATASET       = [y, Xr, Z];
        if strcmp(family, 'binomial')
            boot          = @(Dataset) boot_griPEER_tau_binomial(L_final_inv, colsNorm, Dataset);
        else
            boot          = @(Dataset) boot_griPEER_tau_poisson(L_final_inv, colsNorm, Dataset);
        end
        CI_boot       = ( bootci(nboot,{boot, DATASET}, 'Options', options, 'type', type, 'alpha', alpha) )'; 
    end
    
    % WHILE: STOP
end 

%%
% ______________________________________________________
%|                                                      |
%|--------------  Output of the function  --------------|
%|______________________________________________________|
%
%% OUTPUTS
% standard
out          = struct;
out.beta     = Beta_est;
out.b        = b_est;
out.band     = estimBand;
out.lambQ    = lamQ_final;
out.lambR    = lamR_final;

% additional
if strcmp(ciType, 'both')
    out.ciboot   = CI_boot;
end

%______________________________________________________________

end
% ______________________________________________________
%|                                                      |
%|----------------- SUBROUTINES ------------------------|
%|______________________________________________________|
%

%-------------------------------------------------------------

function [F, H] = fbnd_o(Q, Omega, q, p, x)
%==========================================
ax    =  abs(x);
D0    =  ( ax(1)*Q + ax(2)*eye(p) )^(-1);
D     =  ( ax(1)*Q + ax(2)*eye(p) + Omega )^(-1);
D0Q   =  D0*Q;
DQ    =  D*Q;
Dq    =  D*q;

%-------------------------------------------
F     =  zeros(2,1);
F(1)  =  sign(x(1))* ( trace(   DQ - D0Q  )   +   Dq'*Q*Dq );
F(2)  =  sign(x(2))* (  trace(   D - D0  )    +    Dq'*Dq  );

%-------------------------------------------
H1(1, 1)  =  - trace(  DQ^2 - D0Q^2  );
H1(1, 2)  =  - sign(x(1))*sign(x(2))* trace(  (D^2 - D0^2)*Q  );
H1(2, 1)  =    H1(1, 2);
H1(2, 2)  =  - trace(  D^2 - D0^2  );

H2(1, 1) =  -2*q'*DQ^2*D*q;
H2(1, 2) =  -sign(x(1))*sign(x(2))*( Dq'*DQ*Dq + Dq'*Q*D*Dq );
H2(2, 1) =   H2(1, 2);
H2(2, 2) =  -2*Dq'*D*Dq;

H(1,1)   =  H1(1, 1) + H2(1, 1);
H(2,1)   =  H1(2, 1) + H2(2, 1);
H(1,2)   =  H(2,1);
H(2,2)   =  H1(2, 2) + H2(2, 2);

%==========================================
end

%-------------------------------------------------------------

function [F, H] = fbnd_o_fix_lambR(Q, Omega, q, p, lambR, x)
%==========================================
ax    =  [abs(x); lambR];
D0    =  ( ax(1)*Q + ax(2)*eye(p) )^(-1);
D     =  ( ax(1)*Q + ax(2)*eye(p) + Omega )^(-1);
D0Q   =  D0*Q;
DQ    =  D*Q;
Dq    =  D*q;

%-------------------------------------------
F     =  sign(x(1))* ( trace(   DQ - D0Q  )   +   Dq'*Q*Dq );

%-------------------------------------------
H1  =  - trace(  DQ^2 - D0Q^2  );
H2  =  -2*q'*DQ^2*D*q;
H   =  H1 + H2;
%==========================================
end

%-------------------------------------------------------------

function [F, H] = fbnd_zero_lambQ(Omega, q, p, x)
%==========================================
ax    =  abs(x);
D0    =  (ax*eye(p) )^(-1);
D     =  (ax*eye(p) + Omega )^(-1);
Dq    =  D*q;

%-------------------------------------------
F     =  sign(x)* (  trace(   D - D0  )    +    Dq'*Dq  );

%-------------------------------------------
H1  =  -trace(  D^2 - D0^2  );
H2  =  -2*Dq'*D*Dq;
H   =  H1 + H2;

%==========================================
end
%-------------------------------------------------------------

function estimate = boot_griPEER_binomial(L_final_inv, m, colsNorm, Dataset)
y              =   Dataset( :, 1 );
n              =   length(y);
X              =   Dataset( :, 2:(m+1)   );
Z              =   Dataset( :, (m+2):end );
Z              =   zscore(Z)*colsNorm/sqrt(n-1);
X(:,2:end)     =   zscore(X(:,2:end))*colsNorm/sqrt(n-1);
p              =   size(Z,2);

%------------------   optimization problem   ------------------
ZL             =  Z * L_final_inv;
model          =  glm_logistic(y, [X, ZL], 'nointercept');
fit            =  penalized(model, @p_ridge, 'penaltywt', [zeros(m,1); ones(p,1)], 'lambda', 1/n, 'standardize', false); 
estimate       =  [fit.beta( 1:m ); L_final_inv * (fit.beta( (m+1):end) )];

end

%-------------------------------------------------------------

function estimate = boot_griPEER_poisson(L_final_inv, m, colsNorm, Dataset)
y              =   Dataset( :, 1 );
n              =   length(y);
X              =   Dataset( :, 2:(m+1)   );
Z              =   Dataset( :, (m+2):end );
Z              =   zscore(Z)*colsNorm/sqrt(n-1);
X(:,2:end)     =   zscore(X(:,2:end))*colsNorm/sqrt(n-1);
p              =   size(Z,2);

%------------------   optimization problem   ------------------
ZL             =  Z * L_final_inv;
model          =  glm_poisson(y, [X, ZL], 'nointercept');
fit            =  penalized(model, @p_ridge, 'penaltywt', [zeros(m,1); ones(p,1)], 'lambda', 1/n, 'standardize', false); 
estimate       =  [fit.beta( 1:m ); L_final_inv * (fit.beta( (m+1):end) )];

end

%-------------------------------------------------------------

function estimate = boot_griPEER_tau_binomial(L_final_inv, colsNorm, Dataset)
y              =   Dataset( :, 1 );
n              =   length(y);
XZ             =   Dataset( :, 2:end );
XZ             =   zscore(XZ)*colsNorm/sqrt(n-1);
p              =   size(XZ,2);

%------------------   optimization problem   ------------------
XZL            =  XZ*L_final_inv;
model          =  glm_logistic(y, [ones(n,1), XZL], 'nointercept');
fit            =  penalized(model, @p_ridge, 'penaltywt', [zeros(1,1); ones(p,1)], 'lambda', 1/n, 'standardize', false); 
estimate       =  [ fit.beta(1); L_final_inv * (fit.beta( 2:end)) ];

end
%-------------------------------------------------------------

function estimate = boot_griPEER_tau_poisson(L_final_inv, colsNorm, Dataset)
y              =   Dataset( :, 1 );
n              =   length(y);
XZ             =   Dataset( :, 2:end );
XZ             =   zscore(XZ)*colsNorm/sqrt(n-1);
p              =   size(XZ,2);

%------------------   optimization problem   ------------------
XZL            =  XZ*L_final_inv;
model          =  glm_poisson(y, [ones(n,1), XZL], 'nointercept');
fit            =  penalized(model, @p_ridge, 'penaltywt', [zeros(1,1); ones(p,1)], 'lambda', 1/n, 'standardize', false); 
estimate       =  [ fit.beta(1); L_final_inv * (fit.beta( 2:end)) ];

end
