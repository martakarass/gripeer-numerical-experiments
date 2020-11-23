
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------%
%         Authors:    Marta Karas <mkaras2@jhu.edu>       %
%         Date:       Aug 29, 2020                        %
%---------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: The toolbox 'penalized' should be installed before running griPEER function.
% See "README.md" for instructions.

function [] = run_scenario_1c_p528_CLUSTER(pathBase, AobsFname, nIter, nObs)

    % Project working directory
    %pathBase = '/Users/martakaras/Dropbox/_PROJECTS/OLD/LogisticPEER/gripeer-numerical-experiments';
    %AobsFname= 'A1_diffidx_4.txt'; 
    %nIter    = 2;
    %nObs     = 100;

    % Add path to software 
    addpath(strcat(pathBase, '/MATLAB/software'));
    %thedir = strrep(mfilename('fullpath'),[filesep,'install_penalized'],'');
    thedir = strcat(pathBase, '/MATLAB/penalized');
    addpath(thedir);
    filesep_here = '/'; 
    addpath([thedir,filesep_here,'models'])
    addpath([thedir,filesep_here,'internals'])
    addpath([thedir,filesep_here,'penalties'])
    addpath([thedir,filesep_here,'helpfiles'])
    addpath([thedir,filesep_here,'jss'])
    addpath([thedir,filesep_here,'jss/jsstable'])

    % Fixed simulation parameters 
    p        = 528;
    
    % Path to A, Aobs
    AFname = strcat(extractBetween(AobsFname, 1, 11), '1.txt'); 
    AFname = AFname{1};
    AFpath    = strcat(pathBase,'/data/sc_matrices_used_scenario_1c_p528/', AFname); 
    AobsFpath = strcat(pathBase,'/data/sc_matrices_used_scenario_1c_p528/', AobsFname); 

    % Path to save results
    resTableFname = strcat('results_scenario_1c_', 'nObs_', num2str(nObs), '_AFname_', extractBetween(AobsFname, 1, 12), '.txt'); 
    resTableFname = resTableFname{1}; 
    resTableFpath = strcat(pathBase, '/results/dir_results_scenario_1c/', resTableFname);

    % Read A, Aobs
    A    = importdata(AFpath);     
    Aobs = importdata(AobsFpath);  

    % Seed the random number generator
    rng(1);

    % Objects to store simulation results
    scMatrixNameVec        = []; 
    scMatrixPermParVec     = []; 
    nVec                   = [];
    iVec                   = [];
    bIdxVec                = []; 
    bTrueVec               = []; 
    bEstGripeerVec         = []; 
    bEstGripeer_c_ridgeVec = []; 
    lambQGripeerVec        = [];
    lambRGripeerVec        = [];
    lambQGripeer_c_ridgeVec= [];
    lambRGripeer_c_ridgeVec= [];
    timeGripeerVec         = [];
    timeGripeer_c_ridgeVec = [];

    % Data objects (fixed across simulations)
    Zsigma     = importdata(strcat(pathBase,'/data/Z_cov.txt'));
    % Zsigma correction for p=528
    Zsigma     = blkdiag(Zsigma,Zsigma,Zsigma,Zsigma,Zsigma,Zsigma,Zsigma,Zsigma); 
    Zmu        = zeros(p,1); 

    % Q matrix
    D                  = diag(sum(abs(A),1));
    Q                  = (D^(-0.5))*(D - A)*(D^(-0.5)); 

    % Qtilde, SIGMA
    [AA, BB, CC]       = svd(Q);
    Qrank              = rank(Q);
    BB2                = BB;
    BB2((Qrank+1):end,(Qrank+1):end) = 0.01 * eye(p-Qrank) * BB(Qrank, Qrank);
    BB2_diag = diag(BB2);
    Qtilde_pre = AA * diag(sqrt(1./BB2_diag)); 
    SIGMA = Qtilde_pre * Qtilde_pre';

    % Z, X
    Z         = mvnrnd(Zmu,Zsigma,nObs);
    Z         = zscore(Z);     
    Z         = Z/sqrt(nObs-1); 
    X         = ones(nObs,1);  

    % Qobs
    Dobs      = diag(sum(abs(Aobs),1));   
    Qobs      = (Dobs^(-0.5))*(Dobs - Aobs)*(Dobs^(-0.5));

    % Iterate over experiment runs  
    for ii = 1:nIter  
        fprintf('\n\niter #%d', ii);

        % Generate bTrue         
        bTrue     = (mvnrnd(zeros(p,1), SIGMA))';
        betaTrue  = 0;
        thetaTrue = X*betaTrue + Z*bTrue;
        BernPr    = exp(thetaTrue)./(1+exp(thetaTrue));
        y = binornd(1,BernPr);

        % griPEER family estimation parameters
        GRIPEER_stopCrit    = 1e-6; 
        GRIPEER_initLambs   = [0.0001; 0.00001];
        GRIPEER_initLambs_C = [0.0001];
        GRIPEER_maxIter     = 20; 
        GRIPEER_lambsGrid   = [0.001, 0.001, 0.001, 0.0001, 0.0001, 0.0001, 0.00001, 0.00001;   0.001, 0.0001, 0.00001, 0.001, 0.0001, 0.00001, 0.001, 0.0001];
        GRIPEER_lambsGrid_C = [0.001, 0.001, 0.001, 0.0001, 0.0001, 0.0001, 0.00001, 0.00001];

        % Estimation: griPEER
        fprintf('\nEstimation: griPEER...');
        ibEstGripeer  = NaN(p, 1); 
        ilambQGripeer = NaN; 
        ilambRGripeer = NaN; 
        tStartGripeer = tic;
        tEndGripeer   = NaN; 
        try
            OutGripeer    = griPEER(y, X, Z, Qobs, 'lambsGrid', GRIPEER_lambsGrid, 'maxIter', GRIPEER_maxIter, 'stopCrit', GRIPEER_stopCrit, 'initLambs', GRIPEER_initLambs);
            ibEstGripeer  = OutGripeer.b; 
            ilambQGripeer = OutGripeer.lambQ;
            ilambRGripeer = OutGripeer.lambR;
            tEndGripeer = round(toc(tStartGripeer)); 
            fprintf(strcat('\nEstimation time griPEER: ', num2str(tEndGripeer)));
        catch err
            fprintf('\ngriPEER estimation failed. Returning NaN(p, 1).');
        end 

        % Estimation: griPEER_c_ridge
        fprintf('\nEstimation: griPEER_c_ridge...');
        ibEstGripeer_c_ridge  = NaN(p, 1); 
        ilambQGripeer_c_ridge = NaN; 
        ilambRGripeer_c_ridge = NaN;
        tStartGripeer_c_ridge = tic;
        tEndGripeer_c_ridge   = NaN; 
        try
            OutGripeer_c_ridge    = griPEER(y, X, Z, eye(p), 'lambsGrid', GRIPEER_lambsGrid_C, 'maxIter', GRIPEER_maxIter, 'stopCrit', GRIPEER_stopCrit, 'lambdaR', 0, 'initLambs', GRIPEER_initLambs_C);
            ibEstGripeer_c_ridge  = OutGripeer_c_ridge.b; 
            ilambQGripeer_c_ridge = OutGripeer_c_ridge.lambQ;
            ilambRGripeer_c_ridge = OutGripeer_c_ridge.lambR;
            tEndGripeer_c_ridge   = round(toc(tStartGripeer_c_ridge)); 
            fprintf(strcat('\nEstimation time griPEER_c_ridge: ', num2str(tEndGripeer_c_ridge)));
        catch err
            warning('\ngriPEER_c_ridge estimation failed. Returning NaN(p, 1).');
        end

        % Store simulation results - current state
        scMatrixNameVec        = [scMatrixNameVec; cellstr(repmat(AobsFname,p,1))]; 
        scMatrixPermParVec     = [scMatrixPermParVec; repmat(-1,p,1)]; 
        nVec                   = [nVec; repmat(nObs,p,1)];
        iVec                   = [iVec; repmat(ii,p,1)]; 
        bIdxVec                = [bIdxVec; (1:1:p)']; 
        bTrueVec               = [bTrueVec; bTrue]; 
        bEstGripeerVec         = [bEstGripeerVec; ibEstGripeer];
        bEstGripeer_c_ridgeVec = [bEstGripeer_c_ridgeVec; ibEstGripeer_c_ridge];
        lambQGripeerVec        = [lambQGripeerVec; repmat(ilambQGripeer,p,1)];
        lambRGripeerVec        = [lambRGripeerVec; repmat(ilambRGripeer,p,1)]; 
        lambQGripeer_c_ridgeVec= [lambQGripeer_c_ridgeVec; repmat(ilambQGripeer_c_ridge,p,1)];
        lambRGripeer_c_ridgeVec= [lambRGripeer_c_ridgeVec; repmat(ilambRGripeer_c_ridge,p,1)]; 
        timeGripeerVec         = [timeGripeerVec; repmat(tEndGripeer,p,1)];
        timeGripeer_c_ridgeVec = [timeGripeer_c_ridgeVec; repmat(tEndGripeer_c_ridge,p,1)];     

        % Wrapping up the results
        resTable = table(nominal(scMatrixNameVec),scMatrixPermParVec,nVec,iVec,bIdxVec,bTrueVec,...
                         bEstGripeerVec,bEstGripeer_c_ridgeVec,...
                         lambQGripeerVec,lambRGripeerVec,lambQGripeer_c_ridgeVec,lambRGripeer_c_ridgeVec,...
                         timeGripeerVec,timeGripeer_c_ridgeVec,...
                         'VariableNames',{'SC_mat' 'SC_mat_perm_val' 'n' 'i' 'b_idx' 'b_true' 'b_est_griPEER' 'b_est_griPEER_c_ridge' 'lambQ_Gripeer' 'lambR_Gripeer' 'lambQ_Gripeer_c_ridge' 'lambR_Gripeer_c_ridge' 'exectime_Gripeer' 'exectime_Gripeer_c_ridge' });

        % Save results to file
        writetable(resTable, resTableFpath);

    end

end
          
 
 
 
