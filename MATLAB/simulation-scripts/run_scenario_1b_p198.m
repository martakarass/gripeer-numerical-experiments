
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------%
%         Authors:    Marta Karas <mkaras2@jhu.edu>       %
%         Date:       Oct 13, 2019                        %
%         Most recent update: Aug 29, 2020                %
%---------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: The toolbox 'penalized' should be installed before running griPEER function.
% See "README.md" for instructions.


%% SET WORKING DIRECTORY AND OTHER PATHS

clear

% Project working directory
pathBase = '/Users/martakaras/Dropbox/_PROJECTS/OLD/LogisticPEER/gripeer-numerical-experiments';

% Directory with simulation scenario-specific connectivity matrices
scmatPath = strcat(pathBase,'/data/sc_matrices_used_scenario_1b_p198/'); 
% Add path to software 
addpath(strcat(pathBase,'/MATLAB/software')); 
% Define path to save results 
objsavepath = strcat(pathBase, '/results/results_scenario_1b.txt');
% Define path to save simulation status  
statusFilePath = strcat(pathBase, '/results/status_file_scenario_1b.txt');


%% LOAD PRE-COMPUTED CONNECTIVITY MATRICES

% File names of "base" connectivity matrices used to pre-compute 
% scenario-specific connectivity matrices
scMatrixNameGrid = {'A1.txt', 'A2.txt', 'A3.txt', 'A4.txt'};
% Values of "maximum" obtainable difference between a "base" connectivity
% matrix and matrix created via shuffling
scMatrixDiffMaxVals  = [1.0, 0.95, 0.67, 0.85]; 

% Array to store loaded connectivity matrices
scMatrixAobsList = {}; 
% Array to store values of "maximum" obtainable difference  
scMatrixDiffGridList = {};

% Iterate over number of "base" connectivity matrices
for iScMatrixIdx = 1:length(scMatrixNameGrid)   % iScMatrixIdx = 3; 
    % Get name of "base" connectivity matrix 
    iScMatrixName      = scMatrixNameGrid{iScMatrixIdx};
    % Get value of "maximum" obtainable difference between a "base" connectivity
    % matrix and matrix created via shuffling  
    iScMatrixDiffMax   = scMatrixDiffMaxVals(iScMatrixIdx);
    % Get grid of values of difference between a "base" connectivity
    % matrix and matrix created via shuffling
	iScMatrixDiffGrid  = unique([linspace(0,iScMatrixDiffMax/2,6), linspace(iScMatrixDiffMax/2,iScMatrixDiffMax,3)]); 
    % Array to store loaded connectivity matrices correcponding to "base"
    % connectivity matrix in this loop iteration
    scMatrixAobsList_TMP = {}; 
    
    % Iterate over number of values of difference between a "base" connectivity
    % matrix and matrix created via shuffling
    for iScMatrixDiffIdx = 1:length(iScMatrixDiffGrid)     % iScMatrixDiffIdx = 3
        iScMatrixDiff  = iScMatrixDiffGrid(iScMatrixDiffIdx);
        Aobs_file_name = strcat(strrep(iScMatrixName,'.txt',''), '_diffidx_', num2str(iScMatrixDiffIdx), '.txt');
        % Load precomputed connectivity matrix
        Aobs = importdata(strcat(scmatPath, Aobs_file_name));   
        scMatrixAobsList_TMP{iScMatrixDiffIdx} = Aobs; 
    end
    
    scMatrixAobsList{iScMatrixIdx}     = scMatrixAobsList_TMP;
    scMatrixDiffGridList{iScMatrixIdx} = iScMatrixDiffGrid;
end 

% Check correctness of array object characteristics
length(scMatrixAobsList) == 4
length(scMatrixDiffGridList) == 4
sum(cellfun('length', scMatrixAobsList)) == 32
sum(cellfun('length', scMatrixDiffGridList)) == 32


%% SIMULATION SETUP 

% Seed the random number generator
rng(1);

% Simulation parameters 
m         = 1;
p         = 198;
nIter     = 100;
lamQtr    = 1;  
nGrid     = [100]; 

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

% Fails counters 
fails_cnt_griPEER           = 0;
fails_cnt_griPEER_c_ridge   = 0;

% Data objects (fixed across simulations)
Zsigma     = importdata(strcat(pathBase,'/data/Z_cov.txt'));
% Zsigma correction for p=198
Zsigma     = blkdiag(Zsigma,Zsigma,Zsigma);
Zmu        = zeros(p,1); 


%% SIMULATION RUN 

% Create simulation status file which stores information about simulation
% progress
if exist(statusFilePath, 'file')==2
    delete(statusFilePath);
    fileID = fopen(statusFilePath,'w');
    fclose(fileID);
end

% Define values for local code tetsing
% iScMatrixIdx = 1; in = 100; iScMatrixDiffIdx = 1; ii =1; 

tic()
% Iterate over SC matrices 
for iScMatrixIdx = 1:length(scMatrixAobsList)   
    fprintf('\n\n\niScMatrixIdx: %d', iScMatrixIdx);
    
    % A,Q  
    iScMatrixName      = scMatrixNameGrid{iScMatrixIdx};
    iScMatrixAobsGrid  = scMatrixAobsList{iScMatrixIdx};
    iScMatrixDiffGrid  = scMatrixDiffGridList{iScMatrixIdx}; 
    A                  = iScMatrixAobsGrid{1};    
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
                
    % Iterate over n number of observations
    for in = nGrid % iterate over n number of observations  
        fprintf('\n\n\nin: %d', in);

        Z         = mvnrnd(Zmu,Zsigma,in);
        Z         = zscore(Z);     
        Z         = Z/sqrt(in-1);
        X         = ones(in,1);    
        
        % Iterate over SC matrix permutation parameter
        for iScMatrixDiffIdx = 1:length(iScMatrixAobsGrid)  

            fprintf('\n\niScMatrixDiffIdx: %d', iScMatrixDiffIdx);
            iScMatrixDiff  = iScMatrixDiffGrid(iScMatrixDiffIdx); 
            Aobs           = iScMatrixAobsGrid{iScMatrixDiffIdx}; 
            Dobs           = diag(sum(abs(Aobs),1));   
            Qobs           = (Dobs^(-0.5))*(Dobs - Aobs)*(Dobs^(-0.5));
            
            % Iterate over experiment runs  
            for ii = 1:nIter  
                fprintf('\niter #%d', ii);
                
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
                    tEndGripeer = toc(tStartGripeer); 
                catch err
                    fails_cnt_griPEER = fails_cnt_griPEER + 1; 
                    fprintf('\n\ngriPEER estimation failed. Returning NaN(p, 1).');
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
                    tEndGripeer_c_ridge   = toc(tStartGripeer_c_ridge); 
                catch err
                    fails_cnt_griPEER_c_ridge = fails_cnt_griPEER_c_ridge + 1;
                    warning('\n\ngriPEER_c_ridge estimation failed. Returning NaN(p, 1).');
                end
                
                % Store simulation results 
                scMatrixNameVec        = [scMatrixNameVec; cellstr(repmat(iScMatrixName,p,1))]; 
                scMatrixPermParVec     = [scMatrixPermParVec; repmat(iScMatrixDiff,p,1)]; 
                nVec                   = [nVec; repmat(in,p,1)];
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
            end
            
            % Append status to a file
            fileID = fopen(statusFilePath,'a');
            fprintf(fileID,'iScMatrixIdx: %f, n: %f, iScMatrixDiffIdx: %f\n',iScMatrixIdx,in,iScMatrixDiffIdx);
            fclose(fileID);
            
        end
        
    end
end
time = toc();
fprintf('time: %d\n', time); % time: 

% Wrapping up the results
resTable = table(nominal(scMatrixNameVec),scMatrixPermParVec,nVec,iVec,bIdxVec,bTrueVec,...
                 bEstGripeerVec,bEstGripeer_c_ridgeVec,...
                 lambQGripeerVec,lambRGripeerVec,lambQGripeer_c_ridgeVec,lambRGripeer_c_ridgeVec,...
                 timeGripeerVec,timeGripeer_c_ridgeVec,...
                 'VariableNames',{'SC_mat' 'SC_mat_perm_val' 'n' 'i' 'b_idx' 'b_true' 'b_est_griPEER' 'b_est_griPEER_c_ridge' 'lambQ_Gripeer' 'lambR_Gripeer' 'lambQ_Gripeer_c_ridge' 'lambR_Gripeer_c_ridge' 'exectime_Gripeer' 'exectime_Gripeer_c_ridge' });
         
% Save results to file
writetable(resTable, objsavepath); 
