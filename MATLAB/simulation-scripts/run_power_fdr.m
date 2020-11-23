
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%---------------------------------------------------------%
%         Authors:    Marta Karas <mkaras2@jhu.edu>       %
%         Date:       Oct 13, 2019                        %
%---------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: The toolbox 'penalized' should be installed before running griPEER function.
% See "README.md" for instructions.


%% SET WORKING DIRECTORY AND OTHER PATHS

clear

% Project working directory
pathBase = '/Users/martakaras/Dropbox/_PROJECTS/LogisticPEER/gripeer-numerical-experiments';

% Directory with simulation scenario-specific connectivity matrices
scmatPath = strcat(pathBase,'/data/sc_matrices_used_scenario_1a/'); 
% Add path to software 
addpath(strcat(pathBase,'/MATLAB/software')); 
% Define path to save results 
objsavepath = strcat(pathBase, '/results/results_power_fdr.txt');
% Define path to save simulation status  
statusFilePath = strcat(pathBase, '/results/status_file_power_fdr.txt');


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
p         = 66;
nIter     = 100;
lamQtr    = 1;  
nGrid     = 1000; 
nGrid2    = 150;
      
% Objects to store simulation results
scMatrixNameVec        = []; 
scMatrixPermParVec     = []; 
nVec                   = [];
iVec                   = [];
POW_ASMP = [];
POW_BOOT = [];
FDR_ASMP = [];
FDR_BOOT = [];
RelCoeff_size = [];
IrrCoeff_size = [];

% Fails counters 
fails_cnt_griPEER           = 0;
fails_cnt_griPEER_c_ridge   = 0;

% Data objects (fixed across simulations)
Zsigma     = importdata(strcat(pathBase,'/data/Z_cov.txt'));
Zmu        = zeros(p,1); 


%% SIMULATION RUN 

% Create simulation status file which stores information about simulation
% progress
if exist(statusFilePath, 'file')==2
    delete(statusFilePath);
    fileID = fopen(statusFilePath,'w');
    fclose(fileID);
end

tic()
% Iterate over SC matrices 
for iScMatrixIdx = 2  
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
    for in = nGrid 
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

                % CI "true"
                XZ            =  [X, Z];
                Psi           =  diag(BernPr);
                estimVar      =  (XZ'*Psi*XZ)^(-1);
                CI_delta0     =  norminv(1-0.05/2) * sqrt(diag(estimVar));
                CI_delta      =  CI_delta0(2:end);
                
                % griPEER family estimation parameters
                GRIPEER_stopCrit    = 1e-6; 
                GRIPEER_initLambs   = [0.0001; 0.00001];
                GRIPEER_maxIter     = 20; 
                GRIPEER_lambsGrid   = [0.001, 0.001, 0.001, 0.0001, 0.0001, 0.0001, 0.00001, 0.00001;   0.001, 0.0001, 0.00001, 0.001, 0.0001, 0.00001, 0.001, 0.0001];

                % Results
                TrueBand      =  [bTrue - CI_delta, bTrue + CI_delta ];
                RelCoeff      =  find(TrueBand(:,1).*TrueBand(:,2)>0)';
                IrrCoeff      =  setdiff(1:p, RelCoeff);
                ygP           =  y(1:nGrid2);
                XgP           =  X(1:nGrid2,:);
                ZgP           =  Z(1:nGrid2,:);
                ZgP           =  zscore(ZgP);     
                ZgP           =  ZgP/sqrt(nGrid2-1);
                OutGripeer    =  griPEER(ygP, XgP, ZgP, Qobs, 'lambsGrid', GRIPEER_lambsGrid, 'maxIter', GRIPEER_maxIter, 'stopCrit', GRIPEER_stopCrit, 'initLambs', GRIPEER_initLambs, 'ciType', 'both');

                EstRelCoeff_asymp = find(OutGripeer.band(2:end,1).*OutGripeer.band(2:end,2)>0)';
                EstIrrCoeff_asymp = setdiff(1:p, EstRelCoeff_asymp);
                TrueDiscAsmp      = intersect(EstRelCoeff_asymp, RelCoeff);
                FalseDiscAsmp     = intersect(EstRelCoeff_asymp, IrrCoeff);

                EstRelCoeff_boots = find(OutGripeer.ciboot(2:end,1).*OutGripeer.ciboot(2:end,2)>0)';
                EstIrrCoeff_boots = setdiff(1:p, EstRelCoeff_boots);
                TrueDiscBoot      = intersect(EstRelCoeff_boots, RelCoeff);
                FalseDiscBoot     = intersect(EstRelCoeff_boots, IrrCoeff);

                % Store simulation results 
                scMatrixNameVec        = [scMatrixNameVec; cellstr(repmat(iScMatrixName,1,1))]; 
                scMatrixPermParVec     = [scMatrixPermParVec; repmat(iScMatrixDiff,1,1)]; 
                nVec                   = [nVec; repmat(in,1,1)];
                iVec                   = [iVec; repmat(ii,1,1)]; 
                POW_ASMP = [POW_ASMP, size(TrueDiscAsmp)/size(RelCoeff)];
                POW_BOOT = [POW_BOOT, size(TrueDiscBoot)/size(RelCoeff)];
                FDR_ASMP = [FDR_ASMP, size(FalseDiscAsmp)/size(EstRelCoeff_asymp)];
                FDR_BOOT = [FDR_BOOT, size(FalseDiscBoot)/size(EstRelCoeff_boots)];
                RelCoeff_size = [RelCoeff_size, size(RelCoeff,2)];
                IrrCoeff_size = [IrrCoeff_size, size(IrrCoeff,2)];
                    
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
resTable = table(nominal(scMatrixNameVec),scMatrixPermParVec,nVec,iVec,...
                 POW_ASMP.',POW_BOOT.',FDR_ASMP.',FDR_BOOT.',...
                 RelCoeff_size.',IrrCoeff_size.',...
                 'VariableNames',{'SC_mat' 'SC_mat_perm_val' 'n' 'i' 'POW_ASMP' 'POW_BOOT' 'FDR_ASMP' 'FDR_BOOT' 'RelCoeff_n' 'IrrCoeff_n'});
              
% Save results to file
writetable(resTable, objsavepath); 

