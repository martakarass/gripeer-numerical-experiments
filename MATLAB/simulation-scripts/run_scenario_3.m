
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
scmatPath = strcat(pathBase,'/data/sc_matrices_used_scenario_3/'); 
% Add path to software 
addpath(strcat(pathBase,'/MATLAB/software'));
% Define path to save results 
objsavepath = strcat(pathBase, '/results/results_scenario_3.txt');
% Define path to save simulation status  
statusFilePath = strcat(pathBase, '/results/status_file_scenario_3.txt');


%% SIMULATION SETUP

% Seed the random number generator
rng(1);

% Simulation parameters 
m         = 1;
p         = 66;
nIter     = 100;
lamQtr    = 1;  
nGrid     = [100]; 
      
scMatrixBaseNameGrid = {'A1.txt',...  
                        'A2.txt',...  
                        'A3.txt',...  
                        'A4.txt'};
iScMatrixDerivNameSuff = {'50','62','75','88','100','112','125','138','150'};

% Objects to store simulation results
scMatrixBaseNameVec    = []; 
scMatrixDerivValVec    = []; 
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
% iScMatrixBaseIdx = 1; in = 100; iScMatrixDerivIdx = 1; ii =1; 

tic()
% Iterate over SC matrices
for iScMatrixBaseIdx = 1:length(scMatrixBaseNameGrid)  
    fprintf('\n\n\niScMatrixBaseIdx: %d', iScMatrixBaseIdx);
    
    % A,Q  
    iScMatrixBaseName  = scMatrixBaseNameGrid{iScMatrixBaseIdx};
    iScMatrixDerivNamePref = strcat(iScMatrixBaseName(1:2), '_');
    A                  = importdata(strcat(scmatPath, iScMatrixBaseName));
    D                  = diag(sum(abs(A),1));
    Q                  = (D^(-0.5))*(D - A)*(D^(-0.5)); 
    
    % Qtilde, SIGMA
    [AA, BB, CC]       = svd(Q);
    Qrank              = rank(Q);
    BB2                = BB;
    BB2((Qrank+1):end,(Qrank+1):end) = 0.01 * eye(p-Qrank) * BB(Qrank, Qrank);
    BB2_diag           = diag(BB2);
    Qtilde_pre         = AA * diag(sqrt(1./BB2_diag)); 
    SIGMA              = Qtilde_pre * Qtilde_pre';
    
    % Iterate over n number of observations
    for in = nGrid 
        fprintf('\n\n\nin: %d', in);

        Z         = mvnrnd(Zmu,Zsigma,in);
        Z         = zscore(Z);     
        Z         = Z/sqrt(in-1);
        X         = ones(in,1);    
        
        % Iterate over SC matrix dissimilarity parameter
        for iScMatrixDerivIdx = 1:length(iScMatrixDerivNameSuff)
            fprintf('\n\niScMatrixDerivIdx: %d', iScMatrixDerivIdx);
            iScMatrixDerivName = strcat(iScMatrixDerivNamePref, iScMatrixDerivNameSuff{iScMatrixDerivIdx}, '.txt');
            Aobs               = importdata(strcat(scmatPath, iScMatrixDerivName));
            Dobs               = diag(sum(abs(Aobs),1) + repmat(0.01,p,1)');    
            Qobs               = (Dobs^(-0.5))*(Dobs - Aobs)*(Dobs^(-0.5));
            
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
                    warning('\n\ngriPEER_c_ridge estimation failed. Returning NaN(p, 1).');
                end
                
                % store simulation results 
                scMatrixBaseNameVec    = [scMatrixBaseNameVec; cellstr(repmat(iScMatrixBaseName,p,1))]; 
                scMatrixDerivValVec    = [scMatrixDerivValVec; cellstr(repmat(iScMatrixDerivNameSuff{iScMatrixDerivIdx},p,1))];
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
            
            % append status to a file
            fileID = fopen(statusFilePath,'a');
            fprintf(fileID,'n: %f, iScMatrixBaseName: %s, iScMatrixDerivName: %s\n',in,iScMatrixBaseName,iScMatrixDerivName);
            fclose(fileID);
            
        end
        
    end
end
time = toc();
fprintf('time: %d\n', time); % time: 

% Wrapping up the results
resTable = table(nominal(scMatrixBaseNameVec),nominal(scMatrixDerivValVec),nVec,iVec,bIdxVec,bTrueVec,...
                 bEstGripeerVec,bEstGripeer_c_ridgeVec,...
                 lambQGripeerVec,lambRGripeerVec,lambQGripeer_c_ridgeVec,lambRGripeer_c_ridgeVec,...
                 timeGripeerVec,timeGripeer_c_ridgeVec,...
                 'VariableNames',{'iScMatrixBaseName' 'scMatrixDerivValVec' 'n' 'i' 'b_idx' 'b_true' 'b_est_griPEER' 'b_est_griPEER_c_ridge' 'lambQ_Gripeer' 'lambR_Gripeer' 'lambQ_Gripeer_c_ridge' 'lambR_Gripeer_c_ridge' 'exectime_Gripeer' 'exectime_Gripeer_c_ridge' });

% Save results to file
writetable(resTable, objsavepath);  