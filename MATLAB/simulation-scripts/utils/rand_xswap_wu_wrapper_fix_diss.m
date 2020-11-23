
function [rug,ii,diss,maxDiss] = rand_xswap_wu_wrapper_fix_diss(ug,targetDiss,maxIter,firstSeed,maxIterWrap)

if (firstSeed == 0) iSeed = 1; else iSeed = firstSeed; end
L       = full(sum(abs(ug(:))));
%L       = full(sum(ug(:)));
ii      = 0;
iRug    = ug;
iDiss   = 0;
maxDiss = 0; 

while ((iDiss < targetDiss) && (ii < maxIterWrap))
   ii           = ii + 1; 
   iSeed        = iSeed + 1; 
   [newRug,~,~] = rand_xswap_wu(iRug,1,maxIter,iSeed); 
   %iDiss        = full(sum(abs(ug(:)-newRug(:)))/(L*2)); 
   iDiss        = full(sum(abs(abs(ug(:))-abs(newRug(:))))/(L*2)); 
   iRug         = newRug;
   fprintf('ii = %f\niDiss = %f\n\n',ii,iDiss);  
   if (maxDiss < iDiss) 
       maxDiss = iDiss;
   end 
   if (iDiss >= targetDiss)
       break;
   end 
end 
if iDiss < targetDiss
    warning('target diss: %f, achieved diss: %f',targetDiss,iDiss)
end
rug  = full(iRug); 
diss = iDiss; 