
function diss = rand_xswap_wu_wrapper_diss(ug,numChanges,maxIter,seed)

if seed == 0
    seed = round(rand(1)*100000); 
end
[rug,eff,diss] = rand_xswap_wu(ug,numChanges,maxIter,seed); 


