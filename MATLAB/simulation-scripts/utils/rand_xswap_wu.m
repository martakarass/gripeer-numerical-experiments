% function [rug,eff,diss] = rand_xswap_wu(ug,numChanges,maxIter,seed)
%
% randomization technique that preserves:
% - size, density and strength
% - degree-sequence (and hence degree distribution)
%
% look for {i,k}{j,l} edges and change them for {i,l}{j,k} edges
% input parameters
% - ug: undirected adjacency matrix. Can be symmetric or upper-triangular,
% binary or weighted.
% - numChanges: number of effective "crossings" desired to be performed
% - maxIter: maximum number of attempts to satisfy numChanges
% - seed: specify seed for random-number-generator
% output paremeters
% - rug: randomized adjacency matrix in sparse upper-triangular format
% - eff: relative frequency of success randomizing
% - diss: dissimilarity between the original and the randomized adj matrices

% Joaquin Goni

function [rug,eff,diss] = rand_xswap_wu(ug,numChanges,maxIter,seed)

tolChange = 10^-4; %tolerance for density change in randomization

if ~issparse(ug)
    if islogical(ug)
        ug = sparse(ug);
    else
        ug = sparse(double(ug));
    end
end

ug = triu(ug);
rug = ug;
L = full(sum(rug(:))); %number of edges of the graph
N = size(ug,1); %number of nodes of the graph

if nargin<3
    maxIter=numChanges*100;
end

if nargin==4
    stream0 = RandStream('mt19937ar','Seed',seed); % example of creating a seed: seed = sum(clock.*([1 1 1 1 1 1000]));
    RandStream.setGlobalStream(stream0);
end

iter = 0;
changes = 0;

[nodesA,nodesB] = find(rug~=0); %original pairs of connected nodes (edges)    

while (changes < numChanges) && (iter<=maxIter)
    i = ceil(rand*N); %select i at random
    j = i; 
    while i==j
        j = ceil(rand*N); %select j at random
    end
    %%% random selection of node k within i-neighborhood
    nn = [nodesA(nodesB==i);nodesB(nodesA==i)];
    nn0 = ((nn==i)+(nn==j))==0; %neighborhood of i without {i,j}
    iNeigh = nn(nn0); %candidates to be k
    
    if ~isempty(iNeigh) %if there are candidates to become k 
        k = iNeigh(ceil(rand*numel(iNeigh)));
    
        %%% random selection of node l within j-neighborhood
        nn = [nodesA(nodesB==j);nodesB(nodesA==j)];
        nn0 = ((nn==i)+(nn==j)+(nn==k))==0; %neighborhood of j without {i,j,k}
        jNeigh = nn(nn0); %candidates to be l

        if ~isempty(jNeigh) %if there are candidates to become l 
            l = jNeigh(ceil(rand*numel(jNeigh)));

            if ~rug(min(i,l),max(i,l)) && ~rug(min(k,j),max(k,j))   %if edges {i,l} and {k,j} do not exist
                option = ceil(rand*2);
                switch option
                    case 1
                        rug(min(i,l),max(i,l)) = rug(min(i,k),max(i,k)); %add edge {i,l} with value of {i,k}
                        rug(min(j,k),max(j,k)) = rug(min(j,l),max(j,l)); %add edge {j,k} with value of {j,l}
                    case 2
                        rug(min(i,l),max(i,l)) = rug(min(j,l),max(j,l)); %add edge {i,l} with value of {j,l}
                        rug(min(j,k),max(j,k)) = rug(min(i,k),max(i,k)); %add edge {j,k} with value of {i,k}
                end
                    
                rug(min(i,k),max(i,k)) = 0; %remove edge {i,k} 
                rug(min(j,l),max(j,l)) = 0; %remove edge {j,l} 

                %% {i,j} is substituted by {i,l} in lists nodesA and nodesB
                ijpos = (nodesA==min(i,k)) & (nodesB==max(i,k));
                nodesA(ijpos) = min(i,l);
                nodesB(ijpos) = max(i,l);
                %% {j,l} is substituted by {j,k} in lists nodesA and nodes B
                klpos = (nodesA==min(j,l)) & (nodesB==max(j,l));
                nodesA(klpos) = min(j,k);
                nodesB(klpos) = max(j,k);
                if (nnz(ijpos)~=1) || (nnz(klpos)~=1)
                   warning('lists of connected pairs corrupted during randomization!')
                end
                %[nodesA,nodesB] = find(rug~=0); %update pairs of connected nodes (edges)    
                changes = changes + 1;
            end   
        end
    end
    iter = iter + 1; 
end

Lrand = full(sum(rug(:)));
if (abs(L-Lrand) > tolChange) %if change in density is great than tolerance
    warning('total strength of the graph altered during randomization!')
    warning('original %0.6f, randomized %0.6f',L,Lrand)
end

if nnz(ug)~=nnz(rug)
    warning('density of the graph altered during randomization!')
end

diss = full(sum(abs(ug(:)-rug(:)))/(L*2)); %relative dissimilarity (value within range [0,1])
eff = changes / iter; %relative frequency of effective changes with respect to attempts
rug = rug + rug';
