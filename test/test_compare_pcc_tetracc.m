%TEST_COMPARE_PCC_TETRACC
%
%   Requires the following modules:
%
%     cpuinfo-m
%     corr-m
%
%   Author: Kristian Loewe

nVoxels = 500;
nTimepts = 200;
res = zeros(nVoxels*(nVoxels-1)/2, 3);

for i = 1:numel(nVoxels)
    nV = nVoxels(i);
    
    for j = 1:numel(nTimepts)
        nT = nTimepts(j);
        
        dataFlt = rand(nT,nV,'single');
        
        tmp = corrcoef(dataFlt);
        res(:,1) = flat(tmp(logical(tril(ones(nVoxels),-1))));

        res(:,2) = pcc(dataFlt);
        
        res(:,3) = tetracc(dataFlt);
        
        figure; plotmatrix(res);
    end
end
