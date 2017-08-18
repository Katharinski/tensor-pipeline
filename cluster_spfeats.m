% run K-means for different numbers of features and different values of the
% threshold; we cluster 5 times, and measure for distance is correlation 
% inputs:   spfeats - spatial features, as cell for different values of K
%                     and threshold; each cell contains an N x k x S  array
%                     with the spatial features returned by the tensor
%                     decomposition (decomp_tens)
%           maxntemplates - max. number of clusters you want to use; min.
%                           will be equal to k in each round
% outputs:  silh_vals - silhouette value for each k, number of clusters,
%                       and threshold (#K x maxntemplates-min(k) x #Thr)
%           clust_memb_IDs - cluster membership of each feature, in cell
%                            form, each cell contains a 1 x k*S-vector with
%                            integers from 1 to ntemplates
% 

function [silh_vals,clust_memb_IDs] = cluster_spfeats(spfeats,maxntemplates)
% 1st part - decompose for all thresholds and numbers of features
[nF,nThr] = size(spfeats);

for f_id=1:nF
    for thr_id=1:nThr
        % cluster and determine quality of clustering
        feats = spfeats{f_id,thr_id};
        if f_id==1 && thr_id==1
            [N,nF,S] = size(feats);
            nK  = length(nF:maxntemplates); % different numbers of clusters/templates
            silh_vals = zeros(nF,nK,nThr);
            clust_memb_IDs = cell(nF,nK,nThr);
        end
        F = size(feats,2);
        feats_flat = reshape(feats,[N,F*S]);
        k_id = f_id-1;
        for K = F:maxntemplates
            k_id = k_id+1;
            [idx,~,~] = kmeans(feats_flat',K,'Distance','correlation','Display','final','Replicates',5);
            clust_memb_IDs{f_id,k_id,thr_id} = idx;
            figure;
            [silh,~] = silhouette(feats_flat',idx,'correlation');
            close all
            %         h = gca;
            %         h.Children.EdgeColor = [.8 .8 1];
            %         xlabel 'Silhouette Value';
            %         ylabel 'Cluster';
            %         tit = sprintf('Mean = %.2f\n',mean(silh));
            %         title(tit);
            silh_vals(f_id,k_id,thr_id) = mean(silh);
        end
    end
end
end