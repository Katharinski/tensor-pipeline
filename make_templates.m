% function to make templates from saved features and silhouette values and
% save them

% input -   fname: file name where features and parameters used for decomposition
%           are stored; it is assumed that a correspondingly named
%           surrogate file exists
%           thr: optional; if only looking for the best results for a
%           certain threshold, give the percentile

function [templates,bestF,bestK,bestThr,corr_feats] = make_templates(spfeats,clust_memb_IDs,silh_vals,silh_vals_surr,F,maxntemplates,Thrs,thr)
% load data
N = size(spfeats{1},1);
if nargin==8
    silh_vals = silh_vals(:,:,Thrs==thr);
    silh_vals_surr = silh_vals_surr(:,:,Thrs==thr);
end
% find best parameter combination
silh_vals_diff = silh_vals - silh_vals_surr;
[~,n] = max(silh_vals_diff(:));
[f_id,k_id,thr_id] = ind2sub(size(silh_vals_diff),n);
bestF = F(f_id);
ntemps = F(1):maxntemplates;
bestK = ntemps(k_id);
if nargin==2
    bestThr = thr;
else
    bestThr = Thrs(thr_id);
end
% load features and cluster memberships corresponding to parameter combi
feats = spfeats{f_id,thr_id};
%feats_flat = reshape(feats,[N,bestF*S]);
memb_IDs = clust_memb_IDs{f_id,k_id,thr_id};
% create templates by taking the mean inside each cluster
templates = zeros(N,bestK);
for t=1:bestK
    templates(:,t) = mean(feats(:,memb_IDs==t),2);
end
% compute correlations between ordered features for visualization purposes
[~,w] = sort(memb_IDs,'ascend');
sort_feats = feats(:,w);
corr_feats = corrcoef(sort_feats);
end