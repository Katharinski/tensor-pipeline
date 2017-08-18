% function to make templates from saved features and silhouette values and
% save them

% input -   fname: file name where features and parameters used for decomposition
%           are stored; it is assumed that a correspondingly named
%           surrogate file exists
%           thr: optional; if only looking for the best results for a
%           certain threshold, give the percentile

function [templates,bestF,bestK,bestThr,corr_feats] = make_templates(fname,thr)
% load data
real = load(fname);
[N,~,S] = size(real.spfeats{1});
surr = load([fname,'_surr']);
if nargin==2
    real.silh_vals = real.silh_vals(:,:,real.thrs==thr);
    surr.silh_vals = surr.silh_vals(:,:,surr.thrs==thr);
end
% find best parameter combination
silh_vals_diff = real.silh_vals - surr.silh_vals;
[~,n] = max(silh_vals_diff(:));
[f_id,k_id,thr_id] = ind2sub(size(silh_vals_diff),n);
bestF = real.F(f_id);
ntemps = real.F(1):real.maxntemplates;
bestK = ntemps(k_id);
if nargin==2
    bestThr = thr;
else
    bestThr = real.Thrs(thr_id);
end
% load features and cluster memberships corresponding to parameter combi
feats = real.spfeats{f_id,thr_id};
%feats_flat = reshape(feats,[N,bestF*S]);
memb_IDs = real.clust_memb_IDs{f_id,k_id,thr_id};
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