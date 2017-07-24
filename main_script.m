% tensor main script that shows all the steps in the workflow

% 1) create a tensor from your fMRI data with dimensions time x brain
% regions/channels x subjects; we use correlation because that's fastest
data = load('ex_data');
[T,N,S] = size(data);
w = 60;
method = 'corr';
tens = make_tensor(data,w,method);
W = size(tens,3);

% 2) decompose all tensors and store resulting features
% use different numbers of features and thresholds to compare later
F = 3:5;
Thrs = [75,90,99];
nonneg = 0; % method = correlation --> not non-negative
% save features, lambdas ("eigenvalues"), and errors
spfeats = cell(length(F),length(Thrs));
tfeats = cell(length(F),length(Thrs));
lambdas = cell(length(F),length(Thrs));
err = cell(length(F),length(Thrs));
% loop over number of features and thresholds
for f=F
    f_count = f_count+1;
    for thr=Thrs
        thr_count = thr_count+1;
        [spfeats{f_count,thr_count},tfeats{f_count,thr_count},lambdas{f_count,thr_count},err{f_count,thr_count}] = decomp_tens(tens,f,thr,nonneg);
    end
end

% 3) use K-means clustering to determine best f and thr