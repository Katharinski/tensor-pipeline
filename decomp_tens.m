% function to decompose tensors created with make_tensor

% see software requirements below!!

% inputs: tens   - tensor, Matlab double, NxNxWxS, where N-#ROIs, 
%                  W-#windows, S-#subjects
%         F      - number of features 
%         thr    - binarization threshold, percentile - i.e. all values
%                  smaller than this percentile are set to 0, the rest to 1
%         nonneg - logical; true if tensor is nonnegative, i.e. created
%                  using mutual information ('MI') or absolute value of 
%                  correlation ('abs'); false if tensor has also negative 
%                  values, i.e. created using correlation ('corr')
% outputs: spfeats - spatial features, i.e. maps, NxFxS, contains
%                    membership weights for each ROI
%          tfeats  - associated time courses, TxNxS, contains activation
%                    values for each map
%          lambdas - spfeats and tfeats are normalized to norm=1, lambda is
%                    the value with which the resulting rank-1 tensor has 
%                    to be multiplied in order to get its overall weight in 
%                    the superposition
%          err     - error between original and reconstructed tensor
%
% SOFTWARE REQUIREMENTS (URLS as of 02/2017)
% MATLAB Tensor Toolbox Version 2.6 
% URL: http://www.sandia.gov/~tgkolda/TensorToolbox/
% tensor_toolbox
% URL: https://sites.google.com/site/jingukim/home#nmfcode
% TENSORBOX
% URL: http://www.bsp.brain.riken.jp/~phan/tensorbox.php 
%
% REFERENCES
% Tensor Toolbox: Brett W. Bader, Tamara G. Kolda and others. MATLAB Tensor
% Toolbox Version 2.6, 2015
% ncp-algorithm: Fast Nonnegative Tensor Factorization with an 
%                Active-set-like Method. Jingu Kim and Haesun Park.
%                In: High-Performance Scientific Computing: Algorithms and 
%                Applications, Springer, 2012, pp. 311-326.
% CPD-algorithm: Anh-Huy Phan, Petr Tichavsky and Andrzej Cichocki, 
%                "TENSORBOX: a Matlab package for tensor decomposition", 
%                2013


function [spfeats,tfeats,lambdas,err] = decomp_tens(tens,F,thr,nonneg) 
    [N,~,~,S] = size(tens);
    % to determine percentile, use only upper triangular of each slice
    unique_pairs_IDs = logical(repmat(triu(ones(N),1),[1,1,T]));
    spfeats = zeros(N,K,S);
    tfeats = zeros(T,K,S);
    lambdas = zeros(K,S);
    err = 0;
    for s=1:S
        loc_tens = squeeze(tens(:,:,:,s));
        if nonneg % -->use nonnegative decomposition (ADD REF)
            % due to numerics in the MI algorithm, some values can be
            % negative, but that's in places with very low MI, set to 0
            loc_tens(loc_tens<0) = 0;
        end
        if thr~=0
            % determine threshold
            tens_flat = loc_tens(unique_pairs_IDs);
            thr_val = prctile(tens_flat,thr);
            % binarize
            loc_tens = double(tens(:,:,:,s)>=thr_val);
        end
        tenstens = tensor(loc_tens); % make tensor object for algorithms
        if nonneg
            if thr==0
                [Yd,~,~,~] = ncp(tenstens,F);
            else
                [Yd,~,~,~] = ncp_hamming(tenstens,k);
            end
        else
            cp_param = cp_fLMa;
            cp_param.init = {'dtld' 'nvec' 'random'};
            if thr==0
                [Yd,~] = cp_fLMa(tenstens,F,cp_param);
            else
                [Yd,~] = cp_fLMa_hamming(tenstens,F,cp_param);
            end
        end
    end
    spfeats(:,:,s) = Yd.U{1};
    tfeats(:,:,s) = Yd.U{3};
    lambdas(:,s) = Yd.lambdas;
end