% MI computation using the algorithm described in Kraskov, A., Stoegbauer,
% H., & Grassberger, P. (2004). Estimating mutual information. Physical
% Review E, 69(6), 066138., using version 1 described in eq. (8)

% input: data of single subject, TxN
%        w - window width in frames
% output: tensor, with MI values

function tens = MI_kraskov(data,w)
% empirical
[T,N] = size(data);
W = T-w+1;
tens = zeros(N,N,W);
fprintf('Computing MI...\n')
for t=1:W
    % display(t);
    chunk = data(t:t+w-1,:);
    if any(isnan(chunk(:)))
        nanrows = isnan(chunk(:,1));
        chunk = chunk(~nanrows,:);
    end
    % get the two TCs between which the MI is to be computed
    for n=1:N
        for m=n+1:N
            inputvec1 = chunk(:,n);
            inputvec2 = chunk(:,m);
            L = length(inputvec1);
            dists_x = nan(L);
            dists_y = nan(L);
            for p=1:L
                inds = p+1:L;
                % find nearest neighbor's distance (dist_z = max(dist_x,dist_y)
                % first compute all differences (symmetric)
                vals_x = abs(inputvec1(p)-inputvec1(inds));
                dists_x(p,inds) = vals_x;
                dists_x(inds,p) = vals_x;
                vals_y = abs(inputvec2(p)-inputvec2(inds));
                dists_y(p,inds) = vals_y;
                dists_y(inds,p) = vals_y;
            end
            % take bigger distance as the norm in z
            dists_z = max(cat(3,dists_x,dists_y),[],3);
            nrst_ngbr = min(dists_z,[],1);
            % count points closer than nearest neighbor
            n_x = sum(dists_x<repmat(nrst_ngbr,L,1));
            n_y = sum(dists_y<repmat(nrst_ngbr,L,1));
            psis = psi(n_x+1)+psi(n_y+1);
            I1 = psi(1)-mean(psis)+psi(L);
            tens(n,m,t) = I1;
            tens(m,n,t) = I1;
        end
    end
end
end