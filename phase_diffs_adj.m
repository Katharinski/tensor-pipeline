% function to get phase differences and turn them into adjacencies
% phase differences are only positive and are between 0 and pi
% adjacency are obtained with a linear function such that adj(pi)=0 
% (farthest) and adj(0)=1 (closest)

function tens = phase_diffs_adj(data)
    [T,N] = size(data);
    % first turn data into phases
    
    tens = zeros(N,N,T);
    pdiff = @(Xi,Xj)(abs(Xi-Xj));
    transf = @(x)(-1/pi*x+1);
    for t=1:T
        vec = pdist(data(t,:)',@(Xi,Xj) pdiff(Xi,Xj));
        vec(vec>pi) = 2*pi-vec(vec>pi);
        vec = transf(vec);
        tens(:,:,t) = squareform(vec);
    end
end
