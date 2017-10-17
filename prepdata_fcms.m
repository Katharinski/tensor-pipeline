% modified version from riken folder Jan 2016

% make tensor with windowed FCms for decomposition
% input: data of one session/subject, time x regions
%        w-width of window in frames
% ouput: tensor with pairwise FC values for each time step


function tens = prepdata_fcms(data,w)
[T,N] = size(data);
W = T-w+1;
tens = zeros(N,N,W);
for t=1:W
    chunk = data(t:t+w-1,:);
    tens(:,:,t) = corrcoef(chunk);
end
end