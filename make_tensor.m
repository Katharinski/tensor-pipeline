% Function to create tensor from multichannel, multisubject data

% inputs: data   - TxNxS matrix (t-time points, N-number of ROIs/channels,
%                  S-number of subjects)
%         w      - window size (in frames)
%         method - 'corr' for correlation coefficient, 
%                  'abs' for absolute value of correlation coefficient, 
%                  'MI' for mutual information according to Kraskov et 
%                  al., 2004, 
%                  'ph' for adjacencies derived from phase differences,
%                  further methods can easily be implemented

function tens = make_tensor(data,w,method)
    [T,N,S] = size(data);
    W = T-w+1; % number of windows/layers in the tensor
    tens = zeros(N,N,W,S);
    for s=1:S
        switch method
            case {'corr','abs'}
                tens(:,:,:,s) = prepdata_fcms(data(:,:,s),w);
                
            case 'MI'
                tens(:,:,:,s) = MI_kraskov(data(:,:,s),w);
            case 'ph'
                tens = zeros(N,N,T,S);
                tens(:,:,:,s) = phase_diffs_adj(data(:,:,s));
        end
    end
    if strcmp(method,'abs')
        tens = abs(tens);
    end
end
