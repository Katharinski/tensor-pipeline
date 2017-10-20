% filter fMRI data
% input: data array with time x ROIs/seeds/voxels x subjects/sessions
% output: signal filtered between 0.04 and 0.07 Hz (Glerean et al. 2012,
% paper 11) and centered around zero

function filt_signal = filter_fMRI(data,flp,fhi)
    filt_signal = zeros(size(data));
    N = size(data,2);
    nsubj = size(data,3);
    % bandpass filter
    delt = 2;            % sampling interval
    k = 2;               % 2nd order butterworth filter
    for i=1:nsubj
        for seed = 1:N
            x = data(:,seed,i);
            fnq = 1/(2*delt);       % Nyquist frequency
            Wn = [flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
            [b,a] = butter(k,Wn);   % construct the filter
      
            fx = filtfilt(b,a,x);    % zero phase filter the data
            filt_signal(:,seed,i) = fx;
        end    
    end
    filt_signal = rm_mean_data(filt_signal,1);
end