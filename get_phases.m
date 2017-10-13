function ph_data = extract_phases(data)
    ph_data = zeros(size(data));
    [T,N,S] = size(data);
    % filter data between 0.04 and 0.07 Hz
    filt_data = filter_fMRI(data,0.04,0.07);
    % obtain phase information
    for s=1:S
        for n=1:N
            Xanalytic = hilbert(filt_data(:,n,s));
            ph_data(:,n,s) = angle(Xanalytic);
        end
    end
    % first and last 10 volumnes removed because of finite sample effects in HT
    ph_data = ph_data(11:T-10,:,:); 
end