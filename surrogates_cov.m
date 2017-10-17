% function written by RIKKERT HINDRIKS
% X_rnd = surrogates(X,value)
% value == 0, then in-coherent (correlations destroyed)
% value == 1, then coherent (correlations preserved)
% Returns a data matrix X_rnd with the same sample covariance matrix as X
% (covariance structure at non-zero lags is not completely preserved). It
% hence generates data under the null hypothesis that the data are linear
% and stationary. Rows of X correspond to observations and columns to
% variables (i.e. TxN). X_rnd is constructed by adding the same random phase vector to
% the fft of the columns of X.

function X_rnd = surrogates_cov(X,value)

[N,M] = size(X);
% linearly interpolate if there are NaNs (same positions for all ROIs)
% go chunk by chunk (chunk=series of NaNs)
while any(isnan(X(:,1)))
    before = find(isnan(X(:,1)),1,'first')-1;
    nNaNs = find(~isnan(X(before+1:end,1)),1,'first');
    rplcmt = zeros(nNaNs,M);
    step = (X(before+nNaNs+1,:)-X(before,:))/(nNaNs+1);
    for m=1:M
        rplcmt(:,m) = step(m)*(1:nNaNs);
    end
    X(before+1:before+nNaNs,:) = rplcmt;
end

s = fft(X); %apply FFT

if value == 0
    phase_rnd=zeros(N,M);
    phase_rnd(1,1:M)=0;
    if (odd(N)==1);
        ph=2*pi.*rand(M,(N-1)/2)-pi;
        phase_rnd(2:N,1:M)=[ph';-flipdim(ph,2)'];
    end
    if (odd(N)==0);
        ph(1:M,1:(N-2)/2)=2*pi.*rand(M,(N-2)/2)-pi;
        phase_rnd(2:N,1:M)=[ph';zeros(1,M);-flipdim(ph,2)'];
    end
    
    %randomize the phases of all columns of s
    for m=1:M
        s_rnd(:,m) = abs(s(:,m)).*exp(i.*(angle(s(:,m)) + phase_rnd(:,m)));
        s_rnd(1,m) = s(1,m); %to make sure that s_rnd(1) is real
    end
end

if value == 1
    %construct (conjugate symmetric) array of random phases
    phase_rnd = zeros(N,1);
    %define first phase
    phase_rnd(1) = 0;
    if (odd(N) == 1);
        ph = 2*pi.*rand(1,(N-1)/2) - pi;
        phase_rnd(2:N) = [ph,-flipdim(ph,2)];
    end
    if (odd(N) == 0);
        ph(1:(N-2)/2) = 2*pi.*rand(1,(N-2)/2) - pi;
        phase_rnd(2:N) = [ph,0,-flipdim(ph,2)];
    end
    
    %randomize the phases of all columns of s
    s_rnd = zeros(N,M); %initialization
    for m=1:M
        %s_rnd(:,m) = abs(s(:,m)).*exp(i.*(angle(s(:,m)) + phase_rnd(:,m)));
        s_rnd(:,m) = abs(s(:,m)).*exp(i.*(angle(s(:,m)) + phase_rnd));
        s_rnd(1,m) = s(1,m); %to make sure that s_rnd(1) is real
    end
end

%apply inverse FFT
X_rnd = ifft(s_rnd,'symmetric'); %use "symmetric" because of round-off errors

%define output
output = X_rnd;

%---------------------------------------------------------
% o = odd(n) equals 1 if n is odd, 0 if n is even
    function outp = odd(n)
        
        for i = 1:round(n/2);
            if (2*i == n);
                outp = 0;
            else
                outp = 1;
            end
        end
        
        if (n == 0);
            outp = 0;
        end
    end
%---------------------------------------------------------

end

