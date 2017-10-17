function [rel_Error,sig_Error] = get_error_hamming(Yd,X)
    F_kten = ktensor(Yd); 
    F_kten = arrange(F_kten);
    % binarize result tensor such that nnz is equal to that in data
    res_tens = double(F_kten);
    res_tens = res_tens(:);
    N = numel(double(X));
    npts = nnz(double(X));
    vals = res_tens~=0;
    [n,x] = hist(res_tens(vals),1000);
    n = n/numel(res_tens);
    n1 = zeros(1,999);
    for i=1:999
        n1(i) = sum(n(end:-1:end-i));
    end
    ithr = find(n1>=(npts/N),1,'first');
    thr = x(1000-ithr);
    res_tens_bin = double(res_tens>=thr);
    % calculate Hamming distance between result and data tensor
    N_c1 = nnz(res_tens_bin.*double(X(:))); % # of coinciding ones
    npts_res = nnz(res_tens_bin);
    E_N_c1 = 1/N * npts_res * npts; % expected # of coinciding ones
    rel_Error = 1-N_c1/((nnz(res_tens_bin)+nnz(double(X)))/2);
    sig_Error = 1-(N_c1-E_N_c1)/((nnz(res_tens_bin)+nnz(double(X)))/2);
end