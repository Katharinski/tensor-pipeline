% function to remove mean from xdim data along dimension dim
function zero_mean_data = rm_mean_data(data,dim)

if(nargin==1)
   dim = 1;
   if(size(data,1) > 1)
      dim = 1;
   elseif(size(data,2) > 1)
      dim = 2;
   end
end

dims = size(data);
dimsize = size(data,dim);
dimrep = ones(1,length(dims));
dimrep(dim) = dimsize;

zero_mean_data = data - repmat(mean(data,dim),dimrep);
