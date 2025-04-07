function [dists] = beaddist(pos_t, dt)
    [d, p, nt] = size(pos_t);
     
    stride = 1;
    
   if nargin < 2 %plotbeaddist(pos, dt)
       kv = 1:stride:nt;
       vmax = nt;  
   else
       kv = linspace(0, dt*nt, nt);
       vmax = dt*nt;
   end
   dists = zeros(p*(p-1)/2, length(kv));
    
   for c=1:length(kv)
       dists(:, c) = pdist(pos_t(:, :, c)');
   end
   dists';
end         
