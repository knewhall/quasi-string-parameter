function [p] = unbound(pos_t, dt, nt) 

    [dists] = beaddist(pos_t, dt);
    noncluster = zeros(1, nt);


    for i = 1:nt
        if dists(1,i) < 0.3 || dists(2, i) < 0.3 || dists(3, i) < 0.3
            noncluster(i) = 0; %beads in cluster
        else
            noncluster(i) = 1; %not in cluster 
        end 
    end

    p = sum(noncluster); %amount of time spent in non-cluster configuration

end