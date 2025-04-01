function out = averageNearby(pos, dist)
    if nargin < 2
        dist = 100;
    end
    nt = size(pos,1);
    N = size(pos,2);
    counter = 0;
    for i=1:nt
        pts = squeeze(pos(i,:,:));
        D = pdist(pts);
        counter = counter + sum(D<dist);
    end
    % normalizing factor
    Z = nt*N;
    out = counter / Z;
end