function out = mixing(pos, dist)
    if nargin < 2
        dist = 50;
    end
    nt = size(pos,1);  %number of steps
    N = size(pos,2);    %nubmer of beads
    met = logical(zeros(N,N));
    for i=1:nt
        pts = squeeze(pos(i,:,:)); %361x3 position of beads at that step
        D = squareform(pdist(pts));
        met = met | (D<dist);
    end
    % count results
    met = met - eye(N);
    s = sum(met(:));
    s = s/(N*(N-1));
    
    out = s;
end