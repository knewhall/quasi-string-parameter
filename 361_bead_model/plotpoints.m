function plotpoints(pos, k)
    %k = size(pos,1);
    if nargin < 2
        k = 20000;% number of steps
    end
   % if nargin < 2
  %      c = [0 0 1];
   % else
   %     idx = ceil((20000-10000)/10);  
   %     idx = min(idx, length(pos));
    %    c = c(:, idx);
   % end
    
    % x,y,z,size,color,markertype
    scatter3(squeeze(pos(k,:,1)),squeeze(pos(k,:,2)),squeeze(pos(k,:,3)),15,[0 0 1],'filled', 'o');
    
    if 0
        x = squeeze(pos(end,:,:));
        plot3(x(:,1), x(:,2), x(:,3), 'b-');
    end
    
    daspect([1 1 1])
    pbaspect([1 1 1])
end