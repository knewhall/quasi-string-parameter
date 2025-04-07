% Compute the deterministic pseudo-energy by integrating deterministic
% force along a path
function U = deterministic_U(system, x)
    d = system.d;
    if ndims(x) == 2
        x = reshape(x, d, [], size(x, 2));
    end
    [d, p, n1] = size(x);
    dx = (x(:, :,3:n1) - x(:, :,1:(n1-2)))/2;
    dU = zeros(d,p, n1);
    for j = 2:(n1-1)
        dU(:,:,j) = -deterministic_force(system, x(:,:,j));
    end
    % fill in dx at endpoints
    dx = cat(3, x(:,:,2)-x(:,:,1), dx, x(:,:,end)-x(:,:,(end-1)));
    U = squeeze(cumsum(sum(sum(0.5*(dU(:, :, 2:end).*dx(:, :, 2:end) + dU(:, :, 1:(end-1)).*dx(:, :,1:(end-1)))))))';
    U = [0 U];
    %U = squeeze(cumsum(sum(sum(dU.*dx, 1),2)));
end