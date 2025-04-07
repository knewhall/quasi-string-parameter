function d2W = d2W(system, x, delta)

    if nargin < 3
        delta = 1e-6;
    end

    format long 
    x = x';
    d2W = zeros(6, 1);
    dW = compute_dW(system, x, dx, guess);

    
    for i = 1:size(d2W,1)
        deltavec = zeros(1, 6);
        deltavec(i) = delta;
        xi = x + deltavec';
        [~, ri] = switching(system, xi);
        dr(:, i) = ri - r;
    end 
    dr = dr./delta; 
end

