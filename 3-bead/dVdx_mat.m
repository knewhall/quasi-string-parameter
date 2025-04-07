function [dVdx] = dVdx_mat(system, x, configs)
    d = system.d;
      phi = system.phi; %bead drag coefficient? 
    
    x = reshape(x, d, []);
    [d, p] = size(x);
    % binding spring force
    kv = system.kv;
    % confinement spring force
    kc = system.kc;
    % excluded volume strength
    a_ev = system.a_ev;
    % excluded volume distance
    c_ev = system.c_ev;
    phi = system.phi;
    if nargin < 3
        sc = memoize(@state_configurations);
        configs = sc(length(x));
    end
    n = size(configs,2);
    p = size(configs,1);
    dVdx = zeros([d p n]);

    for i=1:n
        state = configs(:, i);
        for k = 1:p
            if state(k) == k
                for j = 1:2
                    if j == 1
                        dVdx(j, k, i) = -kc*(3*x(j,k)^2+x(j+1,k)^2) ;
                    else 
                        dVdx(j, k, i) = -kc*(3*x(j,k)^2-x(j-1,k)^2) ;
                    end
                end
            else 
                for j = 1:2
                    if j == 1
                       dVdx(j, k, i) = kv - kc*(3*x(j,k)^2+x(j+1,k)^2) ;
                    else
                       dVdx(j, k, i) = kv - kc*(3*x(j,k)^2-x(j-1,k)^2) ;
                    end
                end
            end 
            for l=1:p
                if l==k
                    continue
                end
                for jj = 1:2
                    dVdx(jj, k, i) = dVdx(jj, k, i) + ...
                            a_ev * exp(-norm(x(:,k)-x(:, l)).^2 / c_ev) + ...
                            a_ev * (x(jj,k) - x(jj,l)) * exp(-norm(x(:,k)-x(:, l)).^2 / c_ev) * (-2*(x(jj,k)-x(jj,l))/ c_ev) ;
                end
            end
        end
    end 
    dVdx = reshape(dVdx, [d*p n])./phi;
end