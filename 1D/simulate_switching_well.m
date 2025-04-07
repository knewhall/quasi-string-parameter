%function contains system which includes transition rate matrix and 
%drift matrix
%x: initial condition
%epsilon: additive Brownian noise
%nt: dtmax
%alpha: binding rate parameter
function [pos_t, state_t] = simulate_switching_well(system, x, dt, nt, eps, alpha, to_swap)
    
    assert(isscalar(x));
    
    if nargin < 5
        to_swap = 0;
    end
    
    c = 0.5; %unbinding rate
    
    t = 0;
    
    [S, r] = system.S(x); %transition rate matrix 
    break_mean = eps/(alpha*c);
    
    
    % The current binding configuration state (here, just 1 or 0)
    % state = 1 potential on 'bound'
    % state = 0 potential off 'unbound'
    state = 0;
    
    switchtime = 0;
    
    
    % initialize beads based on a draw from steady state
    state = rand < r(1);
    
    pos_t = zeros(1, nt+1);
    state_t = zeros(1, nt+1);
    pos_t(:, :, 1) = x;
    state_t(:, 1) = state;
    
    
    for step=1:nt
        
        
        % Update bead positions
        V = system.V(x);
        v = V(state+1);
        
        % update variables
        drift = dt*v;
        noise = sqrt(2*eps*dt)*normrnd(0,1);
        x = x + drift + noise;
        
        
        pos_t(step+1) = x;
        state_t(step+1) = state;
        
        t = t+dt;
        
        % check halting condition if applicable
        if to_swap
            % this assumes beads 1 and 2 start close to each other
            if abs(x) > 0.922633
            %if abs(x) > 1
                pos_t = pos_t(1:step+1);
                state_t = state_t(1:step+1);
                return
            end
        end
        
        if state && t >= switchtime
            % bond is ready to break
            state = 0;
        elseif ~state
            % see if bond forms
            Sv = system.S(x);

            sub_switchtime = exprnd(eps/(Sv(2,1)*alpha));
            if sub_switchtime < dt
                % bond forms
                state = 1;
                switchtime = t + exprnd(break_mean);
            end
        end
    end
end