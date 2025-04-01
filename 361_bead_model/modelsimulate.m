function [t_arr, pos_t, bonds_t] = modelsimulate(N, T, affinity, meanon, param, varargin)
    % N: number of beads
    % T: temperature (such as: 297)
    % affinity: either the affinity function defining the CTMC switching
    %   matrix, or anything other than a function handle and it'll do the
    %   old switching with normal distributions

    % note: random hardcoded numbers are from the Hult model with reduced
    %   signicant figures

    % default parameters that can get overwritten by varargin
    pos = 60*randn(N,3);
    steps = 100000; %increase steps so that beads have time to interact 
                    %with every other bead
 
    dt = 0.001;
    confinement = 0.002;
    springScale = 1;
    maxdist = 90;
    
    use_waitbar = 1;

    % you can go back to making this configurable if you want to:
    %   0 - no dynamic crosslinkings
    %   1 - dynamic crosslinks with wormlike (nonlinear) springs
    %   2 - dynamic crosslinks with linear springs - if you use this 
    crosslinks = 2;

    if isa(affinity, 'function_handle')
        crosslink_style = 1;
    else
        crosslink_style = 0;
    end
    
    % handle optional arguments in varargin
    for i=1:2:length(varargin)
        if isequal(varargin{i},"pos")
            % starting bead positions
            pos = varargin{i+1};
        elseif isequal(varargin{i},"steps")
            % number of timesteps to take in simulation     idx = ceil((20000-10000)/10);  
   
            steps = varargin{i+1};
        elseif isequal(varargin{i},"dt")
            dt = varargin{i+1};
        elseif isequal(varargin{i},"confinement")
            % controls strength of confinement force keeping beads near the
            % origin
            confinement = varargin{i+1};
        elseif isequal(varargin{i},"springscale")
            % in the case that harmonic springs are used for the dynamic
            % crosslinking, this controls their strength
            springScale = varargin{i+1};
        elseif isequal(varargin{i},"maxdist")
            % Maximum distance at which binding probabilities are computed
            maxdist = varargin{i+1};
        end
    end
    
    % how many timesteps are actually getting recorded
    %recordsteps = N - recordstart + 1;

    %pos = ones(N,3).*linspace(-5000,5000,N)';
    
    % whether or not non-thermodynamic noise should scale with T
    % may be partially implemented
    scaleWithT = 1;

    
    break_mean = meanon/param; %divided by param value for switching noise 
    pos_t = zeros(steps+1,N,3);
    pos_t(1,:,:) = pos;
    std_t = zeros(steps+1,1);
    std_t(1)= mean(std(pos));
    
    bonds_t = zeros(N, steps+1);

    t_arr = linspace(0,steps*dt,steps+1);

    active = rand(N,1)>0.1;
    %meanon m=0.009;
    stdon = meanon/5;
    meanoff = meanon/9;
    stdoff = meanoff/5;

    if crosslink_style == 0
        switchtime = randn(N,1).*(stdon*active + stdoff*(1-active)) + rand(N,1).*(meanon*active + meanoff*(1-active));
    else
        switchtime = zeros(N,1);
    end
    bonds = zeros(N,1);
    
    dynScale = 10; %Hult model: 10.8794
    baseScale = dynScale/50;
    
    % decide which spring form to use: 1 -> wormlike, 2 -> harmonic
    if crosslinks == 1
        springScale = @(rn) dynScale * (-1 + 1/(1 - rn/45)^2 + (4*rn)/45);
    elseif crosslinks == 2
        springScale = @(rn) springScale*dynScale*rn;
    end
    
   
    if T>0
        scaleForNoise = 3.3234301797991788*sqrt(dt*T);
    end

    t = 0;
   
    if use_waitbar
        h = waitbar(0, 'Simulating');
    end
    for i=1:steps
        % compute excluded volume force
        pd = pdist(pos);
        D = squareform(pd);
        % 2.95294 dx E^(-((9 dx^2)/340000))
        %FEV_mag = 0.0332206*exp(-0.0000326797*(D.^2));
        FEV_mag = 0.03*exp(-0.000033333*(D.^2)); 
        % Hult model: FEV_mag = 0.0332206*exp(-0.0000326797*(D.^2))
        % FEV_mag = squareform(FEV_mag);
        fev = zeros(N,3);
        for j=1:N
            % Can this be matrixized?
            rel_pos = pos-pos(j,:);
            fev(j,:) = -FEV_mag(j,:)*rel_pos;
        end

        attr = zeros(N,3);
        if crosslinks
            % figure out attractions
            if crosslink_style == 0
                % bind any unbound beads
                bindable = (D > 0 & D<=90);
    %             create list of all possible binding pairs and distances
                dists = [];
                for j=1:N
                    if active(j)&&~bonds(j)
                        possible_targets = find(bindable(j,:) & ~bonds');
                        if isempty(possible_targets)
                            % nobody to bind to
                            continue
                        end

                        l = length(possible_targets);
                        d = D(j,possible_targets)';
                        dists = [dists; d, j*ones(l,1), possible_targets'];
                    end
                end
                x=1;
                dists = sortrows(dists);
                % pair up beads in order of ascending distance
                for k=1:size(dists,1)
                    a = dists(k,2);
                    b = dists(k,3);
                    if active(a) && active(b) && ~bonds(a) && ~bonds(b)
                        %fprintf("Forming a bond between %d and %d.\n", a, b);
                        bonds(a) = b;
                        bonds(b) = a;
                    end
                end
            else
                % New CTMC-based crosslinking

                % Update the crosslinking status
                
                % Compute the distance-based crosslinking expected time matrix
                
                % Figure out all possible pairs of bindable beads (how many
                % random variables do we need to draw?)
                bindable = (D > 0 & D<=maxdist);
                
                % convert affinity function to CTMC rates
                % we want to scale this so that affinity=1 leads to an
                % expected waiting time of meanoff
                mean_form_time = meanoff./affinity(D(bindable));
                
                % todo: filter out beads that are currently bound for
                % efficiency
                
                %error('check relationship between times and rates, also squareness')
                num_bindable = sum(sum(bindable));
                % create the exponential variables
                %bonds = zeros(N, N);
                form_times = exprnd(mean_form_time);
                [row, col] = ind2sub(size(bindable),find(bindable));
                [sorted, I] = sort(form_times);
                % mirror sorting onto row and col vectors
                row = row(I);
                col = col(I);
                
                for j=1:length(sorted)
                    % only go through times that fall in the current
                    % timestep
                    form_time = sorted(j);
                    if form_time > dt
                        break;
                    end
                    
                    b1 = row(j);
                    b2 = col(j);
                    if bonds(b1) == 0 && bonds(b2) == 0
                        % the if statement double checks that neither bead
                        % has already bound during this timestep
                        bonds(b1) = b2;
                        bonds(b2) = b1;
                        % draw a breaking time for the bond
                        
                        break_time = t + form_time + exprnd(break_mean);
                        switchtime(b1) = break_time;
                        switchtime(b2) = break_time;

                        % fprintf("Forming bond between %d and %d (%f - %f)\n", b1, b2, t, break_time);
                    end
                end
            end
            % store current binding information
            bonds_t(:, i) = bonds;
            
            
            % compute attractive forces

            for j=1:N
                if bonds(j) ~= 0
                    r = pos(bonds(j),:)-pos(j,:);
                    rn = norm(r);
                    if ~(rn > 0)
                        continue
                    end
                    if crosslinks == 1 && rn > 44
                        rn = 44;
                    end
                    n = r./norm(r); 

                    scale = springScale(rn);
                    
                    if dt*scale > rn*0.1
                        scale = rn*0.1/dt;
                    end
                    
                    v = scale*n;
                    
                    attr(j,:) = v;
                end
            end
        end
        
        
        if 0
            spring = zeros(N,3);
            dists = vecnorm(pos(2:end,:)-pos(1:(end-1),:), 2, 2);
            dists(dists>1650)=1650;
            baseforce = baseScale* (-1 + 1./(1 - dists/1700).^2 + (4*dists)/1700);
 
            spring(1:(end-1),:) = spring(1:(end-1),:) + baseforce.*(pos(2:end,:)-pos(1:(end-1),:))./dists;
            spring(2:(end),:) = spring(2:(end),:) - baseforce.*(pos(2:end,:)-pos(1:(end-1),:))./dists;
        else
            spring = 0;
        end
        
        % factor in external pressure (from EV of unsimulated beads)
        if confinement>0
            %k = .01;
            p = -confinement*pos.*vecnorm(pos,2,2);
        else
            p = zeros(N,3);
        end
        
        % add noise
        if T>0
            noise = scaleForNoise*normrnd(0,1,N,3);
        else
            noise = zeros(N,3);
        end

        % update positions
        dx = fev + attr + p + spring;
        
        % depending on setting, we should scale all force based on T
        if scaleWithT
            dx = (T/297)*dx;
        end
        pos = pos + dt*dx + noise;
        pos_t(i+1,:,:)=pos;
        cur_std = mean(std(pos));
        std_t(i+1) = cur_std;

        t = t + dt;
     

        % check if any beads have switched
        if crosslinks
            if crosslink_style == 0
                for j=1:N
                    if switchtime(j) < t
                        if active(j)
                            % deactivate bead
                            active(j) = 0;
                            % break current bond, if present
                            cur_bond = bonds(j);     idx = ceil((20000-10000)/10);  
   
                            if cur_bond > 0
                                bonds(cur_bond) = 0;
                                bonds(j) = 0;
                            end

                            % compute activation time
                            switchtime(j) = t + stdoff*randn + meanoff;
                        else
                            % activate bead
                            active(j) = 1;

                            % compute deactivation time
                            switchtime(j) = t + stdon*randn + meanon;
                        end
                    end
                end
            elseif crosslink_style == 1
                % Break bonds if expired
                for j=1:N
                    if switchtime(j) < t
                        bonds(j) = 0;
                    end
                end
            end
        end
        if use_waitbar
            waitbar(i/steps, h);
        end
    end
    if use_waitbar
        close(h)
    end
    
    %plot(t_arr,std_t)
    
    
end