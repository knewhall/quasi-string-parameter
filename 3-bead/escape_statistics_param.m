function [mean_times, all_times, stopflags] = escape_statistics_param(system, x, eps_range, trials)
x = find_minima(system, x); 
mean_times = zeros(size(eps_range));
all_times = cell(size(eps_range));
stopflags = cell(size(eps_range));
x2r = reshape(x, system.d, []);

dtmax = 100000000;
%dtmax = 1000000;

if nargin < 4
    trials = 5000;
end

parfor (e_c = 1:length(eps_range), 12)
    eps = eps_range(e_c);
    %disp(eps);
    times = zeros(1, trials);
    flags = zeros(1, trials);
    %h = waitbar(0, 'Simulating');
    for trial=1:trials
        %disp(trial)
        [~, ~, time,flag, ~] = simulate_system_param(system, x2r, 0.001, dtmax, eps, 1); 
        times(trial) = time;
        %stopflags = flag;
        flags(trial)= flag;
        %waitbar(trial/trials, h);
        if mod(trial, 20) == 0
            disp(trial);
        end
    end
    %close(h);
    avg_time = sum(times)/sum(times~=(dtmax));
    %disp(avg_time)
    stopflags{e_c} = flags;
    mean_times(e_c) = avg_time;
    all_times{e_c} = times;
    %hist(times) 
end
end