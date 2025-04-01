function [an_coeff, mix_coeff] = mixing_avgnearby_beta(param, trials)
    
    %run for a range of beta values 
    mix_coeff = zeros(1, length(param));
    an_coeff = zeros(1, length(param));

    for i = 1:length(param)
        beta = param(i);
        affinity = @(x) beta*2./(1+exp(20*(x-75)));

        mixingcoeff = zeros([1, trials]);
        avgnearby = zeros([1, trials]);

        for j = 1:trials 

            [~,pos,~] = modelsimulate(361, 297, affinity, 0.01, beta);
            mixingcoeff(j) = mixing(pos, 50);
            avgnearby(j) = averageNearby(pos, 50);

        end 

        mix_coeff(i) = sum(mixingcoeff)/trials;
        an_coeff(i) = sum(avgnearby)/trials;

    end 