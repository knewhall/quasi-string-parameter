function [out, flag] = compute_dW_1d(system, x, alpha, ~, guess, verbose, mult, shift)
    alpha = 1;
    S = alpha*system.S(x);
    V = system.V(x);
    
    maxiters = 20;
    
    p = guess;
    
    flag = 1;
    
    tol = 1e-10;
    
    for iter=1:maxiters
        Hv = H(p, S, V);
        Hpv = gradH(p, S, V);
        
        pnew = p - Hv/Hpv;
        err = abs(pnew-p);
        
        p = pnew;
        
        if err < tol
            flag = 0;
            break;
        end
    end
    
    out = p;
end