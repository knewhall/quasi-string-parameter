function da = grada(system, x, i, j, d_wrt) %d_wrt is which x with are differentiating wrt
    param = system.param;

    x = reshape(x, 2, []);

    if (i == 1 && j == 2) || (i == 2 && j == 1)
    %distance between bead 1 and bead 2
        if d_wrt == 1 || d_wrt == 3
            dist = x(1,2) - x(1,1);
        elseif d_wrt == 2 || d_wrt == 4
            dist = x(2,2)-x(2,1);
        else 
            dist = 0;
        end 

        da = (40*param*exp(20*(norm(x(:,2)-x(:,1))))*dist)/(norm(x(:, 2)-x(:,1))*(1+exp(20*(norm(x(:, 2)-x(:,1)))))^2);

    elseif (i == 1 && j == 3) || (i == 3 && j == 1)
    %distance between bead 1 and bead 3
        if d_wrt == 1 || d_wrt == 5
            dist = x(1,3) - x(1,1);
        elseif d_wrt == 2 || d_wrt == 6
            dist = x(2,3) - x(2,1);
        else
            dist = 0; 
        end 

        da = (40*param*exp(20*(norm(x(:,3)-x(:,1))))*dist)/(norm(x(:, 3)-x(:,1))*(1+exp(20*(norm(x(:, 3)-x(:,1)))))^2);
    
    elseif(i == 2 && j == 3) || (i == 3 && j == 2)
    %distance between bead 2 and bead 3
        if d_wrt == 3 || d_wrt == 5
            dist = x(1,3) - x(1,2);
        elseif d_wrt == 4 || d_wrt == 6
            dist = x(2,3) - x(2,2);
        else
            dist = 0; 
        end 

        da = (40*param*exp(20*(norm(x(:,3)-x(:,2))))*dist)/(norm(x(:, 3)-x(:,2))*(1+exp(20*(norm(x(:, 3)-x(:,2)))))^2);

    end 
end 

    
