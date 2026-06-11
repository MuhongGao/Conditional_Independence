function [Z,Y,X] = data_gener(n,q,p, model, H, rho)

% Function: generate data in simulation

switch model

    case 1

        Z = rand(n,q);

        if p <= q
            X_0 = Z(:, 1:p);
        else
            X_1 = Z(:, 1) * ones(1,p-q);
            X_0 = [Z, X_1];
        end

        Y_0 = Z(:, q);

        e1 = rand(n,1); 
        e2 = rand(n,p); 

        if H == 0
            Y = Y_0 + e1;
            X = X_0 + e2;

        elseif H == 1

            Y = Y_0 + rho * e2(:,1) + sqrt(1-rho^2) * e1 ;


            X = X_0 + e2;

        end


    case 2

        Z = randn(n,q);

        if p <= q
            X_0 = Z(:, 1:p);
        else
            X_1 = Z(:, 1) * ones(1,p-q);
            X_0 = [Z, X_1];
        end

        Y_0 = Z(:, q);

        e1 = randn(n,1); 
        e2 = randn(n,p); 

        if H == 0
            Y = Y_0 + e1;
            X = X_0 + e2;

        elseif H == 1

            Y = Y_0 + rho * e2(:,1) + sqrt(1-rho^2) * e1 ;


            X = X_0 + e2;

        end


end

end


