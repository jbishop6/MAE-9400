function [P] = GEM(A, B, desc_X, desc_Y, num_eigenvalues, power)
    fprintf('Using GEM \n');
    n1 = size(A, 1);
    n2 = size(B, 1);

    %% Compute the Laplacian of the graphs
    Dh_A = diag(ones(n1,1)./sqrt(A*ones(n1,1)));
    Dh_A(Dh_A==Inf) = 0;
    Dh_A(Dh_A==-Inf) = 0;
    A = eye(n1) - Dh_A*A*Dh_A;
 
    Dh_B = diag(ones(n2,1)./sqrt(B*ones(n2,1)));
    Dh_B(Dh_B==Inf) = 0;
    Dh_B(Dh_B==-Inf) = 0;
    B = eye(n2) - Dh_B*B*Dh_B;  

    %% Initializing feature matrix if necessary
    if nargin < 4
         F1 = ones(n1,n2)';
    else  
        F1 = desc_X*desc_Y';  % Use shots as feature descriptors if needed
    end

    %% Eigen decomposition
    [U, Lambda] = eigs(A, num_eigenvalues, 1e-6);  % Use num_eigenvalues instead of hardcoding 40
    [V, Mu] = eigs(B, num_eigenvalues, 1e-6);  % Same for B

    lambda = diag(Lambda);
    mu = diag(Mu);

    %% Compute filtering coefficients
    if power
        Coeff = 1 ./ (((lambda.^10) - (mu.^10)').^2 + 1);
    else
        Coeff = 1 ./ (((lambda) - (mu)').^2 + 1);
    end

    %% Apply the filtering coefficients
    Coeff = Coeff .* (U' * F1 * V);
    Coeff = U * Coeff * V';

    %% Solve the matching problem
    M = matchpairs(Coeff', -99999, 'max');
    P = full(sparse(M(:, 1), M(:, 2), 1, n1, n1));

    return;
end
