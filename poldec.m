function [U, H, its] = poldec(A, tol)
%POLDEC Polar decomposition.
% [U, H, ITS] = poldec(A) computes the polar decomposition A = U*H
% of the square, nonsingular matrix A. ITS is the number of
% iterations for convergence.

    % Initialize X (step 0) as A.
    X = A;
    % Start iteration counter.
    its = 0;
    inf_norm = norm(eye(size(X, 1)) - X'*X, Inf);
    while inf_norm>tol
        % Check whether X has converged to U by checking norm of I - X*X
        if (norm(X, 1)*norm(X, Inf))^0.5 < 3^0.5
            % If the convergence is guaranteed for Newton-Schulz, we 
            % prefer the Newton-Schulz iteration. Otherwise we proceed
            % with another step of Newton.
            % The real test for guaranteed convergence involves the 2-norm
            % of X, but it is expensive to calculate, so we decide to use 
            % the property that the square root of the product of 1-norm 
            % and inf-norm are greater than 2-norm for the problem. 
            X_1 = X*(3*eye(size(X, 1)) - X'*X)*0.5;
            fprintf('%i', its);
            fprintf('%s', ' N-S: ');
        else
            X_1 = (X + inv(X'))/2;
            fprintf('%i', its);
            fprintf('%s', ' N: ');
        end
        delta = norm(X_1-X, Inf)/norm(X, Inf);
        fprintf('%f', delta);
        fprintf('%s', ', ');
        X = X_1;
        inf_norm = norm(eye(size(X, 1)) - X'*X, Inf);
        fprintf('%f', inf_norm);
        fprintf('%s', '    ');
        its = its+1;
    end
    % When the convergence condition is met, X has converged to U.
    U = X;
    % Since A = U*H => U'*A = U'*U*H => U'*A = H.
    H = U'*A;
    H = 0.5 * (H + H');
end