% main poldec function. takes a matrix A and a tolerance factor eps
function [U, H, its] = poldec(A, eps)
    % Initialize X (step 0) as A.
    X = A;
    % Start iteration counter.
    its = 0;
    Inf = infi_norm(eye(size(X, 1)) - X'*X);
    while Inf>eps
        % Check whether X has converged to U by checking norm of |I - X*X|
        if (one_norm(X)*infi_norm(X))^0.5 < 3^0.5
            % If the convergence is guaranteed for Newton-Schulz, we 
            % prefer the Newton-Schulz iteration. Otherwise we proceed
            % with another step of Newton.
            % The real test for guaranteed convergence involves the 2-norm
            % of X, but it is expensive to calculate, so we decide to use 
            % the property that the square root of the product of 1-norm 
            % and infni-norm are greater than 2-norm for the problem. 
            X_1 = n_schulz_step(X);
            fprintf('%i', its);
            fprintf('%s', ' N-S: ');
        else
            X_1 = newton_step(X);
            fprintf('%i', its);
            fprintf('%s', ' N: ');
        end
        delta = infi_norm(X_1-X)/infi_norm(X);
        fprintf('%f', delta);
        fprintf('%s', ', ');
        X = X_1;
        Inf = infi_norm(eye(size(X, 1)) - X'*X);
        fprintf('%f', Inf);
        fprintf('%s', '    ');
        its = its+1;
    end
    % When the convergence condition is met, X has converged to U.
    U = X;
    % Since A = U*H => U'*A = U'*U*H => U'*A = H.
    H = U'*A;
end

% Helper functions are defined here...
function Xk = newton_step(X)
    % performs one iteration of the newton iterations.
    Xk = (X + inv(X'))/2;
end

function Xk = n_schulz_step(X)
    % performs one iteration of the newton-schulz iterations.
    Xk = X*(3*eye(size(X, 1)) - X'*X)*0.5;
end

% I have defined all the norms i have used
function n = infi_norm(X)
    n = max(sum(abs(X')));
end

function n = one_norm(X)
    n = max(sum(abs(X)));
end