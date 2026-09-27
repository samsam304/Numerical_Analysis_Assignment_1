% Inputs
f = @(x) -exp(-x) + sin(x);

p0 = 0;         % First initial approximation
p1 = 1;         % Second initial approximation
tol = 1e-6;     % Tolerance for convergence
n = 100;        % Maximum iterations

[root, iter] = secantMethod(f, p0, p1, tol, n);

% Output
fprintf('\nThe root found is: %.10f\n', root);
fprintf('Number of iterations: %d\n', iter);

%% Secant Method
function [root, iter] = secantMethod(f, p0, p1, tol, n)
    
iter = 1;

q0 = f(p0);
q1 = f(p1);

vals = ones(1,n);

fprintf('Iter    p0          p1          p           f(p)          |p - p1|\n');
fprintf('-----------------------------------------------------------------------\n');

while iter < n

    iter = iter + 1;

    % Check that denominator is not zero
    if q1 - q0 == 0
        error('Secant method fails because f(p1) - f(p0) = 0.');
    end

    % Compute next approximation
    p = p1 - q1 * (p1 - p0) / (q1 - q0);

    fp = f(p);
    
    % Store approximations
    vals(iter) = p;

    fprintf('%4d  %10.6f  %10.6f  %10.6f  %12.4e  %12.4e\n', ...
        iter, p0, p1, p, fp, abs(p - p1));

    % Check convergence
    if abs(p - p1) < tol
        root = p;
        fprintf('\nSecant method converged successfully.\n');
        if iter > 3
            alpha = log(abs(vals(iter) - vals(iter-1))) / log(abs(vals(iter-1) - vals(iter-2)));
            fprintf('Order of Convergence: %.10f', alpha);
        end
        return
    end

    % Update approximations
    p0 = p1;
    q0 = q1;
    p1 = p;
    q1 = fp;
end

% Maximum iterations reached
root = p;
fprintf('\nMaximum number of iterations reached.\n');

end