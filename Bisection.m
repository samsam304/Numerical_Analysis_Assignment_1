% Input 
f = @(x) -exp(-x) + sin(x);

interval = [0.5; 1];

a = interval(1);          % Lower bound
b = interval(2);          % Upper bound
tol = 1e-10;    % Tolerance for convergence
n = 100;        % Maximum iterations

[root, iter] = bisectionMethod(f, a, b, tol, n);

% Output
fprintf('\nThe root found is: %.10f\n', root);
fprintf('Number of iterations: %d\n', iter);

%% Bisection Method
function [root, iter] = bisectionMethod(f, a, b, tol, n)

    fa = f(a);
    fb = f(b);
    iter = 0;

    p = ones(1,30);

    % Check whether either endpoint is already a root
    if fa == 0
        root = a;
        return
    elseif fb == 0
        root = b;
        return
    end

    % Verify Intermediate Value Theorem holds
    if fa * fb > 0
        error('f(a) and f(b) must have opposite signs.');
    end

    % fprintf('Iter    a           b           mid         f(mid)\n');
    % fprintf('--------------------------------------------------------\n');

    while (b - a)/2 > tol && iter < n

        iter = iter + 1;

        midpoint = (a + b)/2;
        fm = f(midpoint);

        % fprintf('%4d  %10.6f  %10.6f  %10.6f  %12.6e\n', ...
        %         iter, a, b, midpoint, fm);

        % Exact root found
        if fm == 0
            root = midpoint;
            return
        end

        % Decide which half contains the root
        if fa * fm > 0
            a = midpoint;
            fa = fm;
        else
            b = midpoint;
        end

        p(iter) = midpoint;
        fprintf('p_n = %.10f\n', midpoint)
    end

    root = (a + b)/2;

    % Order of convergence
    alpha = log(abs(p(31) - p(30))) / log(abs(p(30) - p(29)));
    fprintf('Order of Convergence: %.10f', alpha);

    % Check why the loop stopped
    if (b - a)/2 <= tol
        fprintf('\nBisection converged successfully.\n');
    else
        fprintf('\nMaximum number of iterations reached.\n');
    end

end