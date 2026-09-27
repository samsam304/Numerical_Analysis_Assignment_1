%% Inputs
syms x      % Initialize x as a symbolic variable

f_sym = x^2 - 2*x*exp(-x) + exp(-2*x);

df_sym = diff(f_sym, x);        % First derivative
ddf_sym = diff(f_sym, x, 2);    % Second derivative

f = matlabFunction(f_sym);
df = matlabFunction(df_sym);
ddf = matlabFunction(ddf_sym);

tol = 1e-4;     % Tolerance for convergence
n = 10000;      % Maximum iterations

p0 = 1;

% Standard Newton's Method
[root, iter] = newtonMethod(f, df, p0, tol, n);

fprintf('\nNewton''s Method:\n');
fprintf('The root found is: %.10f\n', root);
fprintf('Number of iterations: %d\n', iter);

% Modified Newton's Method
[root_mod, iter_mod] = modifiedNewtonMethod(f, df, ddf, p0, tol, n);

fprintf('\nModified Newton''s Method:\n');
fprintf('The root found is: %.10f\n', root_mod);
fprintf('Number of iterations: %d\n', iter_mod);

%% Newton's Method
function [root, iter] = newtonMethod(f, df, p0, tol, n)

iter = 0;

% fprintf('Iter    p0          p           f(p0)       |p - p0|\n');
% fprintf('------------------------------------------------------------\n');

vals = ones(1,4);

while iter < n

    iter = iter + 1;

    fp0 = f(p0);
    dfp0 = df(p0);

    % Check that derivative is not zero
    if dfp0 == 0
        error('Derivative is zero at p0 = %.6f. Newton''s method fails.', p0);
    end

    % Compute next approximation
    p = p0 - fp0 / dfp0;

    % Order of convergence
    vals(iter) = p;

    % fprintf('%4d  %10.6f  %10.6f  %12.6e  %12.6e\n', ...
    %     iter, p0, p, fp0, abs(p - p0));

    % Check convergence
    if abs(p - p0) < tol
        root = p;
        fprintf('\nNewton''s method converged successfully.\n');
        if iter >= 3
            alpha = log(abs(vals(iter) - vals(iter-1))) / ...
                log(abs(vals(iter-1) - vals(iter-2)));

            fprintf('Order of Convergence: %.10f\n', alpha);
        end
        return
    end

    % Update approximation
    p0 = p;
end

% Maximum iterations reached
root = p;
fprintf('\nMaximum number of iterations reached.\n');

end

%% Modified Newton's Method
function [root, iter] = modifiedNewtonMethod(f, df, ddf, p0, tol, n)

iter = 0;

% fprintf('Iter    p0          p           f(p0)       |p - p0|\n');
% fprintf('------------------------------------------------------------\n');

vals = zeros(1, n);

while iter < n

    iter = iter + 1;

    fp0 = f(p0);
    dfp0 = df(p0);
    ddfp0 = ddf(p0);

    % Denominator from modified Newton's formula
    denom = dfp0^2 - fp0 * ddfp0;

    % Check that denominator is not zero
    if denom == 0
        error(['Denominator is zero at p0 = %.6f. ' ...
            'Modified Newton''s method fails.'], p0);
    end

    % Compute next approximation
    p = p0 - (fp0 * dfp0) / denom;

    % Store approximations
    vals(iter) = p;

    % fprintf('%4d  %10.6f  %10.6f  %12.6e  %12.6e\n', ...
    %     iter, p0, p, fp0, abs(p - p0));

    % Check convergence
    if abs(p - p0) < tol
        root = p;
        fprintf('\nModified Newton''s method converged successfully.\n');

        % Order of convergence
        if iter >= 3
            alpha = log(abs(vals(iter) - vals(iter-1))) / ...
                log(abs(vals(iter-1) - vals(iter-2)));

            fprintf('Order of Convergence: %.10f\n', alpha);
        end

        return
    end

    % Update approximation
    p0 = p;

end

% Maximum iterations reached
root = p;
fprintf('\nMaximum number of iterations reached.\n');

end