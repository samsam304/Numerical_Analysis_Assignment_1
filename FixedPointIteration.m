%{
Recall:

f(x) must be of the form x = f(x)
i.e. the eq must be manipulated to isolate a dependent variable

%}

%% Inputs
l = 3.1;
f = @(x) l*x*(1-x);    % Fixed-point function
p0 = 0.6;               % Initial approximation
tol = 1e-2;             % Tolerance for convergence
n = 500;                % Maximum iterations

[root, iter] = fixedPointIteration(f, p0, tol, n);

if root ~= Inf
    fprintf('\nThe approximate fixed point is: %.5f\n', root);
    fprintf('Number of iterations: %d\n', iter);

    cobwebPlot(f, p0, tol, n, [0.3 1]);
end

%% Fixed-Point Iteration Method
function [root, iter] = fixedPointIteration(f, p0, tol, n)

iter = 0;

fprintf('Iter    p0          p           |p - p0|\n');
fprintf('---------------------------------------------\n');

while iter < n

    iter = iter + 1;

    % Compute next approximation
    p = f(p0);

    fprintf('%4d  %10.7f  %10.7f  %12.7e\n', ...
        iter, p0, p, abs(p - p0));

    % Check convergence
    if abs(p - p0) < tol
        root = p;
        fprintf('\nFixed-point iteration converged successfully.\n');
        return
    
    elseif isinf(p) || isnan(p)
        root = Inf;
        fprintf('\nError: Iteration diverges.\n');
        return
    end

    % Update approximation
    p0 = p;
end

% Maximum iterations are reached
root = p;
fprintf('\nMethod failed after %d iterations.\n', n);

end

%% Cobweb Diagram
function cobwebPlot(f, p0, tol, n, xRange)

figure

% Plot f(x)
fplot(f, xRange, 'LineWidth', 1.5)
hold on

% Plot y = x
fplot(@(x) x, xRange, '--', 'LineWidth', 1.5)

x = p0;     % Starting point

for i = 1:n

    y = f(x);

    % Vertical line:    (x,x) --> (x,f(x))
    plot([x x], [x y], '-', 'Color', 'r')

    % Horizontal line:  (x,f(x)) --> (f(x),f(x))
    plot([x y], [y y], '-', 'Color', 'r')

    if abs(y - x) < tol
        break
    end

    % New iterate
    x = y;
end

grid on
xlabel('x')
ylabel('y')
title('Cobweb Diagram for Fixed-Point Iteration')

axis equal
xlim(xRange)
ylim(xRange)

hold off
end