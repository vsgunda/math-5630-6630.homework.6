% Author: Vishnu Gunda / vsg0005@auburn.edu
% Date: 2024-12-02
% Assignment Name: hw06


classdef hw06
    methods (Static)
        
        % Problem 1

        function ret = p1(func, a, b, n, option)
    % Implement composite quadrature rules for numerical integration
    % of a function over the interval [a, b] with n subintervals.
    % The option parameter determines the type of quadrature rule to use: 
    % 1 for the midpoint rule, 2 for the trapezoidal rule, and 3 for Simpson's rule.

    % Check for even number of intervals for Simpson's rule
    if option == 3 && mod(n, 2) ~= 0
        error('n must be even for Simpsons rule');
    end

    % Step size
    h = (b - a) / n;

    % Initialize result
    ret = 0;

    if option == 1
        % Composite midpoint rule
        for i = 1:n
            xk = a + (i - 1) * h; % Left endpoint of the subinterval
            xk1 = xk + h;         % Right endpoint of the subinterval
            midpoint = (xk + xk1) / 2;
            ret = ret + h * func(midpoint);
        end
    elseif option == 2
        % Composite trapezoidal rule
        for i = 1:n
            xk = a + (i - 1) * h; % Left endpoint of the subinterval
            xk1 = xk + h;         % Right endpoint of the subinterval
            ret = ret + (h / 2) * (func(xk) + func(xk1));
        end
    elseif option == 3
        % Composite Simpson's rule
        for i = 1:2:n-1
            xk = a + (i - 1) * h;   % First endpoint of the interval
            xk1 = xk + h;           % Midpoint
            xk2 = xk + 2 * h;       % Second endpoint of the interval
            ret = ret + (h / 3) * (func(xk) + 4 * func(xk1) + func(xk2));
        end
    else
        error('Invalid option: %d', option);
    end
        end

        % Problem 2

        function p2()
            % run with the following command: hw06.p2(). Do not edit this function.
            %
            % It checks the convergence of the composite quadrature rules implemented in p1.
            %
            % Here we use some examples, 
            % f_1(x) = exp(x) 
            % f_2(x) = (1 - x^2)^3, this function's first 2 derivatives at the endpoints are zero.
            % f_3(x) = (1 - x^2)^5, this function's first 4 derivatives at the endpoints are zero.
            % f_4(x) = (1 - x^2)^7, this function's first 6 derivatives at the endpoints are zero.

            % Run this function will plot the figures for the convergence of the composite quadrature rules.
            % Make comments about the results obtained from the plots. 
            %
            % > For instance, you can comment on the convergence rates of the quadrature rules, and how they compare to the theoretical rates.
            % > Here are a few example questions you can answer in your comments:
            % > Does the Runge phenomenon of f1 (Runge's function) lower the convergence rate?
            % > Does Simpson's rule have a faster convergence rate than the other two for these examples?
            % > Based on your observations, what kind of functions can have a faster convergence rate with the given composite quadrature rules?

            % Write your comments here.
            %
            % For f1 all the rules converge, but Simpson's rule converges
            % the fastest. The Midpoint and Trapezoidal rules converge at
            % similar rates. Non-zero endpoints slightly hamper the
            % convergence
            %
            % For f2, the results are similar but everything is converging
            % faster when compared to f1. It is smoother at the boundaries
            % because the derivatives vanish there
            %
            % For f3, Simpson's rule converges even faster due to four
            % derivatives vanishing at x=1 and x=-1. It is closer to
            % 6th-order convergence
            %
            % For f4, six derivatives vanish at x=1 and x=-1, making
            % Simpson's rule converge even faster. Nothing else is
            % particularly special
            %
            % Does the Runge phenomenon of f1 (Runge's function) lower the convergence rate?
            %
            % Yes for f1, especially for the midpoint and trapezoidal rules. This
            % is shown from the slower reduction in error compared to f2,
            % f3, and f4. This is due to the non-zero derivatives at the
            % endpoints. Simpson's rule performs better, but is still
            % slower compared to f2, f3, and f4.
            %
            % Does Simpson's rule have a faster convergence rate than the other two for these examples?
            %
            % Yes, it does with every single sample function.
            %
            % Based on your observations, what kind of functions can have a faster convergence rate with the given composite quadrature rules?
            %
            % Smooth functions and ones that have more derivatives
            % vanishing at endpoints are able to achieve a faster
            % convergence rates. For example, the four functions slowly
            % converge faster and faster. Each continuing function has 2,
            % 4, and 6 vanishing derivatives and as the number of vanishing
            % derivatives increases, the rate of convergence increases.
            %

            f = {  @(x)exp(x),  @(x) (1 - x.^2 ).^3, @(x)(1 - x.^2).^5,  @(x) (1 - x.^2).^7} ;  % Define the integrand
            exact = [exp(1) - exp(-1), 32/35, 512/693 , 4096/6435];  % Define the exact integral
            n = 2.^(1:8);  % Define the number of subintervals
            for k = 1 : length(f)

                error = zeros(3, length(n));  % Initialize the error matrix with zeros

                % Calculate the approximate integral and the error for each quadrature rule and number of subintervals
                for i = 1 : length(n)
                    error(1, i) = abs(hw06.p1(f{k},-1, 1, n(i), 1) - exact(k));
                    error(2, i) = abs(hw06.p1(f{k},-1, 1, n(i), 2) - exact(k));
                    error(3, i) = abs(hw06.p1(f{k},-1, 1, n(i), 3) - exact(k));
                end

                % Plot the error against the number of subintervals using a log-log scale
                figure(k);
    
                loglog(n, error(1, :), 'r-+', 'LineWidth', 2);
                hold on;
                loglog(n, error(2, :), 'g-d', 'LineWidth', 2);
                loglog(n, error(3, :), 'b-x', 'LineWidth', 2);

                loglog(n, 1./ n.^2, 'm--', 'LineWidth', 1);
                loglog(n, 1./ n.^4, 'k-.', 'LineWidth', 1);
                loglog(n, 1./ n.^6, 'm--d', 'LineWidth', 1);
                loglog(n, 1./ n.^8, 'k--o', 'LineWidth', 1);

                xlabel('Number of subintervals');
                ylabel('Absolute error');
                title(sprintf('Convergence of composite quadrature rules for %s', functions(f{k}).function));
                legend('Midpoint rule', 'Trapezoidal rule', 'Simpson''s rule', '2nd order convergence', '4th order convergence', '6th order convergence', '8th order convergence', 'Location', 'best');
                grid on;
                hold off;
            end

        end

        
        % Problem 3

        function ret = p3(func, a, b, N, option)
            % Use Richardson extrapolation (hw05.p1) to implement the Romberg integration method.
        
            % Step sizes and quadrature results
            data = zeros(N+1, 1);       % Quadrature results
            powers = [];                % Extrapolation powers
        
            % Populate the quadrature results
            for i = 0:N
                num_intervals = 2^i;          % Number of subintervals
                if option == 3 && mod(num_intervals, 2) ~= 0
                    num_intervals = num_intervals + 1; % Ensure even number for Simpson's rule
                end
                data(i+1) = hw06.p1(func, a, b, num_intervals, option); % Compute composite quadrature
            end
        
            % Determine powers for Richardson extrapolation
            if option == 3
                powers = 4:2:(4 + 2*(N-1));     % Simpson's rule: powers are [4, 6, 8, ...]
            else
                powers = 2:2:(2 + 2*(N-1));     % Midpoint/Trapezoidal: powers are [2, 4, 6, ...]
            end
        
            % Use hw05.p1 for Richardson extrapolation
            ret = hw05.p1(data, powers);
        end


        % Problem 4
        
        function ret = p4()
            % Construct the Gauss quadrature rule using the roots of the Legendre polynomial of degree 6.
        
            % Step 1: Define intervals where roots are known to exist
            intervals = [ -1, -3/4; -3/4, -1/4; -1/4, 0; 0, 1/4; 1/4, 3/4; 3/4, 1 ];
        
            % Initialize arrays for roots and weights
            roots = zeros(6, 1);
            weights = zeros(6, 1);
        
            % Step 2: Find the roots using the bisection method
            for i = 1:6
                roots(i) = find_root(@hw06.legendre_poly_6, intervals(i, 1), intervals(i, 2), 1e-14);
            end
        
            % Step 3: Compute weights for the roots
            for i = 1:6
                weights(i) = 2 / ((1 - roots(i)^2) * (deriv_legendre_poly(6, roots(i))^2));
            end
        
            % Combine roots and weights into the result
            ret = [roots, weights];
        
            % Nested function: Find the root of a function using the bisection method
            function root = find_root(func, a, b, tol)
                % Bisection method to find the root of func in the interval [a, b]
                fa = func(a);
                fb = func(b);
        
                if fa * fb > 0
                    error('Function does not change sign in the interval [%f, %f]', a, b);
                end
        
                while (b - a) / 2 > tol
                    c = (a + b) / 2;
                    fc = func(c);
                    if fc == 0
                        break;
                    elseif fa * fc < 0
                        b = c;
                        fb = fc;
                    else
                        a = c;
                        fa = fc;
                    end
                end
        
                root = (a + b) / 2;
            end
        
            % Nested function: Derivative of Legendre polynomial of degree n
            function val = deriv_legendre_poly(n, x)
                % Compute the derivative of the Legendre polynomial of degree n at x
                if n == 0
                    val = 0;
                elseif n == 1
                    val = 1;
                else
                    Pn = hw06.legendre_poly(n, x);
                    Pn_minus_1 = hw06.legendre_poly(n - 1, x);
                    val = n / (x^2 - 1) * (x * Pn - Pn_minus_1);
                end
            end
        end


        function ret = p5(n)
            % Construct the Gauss quadrature rule using the roots of the Legendre polynomial of degree n.
        
            % Initialize arrays for roots and weights
            roots = zeros(n, 1);
            weights = zeros(n, 1);
        
            % Step 1: Initial guesses for the roots using cosine formula
            initial_guesses = cos((4 * (1:n) - 1) * pi / (4 * n + 2));
        
            % Step 2: Find the roots using Newton's method
            for k = 1:n
                roots(k) = newton_method(@(x) hw06.legendre_poly(n, x), ...
                                         @(x) hw06.deriv_lengendre_poly(n, x), ...
                                         initial_guesses(k), 1e-14, 100);
            end
        
            % Step 3: Compute weights
            for k = 1:n
                weights(k) = 2 / ((1 - roots(k)^2) * (hw06.deriv_lengendre_poly(n, roots(k))^2));
            end
        
            % Combine roots and weights into the result
            ret = [roots, weights];
        
            % Nested function: Newton's method for finding roots
            function root = newton_method(f, df, x0, tol, max_iter)
                % f: function handle for the function
                % df: function handle for the derivative of the function
                % x0: initial guess
                % tol: tolerance for convergence
                % max_iter: maximum number of iterations
        
                x = x0;
                for iter = 1:max_iter
                    fx = f(x);
                    dfx = df(x);
                    if abs(fx) < tol
                        break;
                    end
                    if dfx == 0
                        error('Derivative is zero. Newton''s method failed.');
                    end
                    x = x - fx / dfx;
                end
        
                root = x;
            end
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                                                             %
        % Helper functions below. Do not modify. You can create your own helper functions if needed.                  %
        %                                                                                                             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Helper functions for p4. The following function is used to evaluate the Legendre polynomial of degree 6.
        function val = legendre_poly_6(x)
            % Compute the Legendre polynomial of degree 6 at the point x.
            %
            % :param x: The point at which to evaluate the Legendre polynomial.
            % :return: The value of the Legendre polynomial of degree 6 at the point x.

            val = (231 * x^6 - 315 * x^4 + 105 * x^2 - 5) / 16;
        end

        % Helper functions for p5. The following function is used to evaluate the Legendre polynomial of degree n.
        function val = legendre_poly(n, x)
            % Compute the nth Legendre polynomial P_n at the point x.
            %
            % :param n: The degree of the Legendre polynomial.
            % :param x: The point at which to evaluate the Legendre polynomial.
            % :return: The value of the nth Legendre polynomial at the point x.

            if (n == 0)
                val = 1;
            elseif (n == 1)
                val = x;
            else
                val = hw06.legendre_poly(n-1, x) * x * (2 * n - 1)/n - (n - 1) * hw06.legendre_poly(n - 2, x) / n;
            end
        end

        function val = deriv_lengendre_poly(n, x)
            % Compute the derivative of the nth Legendre polynomial P_n at the point x.
            %   
            % :param n: The degree of the Legendre polynomial.
            % :param x: The point at which to evaluate the derivative of the Legendre polynomial.
            % :return: The value of the derivative of the nth Legendre polynomial at the point x.
            val = n / (x^2 - 1) * (x * hw06.legendre_poly(n, x) - hw06.legendre_poly(n - 1, x));
        end
    end
end
