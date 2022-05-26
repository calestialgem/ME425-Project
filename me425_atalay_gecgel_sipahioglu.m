% me425 spring2022 project
% atalay gecgel sipahioglu

close('all');
clear();
clc();

print("me425 spring2022 project");
print("atalay gecgel sipahioglu");

% ------------------------------------------------------------------------------
% INPUT ------------------------------------------------------------------------
% ------------------------------------------------------------------------------

% Get the n and u from the user.

% Number of Disks (Min: 2, Max: 5)
n = 4;

% Total Houdaille Damper Viscosity (Min: 0.1, Max: 0.3)
u = 0.1;

% Print
print("");
print("Input:");
print("~~~~~~");
print("[-] n = %4.0f", n);
print("[-] u = %4.2f", u);

% ------------------------------------------------------------------------------
% INITIAL CALCULATIONS ---------------------------------------------------------
% ------------------------------------------------------------------------------

% Find I and k.
% Note: using `n`, `u` from the previous part.

% Rotational Inertia of a Disk
I = 100 / n;

% Torsional Stiffness Between Disks
k = 25 * n;

% Print
print("");
print("Initial Caclculations:");
print("~~~~~~~~~~~~~~~~~~~~~~");
print("[-] I = %5.1f", I);
print("[-] k = %5.1f", k);

% ------------------------------------------------------------------------------
% PART A -----------------------------------------------------------------------
% ------------------------------------------------------------------------------

% Find M and K.
% Use modal analysis to find the natural frequencies and mode shapes.
% Note: using `n`, `I`, `k` from the previous part.
% Note: using functions `f_M`, `f_K` that are defined at the end.

% Timer Start
c_start = tic();

% Inertia Matrix
M = f_M(n, 0, I, []);

% Stiffness Matrix
K = f_K(n, 0, k);

% First Transformation
M_ = M^(-1/2);
K_ = M_ * K * M_;

% Eigenvector and Eigenvalue Matrix
% Eigenvectors are the mode shapes.
% Eigenvalues are the squares of the natural frequencies.
[P, L] = eig(K_);

% Natural Frequencies in rad/s
w = zeros(n, 1);
for j = 1:n
    % Natural frequency is found by taking the square root of the eigenvalues.
    w(j) = L(j, j)^(1/2);
    % Mode shapes are found by making the first element of the eigenvectors positive.
    P(:, j) = P(1, j) / abs(P(1, j)) * P(:, j);
end

% Plot
figure();
tiledlayout(n, 1);
xlabel('n');
for j = 1:n
    nexttile();
    hold('on');
    grid('on');
    title(sprintf("\\omega_%.0f=%5.3f rad/s", j, w(j)));
    ylim([-1 1]);
    plot(1:n, P(:, j), '-o', 'LineWidth', 2);
    plot(1:n, zeros(1, n), 'k--', 'LineWidth', 2);
end
export_plot("Part A Mode Shapes", n, u, 0.9, 0.9);

% Elapsed Time
c_elapsed = toc(c_start);

% Print
print("");
print("Part A:");
print("~~~~~~~");
print("[-] Elapsed Time: %5.2f s", c_elapsed);
prmat("[-] M", M, "%9.3f");
prmat("[-] K", K, "%9.3f");
prmat("[-] M_", M_, "%9.3f");
prmat("[-] K_", K_, "%9.3f");
for j = 1:n
    print("[*] w_%1.0f = %5.3f rad/s", j, w(j));
    prvec(sprintf("[*] v_%1.0f", j), P(:, j), "%5.1f");
end

% ------------------------------------------------------------------------------
% PART B -----------------------------------------------------------------------
% ------------------------------------------------------------------------------

% Construct the range of excitation frequencies and find the transmissibilities
% using the repectance matrix.
% Note: using `w`, `n`, `M`, `K`, `k` from the previous part.
% Note: using functions `f_C`, `f_T_peaks` that are defined at the end.

% Timer Start
c_start = tic();

% Excitation Frequency Range
w_e_range = max(w) * (10.^(-1:0.001:log10(1.5)));

% Damping Matrix
% This is just a zero matrix in this part.
C = f_C(n, 0, [], []);

% Plot Transmissibility Range
plot_T_range(n, u, w_e_range, M, C, K, k, "Part B Transmissibility");

% Elapsed Time
c_elapsed = toc(c_start);

% Print
print("");
print("Part B:");
print("~~~~~~~");
print("[-] Elapsed Time: %5.2f s", c_elapsed);
print("[-] T_max = %.3f", max(f_T_peaks(n, w, M, C, K, k)));

% ------------------------------------------------------------------------------
% PART C -----------------------------------------------------------------------
% ------------------------------------------------------------------------------

% Optimize the maximum of the peaks in the transmissibility for the damping
% coefficients of the absorbers which are connected to 1 and 5.
% Note: using `u`, `n`, `I`, `k`, `w`, `w_e_range` from the previous part.
% Note: using functions `f_M`, `f_K`, `f_T_peaks`, `f_C` that are defined at the
% end.

% Timer Start
c_start = tic();

% Number of Absorbers
m = 2;

% Absorber Inertias
Ia = zeros(m, 1);
Ia(1) = u / 2;
Ia(2) = u / 2;

% Absorber Positions (Assumed to be unique for each absorber.)
na = zeros(m, 1);
na(1) = 1;
na(2) = n;

% Inertia Matrix
M = f_M(n, m, I, Ia);

% Stiffness Matrix
K = f_K(n, m, k);

% Optimization Parameter Vector
% [ca1, ca2]
f_ca = @(x) [x(1); x(2)];

% Initial Value
x_0 = [0.5, 0.5];

% Lower Bound
x_lb = [0, 0];

% Optimized Function
x_f = @(x) f_T_peaks(n, w, M, f_C(n, m, na, f_ca(x)), K, k);

% Optimization Options
x_options = optimoptions('fminimax');
x_options.MaxIterations = 100;
x_options.MaxFunctionEvaluations = 1000;
x_options.Display = 'off';

% Optimization Results
[x, ~, x_maxfval] = fminimax(x_f, x_0, [], [], [], [], x_lb, [], [], x_options);

% Absorber Dampings
ca = f_ca(x);

% Minimized Maximum Transmissibility
T_min = x_maxfval;

% Damping Matrix
C = f_C(n, m, na, ca);

% Plot Transmissibility Range
plot_T_range(n, u, w_e_range, M, C, K, k, "Part C Transmissibility");

% Elapsed Time
c_elapsed = toc(c_start);

% Print
print("");
print("Part C:");
print("~~~~~~~");
print("[-] Elapsed Time: %5.2f s", c_elapsed);
prvec("[-] Ia", Ia, "%7.3f");
prvec("[-] na", na, "%7.0f");
prmat("[-] M", M, "%9.3f");
prmat("[-] K", K, "%9.3f");
prvec("[-] x_0", x_0, "%7.3f");
prvec("[-] x", x, "%7.3f");
prvec("[-] ca", ca, "%7.3f");
prmat("[-] C", C, "%9.3f");
print("[-] T_max = %.3f", T_min);

% ------------------------------------------------------------------------------
% PART D -----------------------------------------------------------------------
% ------------------------------------------------------------------------------

% Optimize all the possible combinations of the absorber positions over the
% damping coefficients and absorber inertias.
% Select the one with the smallest transmissibility.
% Note: using `m`, `n`, `w`, `u`, `I`, `k`, `w_e_range` from the previous part.
% Note: using functions `f_T_min`, `f_M`, `f_C`, `f_K` that are defined at the
% end.

% Timer Start
c_start = tic();

% Absorber Positions (Assumed to be unique for each absorber.)
na = zeros(m, 1);

% Absorber Inertias
Ia = zeros(m, 1);

% Absorber Dampings
ca = zeros(m, 1);

% Minimum Maximum Transmissibility
T_min = Inf;

% For all possible combinations...
na_combinations = nchoosek(1:n, m);
for j = size(na_combinations, 1)
    na_j = na_combinations(j, :);

    % Optimize
    [Ia_j, ca_j, T_min_j] = f_T_min(n, m, w, na_j, u, I, k);

    % ... replace if better.
    if T_min > T_min_j
        na = na_j;
        Ia = Ia_j;
        ca = ca_j;
        T_min = T_min_j;
    end
end

if ~isinf(T_min)
    % Inertia Matrix
    M = f_M(n, m, I, Ia);

    % Damping Matrix
    C = f_C(n, m, na, ca);

    % Stiffness Matrix
    K = f_K(n, m, k);

    % Plot Transmissibility Range
    plot_T_range(n, u, w_e_range, M, C, K, k, "Part D Transmissibility");
end

% Elapsed Time
c_elapsed = toc(c_start);

% Print
print("");
print("Part D:");
print("~~~~~~~");
print("[-] Elapsed Time: %5.2f s", c_elapsed);
prvec("[-] na", na, "%7.0f");
prvec("[-] Ia", Ia, "%7.3f");
prvec("[-] ca", ca, "%7.3f");
if ~isinf(T_min)
    prmat("[-] M", M, "%9.3f");
    prmat("[-] C", C, "%9.3f");
    prmat("[-] K", K, "%9.3f");
    print("[-] T_max = %.3f", T_min);
else
    print("[!] Could not found even a single finite solution!");
end

% ------------------------------------------------------------------------------
% CONSTRUCTION FUNCTIONS -------------------------------------------------------
% ------------------------------------------------------------------------------

% For creating the inertia matrix. In the case with no absorbers give `m` as
% `0`, `Ia` as `[]`.
function M = f_M(n, m, I, Ia)
    % Inertia Matrix
    M = zeros(n + m);
    for j = 1:n
        M(j, j) = I;
    end
    for j = 1:m
        M(n + j, n + j) = Ia(j);
    end
end

% For creating the damping matrix. In the case with no absorbers give `m` as
% `0`, `na` and `ca` as `[]`. The resulting matrix will be just zeros.
function C = f_C(n, m, na, ca)
    % Damping Matrix
    C = zeros(n + m);
    for j = 1:m
        C(n + j, n + j) = ca(j);
        C(n + j, na(j)) = -ca(j);
        C(na(j), n + j) = -ca(j);
        C(na(j), na(j)) = ca(j);
    end
end

% For creating the stiffness matrix. In the case with no absorbers give `m` as
% `0`.
function K = f_K(n, m, k)
    % Stiffness Matrix
    K = zeros(n + m);
    for j = 1:n
        if j > 1
            K(j, j - 1) = -k;
        end
        if j < n
            K(j, j + 1) = -k;
            K(j, j) = 2 * k;
        else
            K(j, j) = k;
        end
    end
end

% ------------------------------------------------------------------------------
% TRANSMISSIBILITY FUNCTIONS ---------------------------------------------------
% ------------------------------------------------------------------------------

% For designing the absorbers in part D. Minimizes the maximum of the peaks
% using `fminimax` and returns the optimum parameters. Only needs to know the
% positions of the absorbers. The combination of all possible absorber positions
% can be iterated over to get the best desing. The total absorber inertia is
% equated to viscosity by giving a linear equality constraint to `fminimax`.
function [Ia, ca, T_min] = f_T_min(n, m, w, na, u, I, k)
    % Stiffness Matrix
    K = f_K(n, m, k);

    % Optimization Parameter Vector
    % [Ia1, Ia2, ca1, ca2]
    f_Ia = @(x) [x(1); x(2)];
    f_ca = @(x) [x(3); x(4)];

    % Initial Value
    x_0 = [u / 2, u / 2, 0.5, 0.5];

    % Lower Bound
    x_lb = [0, 0, 0, 0];

    % Optimized Function
    x_f = @(x) f_T_peaks(n, w, f_M(n, m, I, f_Ia(x)), f_C(n, m, na, f_ca(x)), K, k);

    % Linear Equality Constraint
    x_Aeq = [1, 1, 0, 0];
    x_beq = u;

    % Optimization Options
    x_options = optimoptions('fminimax');
    x_options.MaxIterations = 100;
    x_options.MaxFunctionEvaluations = 1000;
    x_options.Display = 'off';

    % Optimization Results
    [x, ~, x_maxfval] = fminimax(x_f, x_0, [], [], x_Aeq, x_beq, x_lb, [], [], x_options);

    % Absorber Inertias
    Ia = f_Ia(x);

    % Absorber Dampings
    ca = f_ca(x);

    % Minimized Maximum Transmissibility
    T_min = x_maxfval;
end

% For finding the peaks in the transmissibility plot very fast and secure. It
% always returns the peaks because it uses the natural frequencies. With the
% damping the peaks shift to left slightly. Looks for the peaks via iterating
% over the natural frequencies and using `fminbnd` starting from the 95% of the
% natural frequency to 100% of it. Fast because it looks for a very small
% window. Guaranteed to find the peaks because it does not do grid search.
function T_peaks = f_T_peaks(n, w, M, C, K, k)
    % Transmissibility Peaks
    T_peaks = zeros(1, n);
    for j = 1:n
        % Optimization Parameter Vector
        % [w_e]
        % Lower Bound
        x_lb = w(j) * 0.95;

        % Upper Bound
        x_ub = w(j);

        % Optimized Function
        x_f = @(x) -f_T(n, x, M, C, K, k);

        % Optimization Options
        x_options = optimset('fminbnd');
        x_options.MaxIterations = 100;
        x_options.MaxFunctionEvaluations = 1000;
        x_options.Display = 'off';

        % Optimization Results
        [~, x_fval] = fminbnd(x_f, x_lb, x_ub, x_options);

        % Maximum Transmissibility
        T_peaks(j) = -x_fval;
    end
end

% For calculating a range of transmissibilities for the given range of
% excitation frequencies.
function T_range = f_T_range(n, w_e_range, M, C, K, k)
    % Transmissibility Range
    T_range = zeros(size(w_e_range));
    for j = 1:length(w_e_range)
        T_range(j) = f_T(n, w_e_range(j), M, C, K, k);
    end
end

% For calculating the transmissibility of the last disk for the given excitation
% frequency. Very fast because it only uses the necessary element of the
% receptance matrix.
function T = f_T(n, w_e, M, C, K, k)
    % Receptance Matrix
    a = (-w_e^2 * M + 1i * w_e * C + K)^ - 1;

    % Transmissibility
    T = abs(a(1, n) * k);
end

% ------------------------------------------------------------------------------
% OUTPUT FUNCTIONS -------------------------------------------------------------
% ------------------------------------------------------------------------------

% For plotting the transmissibility over a range of excitation frequencies.
function plot_T_range(n, u, w_e_range, M, C, K, k, name)
    % Plot
    figure();
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    hold('on');
    grid('on');
    xlabel('\omega (rad/s)');
    ylabel('|\Theta_n/\Phi|');
    plot(w_e_range, f_T_range(n, w_e_range, M, C, K, k), 'LineWidth', 2);
    title(name);
    export_plot(name, n, u, 0.9, 0.5);
end

% For easier general printing. Puts new line at the start.
function print(varargin)
    fprintf('%s\n', sprintf(varargin{:}));
end

% For printing matrices.
function prmat(name, matrix, element)
    print("%s [%.0f, %.0f]: ", name, size(matrix, 1), size(matrix, 2));
    for k = 1:size(matrix, 1)
        for j = 1:size(matrix, 2)
            fprintf(element, matrix(k, j));
        end
        fprintf("\n");
    end
end

% For printing vectors.
function prvec(name, vector, element)
    fprintf("%s [%.0f]: ", name, length(vector));
    for k = 1:length(vector)
        fprintf(element, vector(k));
    end
    fprintf("\n");
end

function export_plot(name, n, u, w, h)
    matlab2tikz(sprintf('%s n=%.0f u=%.2f.tikz', name, n, u), ...
        'height', sprintf('%.2f\\textwidth', h), ...
        'width', sprintf('%.2f\\textwidth', w), ...
        'extraAxisOptions', ...
        {'ylabel style={font=\small}', ...
        'xlabel style={font=\small}'}, ...
        'showInfo', false);
end
