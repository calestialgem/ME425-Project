% gecgelcem 02.05.2022
% me425 spring2022 prj

close('all');
clear();
clc();

file = printer('Output.txt');
file.print("gecgelcem 02.05.2022");
file.print("me425 spring2022 prj");

% Number of Disks
n_min = 2;
n_max = 5;
n = round(ask(sprintf("\nEnter %1s in [%1.0f, %1.0f] %0s: ", "n", n_min, n_max, ""), n_min, n_max));
% Total Houdaille Damper Viscosity
u_min = 0.1;
u_max = 0.3;
u = ask(sprintf("\nEnter %1s in [%3.1f, %3.1f] %0s: ", "u", u_min, u_max, ""), u_min, u_max);

file.print("");
file.print("Input:");
file.print("~~~~~~");
file.print("[-] %1s = %4.0f %0s", "n", n, "");
file.print("[-] %1s = %4.2f %0s", "u", u, "");

% Rotational Inertia of a Disk
I = 100 / n;
% Torsional Stiffness Between Disks
k = 25 * n;

file.print("");
file.print("Initial Caclculations:");
file.print("~~~~~~~~~~~~~~~~~~~~~~");
file.print("[-] %1s = %5.1f %0s", "I", I, "");
file.print("[-] %1s = %5.1f %0s", "k", k, "");

% PART A ----------------------------------------------------------------------

% Timer Start
c_start = tic();
% Inertia Matrix
M = zeros(n);
for j = 1:n
    M(j, j) = I;
end
% Stiffness Matrix
K = zeros(n);
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
% First Transformation
M_ = M^(-1/2);
K_ = M_ * K * M_;
% Eigenvector and Eigenvalue Matrix
[P, L] = eig(K_);
% Natural Frequencies in rad/s
w = zeros(n, 1);
for j = 1:n
    w(j) = L(j, j)^(1/2);
end
% Elapsed Time
c_elapsed = toc(c_start);

file.print("");
file.print("Part A:");
file.print("~~~~~~~");
file.print("[-] Elapsed Time: %5.1f s", c_elapsed);
file.prmat("[-] M", M, "%5.0f");
file.prmat("[-] K", K, "%5.0f");
file.prmat("[-] M_", M_, "%5.1f");
file.prmat("[-] K_", K_, "%5.1f");
for j = 1:n
    file.print("[*] w_%1.0f = %5.3f %5s", j, w(j), "rad/s");
    file.prvec(sprintf("[*] v_%1.0f", j), P(j, :), "%5.1f");
end

% PART B ----------------------------------------------------------------------

% Timer Start
c_start = tic();
% Maximum Excitation Frequency
w_e_max = max(w) * 1.5;
% Excitation Frequency Range
w_e_range = max(w) * (10.^(-1:0.001:log10(1.5)));
% Plot
plot_T_range(n, w_e_range, M, f_C(n, 0, [], []), K, k, 'Part B Transmissibility');
% Elapsed Time
c_elapsed = toc(c_start);

file.print("");
file.print("Part B:");
file.print("~~~~~~~");
file.print("[-] Elapsed Time: %5.1f s", c_elapsed);

% PART C ----------------------------------------------------------------------

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
x_f = @(x) f_T_max(n, w_e_range, M, f_C(n, m, na, f_ca(x)), K, k);
% Optimization Options
x_options = optimoptions('fminimax');
x_options.MaxIterations = 100;
x_options.MaxFunctionEvaluations = 1000;
x_options.Display = 'iter';
% Optimization Results
[x, ~, ~] = fminimax(x_f, x_0, [], [], [], [], x_lb, [], [], x_options);
% Absorber Dampings
ca = f_ca(x);
% Damping Matrix
C = f_C(n, m, na, ca);
% Plot
plot_T_range(n, w_e_range, M, C, K, k, 'Part C Transmissibility');
% Elapsed Time
c_elapsed = toc(c_start);

file.print("");
file.print("Part C:");
file.print("~~~~~~~");
file.print("[-] Elapsed Time: %5.1f s", c_elapsed);
file.prvec("[-] Ia", Ia, "%7.3f");
file.prvec("[-] na", na, "%7.0f");
file.prmat("[-] M", M, "%7.1f");
file.prmat("[-] K", K, "%7.1f");
file.prvec("[-] x_0", x_0, "%7.3f");
file.prvec("[-] x", x, "%7.3f");
file.prvec("[-] ca", ca, "%7.3f");
file.prmat("[-] C", C, "%7.1f");

% PART D ----------------------------------------------------------------------

% Timer Start
c_start = tic();
% Number of Absorbers
m = 2;
% Absorber Positions (Assumed to be unique for each absorber.)
na = zeros(m, 1);
% Absorber Inertias
Ia = zeros(m, 1);
% Absorber Dampings
ca = zeros(m, 1);
% Minimum Maximum Transmissibility
T_minimax = Inf;
% For all possible combinations...
for na1 = 1:n - 1
    for na2 = na1 + 1:n
        % Absorber Positions (Assumed to be unique for each absorber.)
        na_j = zeros(m, 1);
        na_j(1) = na1;
        na_j(2) = na2;
        % Optimize
        [Ia_j, ca_j, T_minimax_j] = f_T_minimax(n, w_e_range, m, na_j, u, I, k);
        % ... replace if better.
        if T_minimax > T_minimax_j
            na = na_j;
            Ia = Ia_j;
            ca = ca_j;
            T_minimax = T_minimax_j;
        end
    end
end
if ~isinf(T_minimax)
    % Inertia Matrix
    M = f_M(n, m, I, Ia);
    % Damping Matrix
    C = f_C(n, m, na, ca);
    % Stiffness Matrix
    K = f_K(n, m, k);
    % Plot
    plot_T_range(n, w_e_range, M, C, K, k, 'Part D Transmissibility');
end
% Elapsed Time
c_elapsed = toc(c_start);

file.print("");
file.print("Part D:");
file.print("~~~~~~~");
file.print("[-] Elapsed Time: %5.1f s", c_elapsed);
file.prvec("[-] na", na, "%7.0f");
file.prvec("[-] Ia", Ia, "%7.3f");
file.prvec("[-] ca", ca, "%7.3f");
if ~isinf(T_minimax)
    file.prmat("[-] M", M, "%7.1f");
    file.prmat("[-] C", C, "%7.1f");
    file.prmat("[-] K", K, "%7.1f");
else
    file.print("[!] Could not found even a single finite solution!");
end

% Inertia Matrix
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

% Damping Matrix
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

% Stiffness Matrix
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

function [Ia, ca, T_minimax] = f_T_minimax(n, w_e_range, m, na, u, I, k)
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
    x_f = @(x) f_T_max(n, w_e_range, f_M(n, m, I, f_Ia(x)), f_C(n, m, na, f_ca(x)), K, k);
    % Linear Equality Constraint
    x_Aeq = [1, 1, 0, 0];
    x_beq = u;
    % Optimization Options
    x_options = optimoptions('fminimax');
    x_options.MaxIterations = 100;
    x_options.MaxFunctionEvaluations = 1000;
    x_options.Display = 'iter';
    % Optimization Results
    [x, ~, x_maxfval] = fminimax(x_f, x_0, [], [], x_Aeq, x_beq, x_lb, [], [], x_options);
    % Absorber Inertias
    Ia = f_Ia(x);
    % Absorber Dampings
    ca = f_ca(x);
    % Minimized Maximum Transmissibility
    T_minimax = x_maxfval;
end

function T_range = f_T_max(n, w_e_range, M, C, K, k)
    % Transmissibility Range
    T_range = zeros(size(w_e_range));
    for j = 1:length(w_e_range)
        T_range(j) = f_T(n, w_e_range(j), M, C, K, k);
    end
end

% Transmissibility
function T = f_T(n, w_e, M, C, K, k)
    % Receptance Matrix
    a = (-w_e^2 * M + w_e * C + K)^ - 1;
    % Transmissibility
    T = abs(a(1, n) * k);
end

% Plot Transmisibility Range
function plot_T_range(n, w_e_range, M, C, K, k, title)
    % Plot
    figure();
    set(gca, 'YScale', 'log');
    set(gca, 'XScale', 'log');
    hold('on');
    grid('on');
    xlabel('\omega');
    ylabel('|\Theta_n/\Phi|');
    plot(w_e_range, f_T_max(n, w_e_range, M, C, K, k), 'LineWidth', 2);
    saveas(gcf, title, 'jpeg');
end
