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
% Base Excitation (Divided by exp(iwt))
P_ = 1;
% Force Vector (Divided by exp(iwt))
F = zeros(n, 1);
F(1) = k * P_;
% Modal Displacement Transformation
S = M_ * P;
% Modal Forcing (Divided by exp(iwt))
f = P' * M_ * F;
% Excitation Frequency Range
w_e_range = max(w) * (0:1.5e-3:1.5);
T_n_range = zeros(size(w_e_range));
for j = 1:length(w_e_range)
    % Base Excitation Frequency
    w_e = w_e_range(j);
    % Normalized Excitation Frequencies
    r = w_e ./ w;
    % Modal Displacement Vector (Divided by exp(iwt))
    R = (f ./ w.^2) ./ (1 - r.^2);
    % Displacement Vector (Divided by exp(iwt))
    T_ = S * R;
    T_n_range(j) = T_(n);
end
% Plot
figure();
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
hold('on');
grid('on');
xlabel('\omega');
ylabel('|\Theta_n/\Phi|');
plot(w_e_range, abs(T_n_range) / P_, 'LineWidth', 2);
for j = 1:n
    xline(w(j), '--');
end
saveas(gcf, 'Transmissibility', 'jpeg');
% Elapsed Time
c_elapsed = toc(c_start);

file.print("");
file.print("Part B:");
file.print("~~~~~~~");
file.print("[-] Elapsed Time: %5.1f s", c_elapsed);
file.print("[-] %2s = %5.3f %0s", "P_", P_, "");
file.prvec("[-] F", F, "%7.2f");
file.prvec("[-] f", f, "%7.2f");
file.prmat("[-] S", S, "%7.2f");

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
M = zeros(n + m);
for j = 1:n
    M(j, j) = I;
end
for j = 1:m
    M(n + j, n + j) = Ia(j);
end
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
% First Transformation
M_ = M^(-1/2);
K_ = M_ * K * M_;
% Eigenvector and Eigenvalue Matrix
[P, L] = eig(K_);
% Natural Frequencies in rad/s
w = zeros(n + m, 1);
for j = 1:n + m
    w(j) = L(j, j)^(1/2);
end
% Force Vector (Divided by exp(iwt))
F = zeros(n + m, 1);
F(1) = k * P_;
% Modal Displacement Transformation
S = M_ * P;
% Modal Forcing (Divided by exp(iwt))
f = P' * M_ * F;
% Excitation Frequency Limit
w_e_max = max(w) * 1.5;
% Optimization Parameter Vector
% [w_e, ca1, ca2]
% Initial Value
x_0 = [max(w), 40, 10];
% Lower Bound
x_lb = [0, 0.1, 0.1];
% Upper Bound
x_ub = [w_e_max, 100, 100];
% Optimized Function
x_f = @(x) f_T(n, m, na, [x(2); x(3)], w, x(1), M, K, f, S, P_);
% Optimization Options
x_options = optimoptions('fminimax');
x_options.MaxIterations = 2000;
x_options.MaxFunctionEvaluations = 2000;
% Optimization Results
[x, ~, ~, flag, output] = fminimax(x_f, x_0, [], [], [], [], x_lb, x_ub, [], x_options);
% Absorber Dampings
ca = [x(2); x(3)];
% Damping Matrix
C = zeros(n + m);
for j = 1:m
    C(n + j, n + j) = ca(j);
    C(n + j, na(j)) = -ca(j);
    C(na(j), n + j) = -ca(j);
    C(na(j), na(j)) = ca(j);
end
% Elapsed Time
c_elapsed = toc(c_start);

file.print("");
file.print("Part C:");
file.print("~~~~~~~");
file.print("[-] Elapsed Time: %5.1f s", c_elapsed);
if flag == 0
    file.print("[!] Number of iterations exceeded options.MaxIterations or the number of function evaluations exceeded options.MaxFunctionEvaluations!");
elseif flag == 4
    file.print("[!] Magnitude of the search direction was less than the specified tolerance, and the constraint violation was less than options.ConstraintTolerance!");
elseif flag == 5
    file.print("[!] Magnitude of the directional derivative was less than the specified tolerance, and the constraint violation was less than options.ConstraintTolerance!");
elseif flag == -1
    file.print("[!] Stopped by an output function or plot function!");
elseif flag == -2
    file.print("[!] No feasible point was found!");
elseif flag ~= 1
    file.print("[!] An unknown error occured!");
end
file.print("[?] Optimizer Output: \n%s", output.message);
file.prvec("[-] Ia", Ia, "%7.3f");
file.prvec("[-] na", na, "%7.0f");
file.prmat("[-] M", M, "%7.1f");
file.prmat("[-] K", K, "%7.1f");
file.prmat("[-] M_", M_, "%7.1f");
file.prmat("[-] K_", K_, "%7.1f");
for j = 1:n + m
    file.print("[*] w_%1.0f = %5.3f %5s", j, w(j), "rad/s");
    file.prvec(sprintf("[*] v_%1.0f", j), P(j, :), "%5.1f");
end
file.prvec("[-] x_0", x_0, "%7.3f");
file.prvec("[-] x_lb", x_lb, "%7.3f");
file.prvec("[-] x_ub", x_ub, "%7.3f");
file.prvec("[-] x", x, "%7.3f");
file.prvec("[-] ca", ca, "%7.3f");
file.prmat("[-] C", C, "%7.1f");

% Transmissibility
function T = f_T(n, m, na, ca, w, w_e, M, K, f, S, P_)
    % Damping Ratios
    z = f_z(n, m, na, ca, w, M, K);
    % Normalized Excitation Frequencies
    r = w_e ./ w;
    % Modal Displacement Vector (Divided by exp(iwt))
    R = (f ./ w.^2) ./ (((1 - r.^2) ./ (2 .* z .* r)).^2 + 1).^(1/2);
    % Displacement Vector (Divided by exp(iwt))
    T_ = S * R;
    % Transmissibility
    T = abs(T_(n)) / P_;
end

% Rayleigh Damping Approximation
function z = f_z(n, m, na, ca, w, M, K)
    % Real Damping Matrix
    C = zeros(n + m);
    for j = 1:m
        C(n + j, n + j) = ca(j);
        C(n + j, na(j)) = -ca(j);
        C(na(j), n + j) = -ca(j);
        C(na(j), na(j)) = ca(j);
    end
    % Root-Mean-Square Difference of Rayleigh Damping Approximation
    rms = @(a, b) sqrt(sum(sum((C - (a * M + b * K)).^2)));
    % Optimization Parameter Vector
    % [a, b]
    % Initial Value
    x_0 = [0.5, 0.1];
    % Optimized Function
    x_f = @(x) rms(x(1), x(2));
    % Optimization Results
    x = fminsearch(x_f, x_0);
    % Rayleigh Damping
    a = x(1);
    b = x(2);
    % Damping Ratios
    z = w / a / 2 + w * b / 2;
end
