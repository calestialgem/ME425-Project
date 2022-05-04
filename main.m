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
I = 100/n;
% Torsional Stiffness Between Disks
k = 25*n;
% Inertia Matrix
M = zeros(n, n);
for j = 1:n
    M(j, j) = I;
end
% Stiffness Matrix
K = zeros(n, n);
for j = 1:n
    if j > 1
        K(j, j-1) = -k;
    end
    if j < n
        K(j, j+1) = -k;
        K(j, j) = 2*k;
    else
        K(j, j) = k;
    end
end

file.print("");
file.print("Initial Caclculations:");
file.print("~~~~~~~~~~~~~~~~~~~~~~");
file.print("[-] %1s = %5.1f %0s", "I", I, "");
file.print("[-] %1s = %5.1f %0s", "k", k, "");
file.prmat("[-] M", M, "%5.0f");
file.prmat("[-] K", K, "%5.0f");

% PART A ----------------------------------------------------------------------

% First Transformation
M_ = M^(-1/2);
K_ = M_*K*M_;

% Eigenvector and Eigenvalue Matrix
[P, L] = eig(K_);

% Natural Frequencies in rad/s
w = zeros(n, 1);
for j = 1:n
    w(j) = L(j, j)^(1/2);
end

file.print("");
file.print("Part A:");
file.print("~~~~~~~");
file.prmat("[-] M_", M_, "%5.1f");
file.prmat("[-] K_", K_, "%5.1f");

for j = 1:n
    file.print("[*] w_%1.0f = %5.3f %5s", j, w(j), "rad/s");
    file.prvec(sprintf("[*] v_%1.0f", j), P(j, :), "%5.1f");
end

% PART B ----------------------------------------------------------------------

% Force Vector (Divided by exp(iwt))
F = zeros(n, 1);
F(1) = k;

w_range = max(w)*(0:0.001:1.5);
T_n_range = zeros(size(w_range));
for j = 1:length(w_range)
    % Base Excitation Frequency
    w = w_range(j);
    % Angular Displacement Vector (Divided by exp(iwt))
    T_ = F'/(K-w^2*M);
    T_n_range(j) = T_(n);
end

figure();
set(gca, 'YScale', 'log');
set(gca, 'XScale', 'log');
hold('on');
grid('on');
xlabel('\omega');
ylabel('|\Theta_n/\Phi|');
plot(w_range, abs(T_n_range), 'LineWidth', 2);
saveas(gcf, 'Transmissibility', 'jpeg');

% PART C ----------------------------------------------------------------------

Ia = u/2;

file.print("");
file.print("Part C:");
file.print("~~~~~~~");
file.print("[-] %2s = %5.3f %0s", "Ia", Ia, "");
