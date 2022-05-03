% gecgelcem 02.05.2022
% me425 spring2022 prj

close('all');
clear();
clc();

file = printer('Output.txt');
file.print("gecgelcem 02.05.2022");
file.print("me425 spring2022 prj");

% Number of Disks
n = ask(sprintf("\nEnter %1s in [%1.0f, %1.0f] %0s: ", "n", 2, 5, ""), 2, 5);

file.print("");
file.print("Input:");
file.print("~~~~~~");
file.print("[-] %1s = %1.0f %0s", "n", n, "");

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
    if j > 1 && j < n
        K(j, j) = 2*k;
    else
        K(j, j) = k;
    end
    if j > 1
        K(j, j-1) = -k;
    end
    if j < n
        K(j, j+1) = -k;
    end
end

file.print("");
file.print("Initial Caclculations:");
file.print("~~~~~~~~~~~~~~~~~~~~~~");
file.print("[-] %1s = %5.3f %0s", "I", I, "");
file.print("[-] %1s = %5.3f %0s", "k", k, "");
file.prmat("[-] M", M, "%5.0f");
file.prmat("[-] K", K, "%5.0f");

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
