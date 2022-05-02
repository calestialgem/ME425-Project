% gecgelcem 02.05.2022
% me425 spring2022 prj

close('all');
clear();
clc();

file = printer('Output.txt');
file.print("gecgelcem 02.05.2022");
file.print("me425 spring2022 prj");

n = ask(sprintf("\nEnter %1s in [%1.0f, %1.0f] %0s: ", "n", 2, 5, ""), 2, 5);

file.print("");
file.print("Input:");
file.print("~~~~~~");
file.print("[-] %1s = %1.0f %3s", "n", n, "");
