clear all; clc; 
syms L12 L21 G1 G2 P1 P2 T
T_arr = zeros(20,1);
y = zeros(20,1);
x = [0.9809, 0.9707, 0.9629, 0.9503, 0.9162, 0.8961, 0.8896, 0.8016, 0.7761, 0.6903, 0.6622, 0.5544, 0.5416, 0.4633, 0.3910, 0.3512, 0.2743, 0.2402, 0.1366, 0.1097];
A12 = 1098.7901;
A21 = 1266.6728;
A1=8.07131; B1=1730.630; C1=233.426;
A2=7.04115; B2=1373.799; C2=214.979;
for index = 1:20
    x1 = x(index);
    x2 = 1-x1;
    Eq1 = L12 == 80.86/18.07*exp(-A12/(1.98*T));
    Eq2 = L21 == 18.07/80.86*exp(-A21/(1.98*T));
    Eq3 = log(G1) == -log(x1+L12*x2)+x2*(L12/(x1+L12*x2) - L21/(L21*x1+x2));
    Eq4 = log(G2) == -log(x2+L21*x1)-x1*(L12/(x1+L12*x2) - L21/(L21*x1+x2));
    Eq5 = log10(P1) == A1 - (B1/(T+C1-273));
    Eq6 = log10(P2) == A2 - (B2/(T+C2-273));
    Eq7 = 760 == P1*x1*G1 + P2*x2*G2;
    sol = vpasolve([Eq1, Eq2, Eq3, Eq4, Eq5, Eq6, Eq7], [L12, L21, G1, G2, P1, P2, T]);
    y(index) = double(sol.P1)*x1*double(sol.G1)/760;
    T_arr(index) = double(sol.T)-273;
    disp(x1);
    disp(double(sol.T)-273);
    disp(y(index));
end
subplot(1,2,1);
plot(x,T_arr);
hold on;
plot(y,T_arr);

subplot(1,2,2);
plot(x,y);