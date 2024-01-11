
% Gibbs free energy minimization algorithm - CO2 methanation

close all, clear all, clc

% chemical properties and in/out molar balance
data = readtable('summary.xlsx','Sheet','properties');

% operational variables and initial guess
Peval = 1; % bar
Teval = 100:20:900; % ºC
ratio = 4; % 4 mol H2 / mol CO2
initial_guess = rand([1,8]); % for Teval(1) / CH4 CO CO2 H2 H2O C_graphitic C_nanotubes C_amorphous

CH4 = NaN(1,length(Teval));
CO = NaN(1,length(Teval));
CO2 = NaN(1,length(Teval));
H2 = NaN(1,length(Teval));
H2O = NaN(1,length(Teval));
C_T = NaN(1,length(Teval));

for i = 1:1:length(Teval)
    
    [x(i,:),Gt] = min_Gt(Peval,Teval(1,i)+273.15,ratio,data,initial_guess);
    
    CH4(1,i) = x(i,1);
    CO(1,i) = x(i,2);
    CO2(1,i) = x(i,3);
    H2(1,i) = x(i,4);
    H2O(1,i) = x(i,5);
    C_T(1,i) = x(i,end) + x(i,end-1) + x(i,end-2);
    
    X(1,i) = (1 - x(i,3)) * 100;
    
    initial_guess = x(i,:);
    
end

set(0,'DefaultAxesFontName', 'Times New Roman')

f = figure();

tiledlayout(2,1,"TileSpacing","compact")

nexttile
plot(Teval,CH4,'o-',Teval,CO,'d-',Teval,CO2,'^-',Teval,C_T,'s-',Teval,H2,'>-',Teval,H2O,'<-')
ylabel('Yield (mol/mol CO_2 in)')
xlim([Teval(1) Teval(end)])
ylim([-0.01 3.5])
grid on
legend('CH_4','CO','CO_2','C','H_2','H_2O')
legend('box','off')
legend('Location','eastoutside')
set(gca, 'FontSize', 13)

nexttile
plot(Teval,X,'o-')
ylabel('Conversion (%)')
xlabel('Temperature (°C)')
grid on
set(gca, 'FontSize', 13)

Reval = 0:5:90; % water removal
[T,R] = meshgrid(Teval,Reval);

for i = 1:1:length(Teval)
    for j = 1:1:length(Reval)
    [x_R,Gt] = min_Gt_insitu(Peval,Teval(1,i)+273.15,data,x(i,:),Reval(j)/100);

    CH4(j,i) = x_R(1);
    CO(j,i) = x_R(2);
    CO2(j,i) = x_R(3);
    H2(j,i) = x_R(4);
    H2O(j,i) = x_R(5);
    C_T(j,i) = x_R(end) + x_R(end-1) + x_R(end-2);
    
    X(j,i) = (1 - x_R(3)) * 100;
    end
end

g = figure;

tiledlayout(3,2,"TileSpacing","compact")

nexttile
surf(T,R,X)
zlabel('CO_2 conv. (%)')
xlabel('T (°C)')
ylabel('Removal (%)')
xlim([Teval(1) Teval(end)])
ylim([Reval(1) Reval(end)])
title('(a)')

nexttile
surf(T,R,CH4)
zlabel('CH_4/CO_2')
xlabel('T (°C)')
ylabel('Removal (%)')
xlim([Teval(1) Teval(end)])
ylim([Reval(1) Reval(end)])
title('(b)')

nexttile
surf(T,R,C_T)
zlabel('C/CO_2')
xlabel('T (°C)')
ylabel('Removal (%)')
xlim([Teval(1) Teval(end)])
ylim([Reval(1) Reval(end)])
title('(c)')

nexttile
surf(T,R,H2)
zlabel('H_2/CO_2')
xlabel('T (°C)')
ylabel('Removal (%)')
xlim([Teval(1) Teval(end)])
ylim([Reval(1) Reval(end)])
title('(d)')

nexttile
surf(T,R,CO)
zlabel('CO/CO_2')
xlabel('T (°C)')
ylabel('Removal (%)')
xlim([Teval(1) Teval(end)])
ylim([Reval(1) Reval(end)])
title('(e)')

nexttile
surf(T,R,H2O)
zlabel('H_2O/CO_2')
xlabel('T (°C)')
ylabel('Removal (%)')
xlim([Teval(1) Teval(end)])
ylim([Reval(1) Reval(end)])
title('(f)')

% main
function [x,Gt] = min_Gt(Peval,T,ratio,data,initial_guess)

% standard properties
P0 = 1; % bar
T0 = 298; % K

% fugacity coefficient
Tr = T./data.Tc_K;
Pr = Peval./data.Pc_bar;
B0 = 0.083 - 0.422./(Tr.^1.6);
B1 = 0.139 - 0.172./(Tr.^4.2);
phi = exp((Pr./Tr) .* (B0 + data.w.*B1));

P = phi.*(Peval/P0); % bar/bar

in = readmatrix('summary.xlsx','Sheet','in');

out = readmatrix('summary.xlsx','Sheet','out'); out = out(:,2:end)';

in(4,:) = ratio * in(4,:);

in = sum(in); in = in(:,2:end)';

u_0 = (T / T0) * data.Gf ...
     + (1 - (T / T0)) * data.Hf ...
     - data.CpA * (T * log(T / T0) - T + T0) ...
     - (data.CpB/2) * (T^2 - 2 * T * (T0^1) + 1 * (T0^1)) ...
     - (data.CpC/6) * ( T^3 - 3 * T * (T0^2) + 2 * (T0^3)) ...
     - (data.CpD/12) * ( T^4 - 4 * T * (T0^3) + 3 * (T0^4)) ...
     - (data.CpE/20) * ( T^5 - 5 * T * (T0^4) + 4 * (T0^5));

NN = length(u_0);

u_0(NN+1,1) = u_0(end) + 1.971e+03 - 2.801 * T; % mwcnt
u_0(NN+2,1) = - (1.392e-12) * (T^4) + (4.725e-9) * (T^3) - (6.602e-6) * (T^2) ...
    + (2.352e-3) * T + 3.560; % amorphous

out(:,NN+1) = out(:,NN);
out(:,NN+2) = out(:,NN);

LB = zeros(1,length(u_0)); % ni > 0
options = optimset('Algorithm','sqp');
[x,Gt] = fmincon(@phase_equilibrium,initial_guess,[],[],out,in,LB,[],[],options,P,T,u_0');

end

function [x_R,Gt] = min_Gt_insitu(Peval,T,data,initial_composition,R)

% standard properties
P0 = 1; % bar
T0 = 298; % K

% fugacity coefficient
Tr = T./data.Tc_K;
Pr = Peval./data.Pc_bar;
B0 = 0.083 - 0.422./(Tr.^1.6);
B1 = 0.139 - 0.172./(Tr.^4.2);
phi = exp((Pr./Tr) .* (B0 + data.w.*B1));

P = phi.*(Peval/P0); % bar/bar

in = readmatrix('summary.xlsx','Sheet','in_insitu');

out = readmatrix('summary.xlsx','Sheet','out'); out = out(:,2:end)';

initial_composition(5) = initial_composition(5) * (1 - R);

for jj = 1:8
    in(jj,:) = in(jj,:) * initial_composition(jj);
end

in = sum(in); in = in(:,2:end)';

u_0 = (T / T0) * data.Gf ...
     + (1 - (T / T0)) * data.Hf ...
     - data.CpA * (T * log(T / T0) - T + T0) ...
     - (data.CpB/2) * (T^2 - 2 * T * (T0^1) + 1 * (T0^1)) ...
     - (data.CpC/6) * ( T^3 - 3 * T * (T0^2) + 2 * (T0^3)) ...
     - (data.CpD/12) * ( T^4 - 4 * T * (T0^3) + 3 * (T0^4)) ...
     - (data.CpE/20) * ( T^5 - 5 * T * (T0^4) + 4 * (T0^5));

NN = length(u_0);

u_0(NN+1,1) = u_0(end) + 1.971e+03 - 2.801 * T; % mwcnt
u_0(NN+2,1) = - (1.392e-12) * (T^4) + (4.725e-9) * (T^3) - (6.602e-6) * (T^2) ...
    + (2.352e-3) * T + 3.560; % amorphous

out(:,NN+1) = out(:,NN);
out(:,NN+2) = out(:,NN);

LB = zeros(1,length(u_0)); % ni > 0
options = optimset('Algorithm','sqp');
[x_R,Gt] = fmincon(@phase_equilibrium,initial_composition,[],[],out,in,LB,[],[],options,P,T,u_0');

end

% Objective function
function G = phase_equilibrium(nj,P,T,u_0)

R = 1.987; % cal/mol/K

Enj = sum(nj(1:end-1));

G = 0;
for i = 1:length(nj)-3
    G = G + nj(i) * ( u_0(i) + R * T * ( log(P(i) * nj(i) / Enj) ) );
end

for j = length(nj)-2:length(nj)
    G = G + nj(j) * u_0(j);
end

end
