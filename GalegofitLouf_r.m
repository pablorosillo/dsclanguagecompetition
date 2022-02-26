
% Dynamical Systems and Chaos
% Final Project on Language Competition
% Pablo Rosillo

% This code fits IGE data to Louf et al. model

% Parameters: theta = [r;s;q]

set(groot,'defaulttextinterpreter','latex'); 
set(groot,'defaultAxesTickLabelinterpreter','latex'); 
set(groot,'defaultlegendinterpreter','latex');


t = [1992;2003;2008;2013;2018];
tf = t(1);
t = t - tf;
c1 = [38;42.9;29.9;30.8;30.3]/100; %Galician speakers
c2 = [10.6;19.6;20;25.9;24.2]/100;
c = [c1 c2];


theta0=rand(3,1);

[theta,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat]=lsqcurvefit(@kinetics,theta0,t,c, [0.001;0.001;0.001], [1;1;1]);

disp(Rsdnrm)

fprintf(1,'\tParameters:\n')
for k1 = 1:length(theta)
    fprintf(1, '\t\tTheta(%d) = %8.5f\n', k1, theta(k1))
end

tv = linspace(min(t), 2030-1989);
Cfit = kinetics(theta, tv);

figure
plot(t+1989, c(:,1), '.')
hold on
plot(t+1989, c(:,2), 'x')
hold on
plot(t+1989, 1-c(:,1)-c(:,2), 'o')
hold on
hlp = plot(tv+1989, Cfit);
hold on
plot(tv+1989, 1-Cfit(:,1)-Cfit(:,2))
hold off
grid
xlabel('Year', 'FontSize', 12)
legend('Galego speakers', 'Spanish speakers', 'Bilingual speakers',...
    '$p_\mathrm{A}$', '$p_\mathrm{B}$', '$p_\mathrm{AB}$')
ylim([0,1])
xlim([1989, 2030])
ytickformat('%.1f')
xtickformat('%.0f')

function C=kinetics(theta,t)
c0=[0.38;0.106];
[T,Cv]=ode45(@DifEq,t,c0);
%
    function dC=DifEq(t,c)
    dcdt=zeros(2,1);
    dcdt(1) = theta(1)*(1-c(1)-c(2))*theta(2)*(c(1)+theta(3)*(1-c(1)-c(2)))-c(1)*(1-theta(2))*(c(2)+(1-theta(3))*(1-c(1)-c(2)));
    dcdt(2) = theta(1)*(1-c(1)-c(2))*(1-theta(2))*(c(2)+(1-theta(2))*(1-c(1)-c(2))-c(2)*theta(2)*(c(1)+theta(3)*(1-c(1)-c(2))));
    dC=dcdt;
    end
C=Cv;
end
