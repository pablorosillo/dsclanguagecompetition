% Dynamical Systems and Chaos
% Final Project on Language Competition
% Pablo Rosillo

% This code fits SIES and Miralles et al. data to Louf et al. model

% Parameters: theta = [r;s;q]

set(groot,'defaulttextinterpreter','latex'); 
set(groot,'defaultAxesTickLabelinterpreter','latex'); 
set(groot,'defaultlegendinterpreter','latex');


t = [1989; 1991; 1995; 2005; 2010];
t = t - 1989;
c1 = [37.2; 34.6; 36.8; 25.9; 23]/100; %Valencian speakers
errorc1 = [1.1;0.9;1.9;0.9;0.8]/100;
c2 = [48; 44; 44.7; 57.3; 60.4]/100;
errorc2 = [1.2; 1; 1.9; 1.2; 1.2]/100;
c = [c1 c2];
errorc = [errorc1 errorc2];


theta0=rand(3,1);

[theta,Rsdnrm,Rsd,ExFlg,OptmInfo,Lmda,Jmat]=lsqcurvefit(@kinetics,theta0,t,c, [0.001;0.001;0.001], [1;1;1]);

disp(Rsdnrm)

fprintf(1,'\tRate Constants:\n')
for k1 = 1:length(theta)
    fprintf(1, '\t\tTheta(%d) = %8.5f\n', k1, theta(k1))
end

tchange = 2030;

tv = linspace(min(t), tchange-1989);
Cfit = kinetics(theta, tv);

figure
errorbar(t+1989, c(:,1), errorc(:,1), '.')
hold on
errorbar(t+1989, c(:,2), errorc(:,2), 'x')
hold on
errorbar(t+1989, 1-c(:,1)-c(:,2), sqrt(errorc(:,2).^2 + errorc(:,1).^2), 'o')
hold on
hlp = plot(tv+1989, Cfit);
hold on
plot(tv+1989, 1-Cfit(:,1)-Cfit(:,2))
hold off
grid
xlabel('Year', 'FontSize', 12)
legend("Valenci\`a speakers", 'Spanish speakers', 'Bilingual speakers',...
    '$p_\mathrm{A}$', '$p_\mathrm{B}$', '$p_\mathrm{AB}$')
ylim([0,1])
xlim([1989, tchange])
ytickformat('%.5f')
xtickformat('%.0f')

function C=kinetics(theta,t, c0)
c0=[0.372;0.48];
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






