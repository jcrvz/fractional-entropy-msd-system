%% System 1: Simple mechanical oscillator
%% Mathematical formulae (ODE and FDE)
% Dynamic system under the traditional derivative definition is given by
% 
% $\frac{d^2x}{dt^2} + 2\zeta \omega_n \frac{dx}{dt} + \omega_n^2 x = 0$,
% 
% and under the fractional-conformable derivative definition by
% 
% $\frac{d^2x}{dt^2} + 2\zeta \omega_n^{\gamma} t^{\gamma-1} \frac{dx}{dt} 
% + \omega_n^{2\gamma} t^{2(\gamma-1)} x = 0$,
% 
% since the parameters for both models, $\zeta$ and $\omega$, are obtained 
% from the system's characteristics: 
 
% the body mass ($m~[\text{kg}]$),
m = 1; % kg

% the spring constant ($k~[\text{N/m}]$),
k = 1; % N/m

% the damping constant ($b~[\text{Ns/m}]$),
b = 1; % N.s/m 1 (underdamped), 2 (critically damped), 4 (overdamped)

% the initial conditions ($x_0~[\text{m}]$ and $v_0~[\text{m/s}]$),
x_0 = 1; % m
v_0 = 0.0; % m/s

% and the ambient temperature ($T_0~[\text{K}]$),
T_0 = 25 + 273.15; % K

% Then, the natural frequency ($\omega_n=\sqrt{k/m}$) and the damping coefficient 
% ($\zeta=b/2\sqrt{mk}$) are obtained with
omega_n = sqrt(k/m);
zeta    = b/(2*sqrt(k*m));

% Assign label for files to save
sz = sprintf('%.3f',zeta);
sz(sz=='.') = 'd';
if zeta < 1
    condLabel = sprintf('_underdamped_z%s',sz);
elseif z == 1
    condLabel = sprintf('_criticallydamped_z%s',sz);
else
    condLabel = sprintf('_overdamped_z%s',sz);
end

%% Varying the fractional-derivative order ($\gamma$)
% Time series for the simulations is defined such as,¡
tseries = linspace(1e-6,20,500);

% and the $\gamma$ values are given by,
gamma = linspace(0.65,1.55 ,5);
ng = numel(gamma);

% It was found the system's behaviour using the Ordinary Differential Equations 
% (ODEs), and the initial conditions prior described,
sys1.ode.func = @(t,x) mechsys_01(t,x,[m,b,k],'traditional');
[sys1.ode.t,sys1.ode.x] = ode45(sys1.ode.func, tseries, [x_0, v_0]);

% Plus, the entropy generation rate ($\dot{S}_{gen}~[\text{W/K}]$) is found 
% through,
sys1.ode.ds = 2*m/T_0*zeta*omega_n*sys1.ode.x(:,2).^2 ;
 
% This curve is displayed as follows,

Colors = [0 0 0;lines(ng)]; % Define colors for curves to be plotted

Graph1 = Graphics(['Sys1_displacement',condLabel]); 
ax1 = axes(Graph1.objID,'NextPlot','add','Box','on');

Graph2 = Graphics(['Sys1_velocity',condLabel]); 
ax2 = axes(Graph2.objID,'NextPlot','add','Box','on');

Graph3 = Graphics(['Sys1_entropy',condLabel]); 
ax3 = axes(Graph3.objID,'NextPlot','add','Box','on');

% Plot ODE's displacement in Graph1
plot(ax1,sys1.ode.t,sys1.ode.x(:,1),'LineWidth',1, ...
        'Color',Colors(1,:),'LineStyle','-','DisplayName','Trad.');
    
% Plot ODE's velocity in Graph2
plot(ax2,sys1.ode.t,sys1.ode.x(:,2),'LineWidth',1, ...
        'Color',Colors(1,:),'LineStyle','-','DisplayName','Trad.'); 
    
% Plot ODE's entropy in Graph3
plot(ax3,sys1.ode.t,sys1.ode.ds,'LineWidth',1, ...
        'Color',Colors(1,:),'LineStyle','-','DisplayName','Trad.');      
%% 
% Likewise, the system's behaviour using the Fractional Differential Equations 
% (FDEs), in the Conformable sense, was determined by varying $\gamma\in[0,1]$, 
% such as,

for iig = ng : -1 : 1
    g = gamma(iig); % Get the gamma value
    
    % Create the corresponding label
    sys1.fde(iig).label = sprintf('Conf., $$\\gamma = %.3f$$',g);
    
    % Read the fractional-conformable model for a given gamma value
    sys1.fde(iig).func = @(t,x) mechsys_01(t,x,[m,b,k,g],'conformable');
    
    % Solve this problem
    [sys1.fde(iig).t,sys1.fde(iig).x] = ode45(sys1.fde(iig).func, tseries, [x_0, v_0]);
    
    % Determine the entropy generation rate
    sys1.fde(iig).ds = 2*m/T_0*zeta*omega_n*sys1.fde(iig).x(:,2).^2 ;

    % Plot FDE's displacement in Graph1
    plot(ax1,sys1.fde(iig).t,sys1.fde(iig).x(:,1),'LineWidth',1, ...
        'Color',Colors(iig+1,:),'LineStyle','-','DisplayName',sys1.fde(iig).label);
    
    % Plot FDE's velocity in Graph2
    plot(ax2,sys1.fde(iig).t,sys1.fde(iig).x(:,2),'LineWidth',1, ...
        'Color',Colors(iig+1,:),'LineStyle','-','DisplayName',sys1.fde(iig).label);
    
    % Plot FDE's velocity in Graph3
    plot(ax3,sys1.fde(iig).t,sys1.fde(iig).ds,'LineWidth',1, ...
        'Color',Colors(iig+1,:),'LineStyle','-','DisplayName',sys1.fde(iig).label);
end

% Additional commands
strS = 'Entropy Gen., $$\dot{S}_{gen}(t)$$~[W/K]';
strV = 'Velocity, $$v(t)$$~[m/s]';
strX = 'Displacement, $$x(s)$$~[m]';
strT = 'Time, $$t$$~[s]';

leg1 = legend(ax1,'show'); 
set(leg1,'interpreter','latex','location','best','box','off','NumColumns',3); 
xlabel(ax1,strT,'Interpreter','LaTeX'); ylabel(ax1,strX,'Interpreter','LaTeX');
setall(Graph1,1,[3 1],12,1);

leg2 = legend(ax2,'show'); 
set(leg2,'interpreter','latex','location','best','box','off','NumColumns',3); 
xlabel(ax2,strT,'Interpreter','LaTeX'); ylabel(ax2,strV,'Interpreter','LaTeX');
setall(Graph2,1,[3 1],12,1);

leg3 = legend(ax3,'show'); 
set(leg3,'interpreter','latex','location','best','box','off','NumColumns',3); 
xlabel(ax3,strT,'Interpreter','LaTeX'); ylabel(ax3,strS,'Interpreter','LaTeX');
setall(Graph3,1,[3 1],12,1);

%% Theoretical expression for the entropy generation rate and its comparison 
%  against the numerical result
% The auxiliar parameters like $\omega_d$, $A$ and $\phi$, are determined with
% 
% $\omega_d = \omega_n \sqrt{1-\zeta^2}$, $A = \sqrt{ \left( \frac{x_0\zeta}{\sqrt{1-\zeta^2}} 
% + \frac{v_0}{\omega_d} \right)^2 +x_0^2 }$, and $\phi = \tan^{-1}\left(\frac{\omega_d 
% x_0 }{v_0 + x_0 \zeta \omega_n} \right)$.
%%
omega_d = omega_n*sqrt(1 - zeta^2);
A = sqrt(( x_0*zeta/sqrt(1-zeta^2) + v_0/omega_d )^2 + x_0^2);
phi = atan(omega_d*x_0/(v_0+x_0*zeta*omega_n));
%% 
% Thence, the entropy generation rate ($\dot{S}_{gen}~[\text{W/K}]$) can 
% be found by using, 
% 
%    

% func_entropy = @(t) (-(m/T_0)*A^2*exp(-2*zeta*omega_n*t)).*...
%     ((omega_d*cos(omega_d*t + phi) - zeta*omega_n*sin(omega_d*t + phi)).*...
%     (omega_n^2*sin(omega_d*t + phi) - zeta*omega_n*(omega_d*cos(omega_d*t + phi) - ...
%     zeta*omega_n*sin(omega_d*t + phi)) - omega_d.*(omega_d*sin(omega_d*t + phi) + ...
%     zeta*omega_n*cos(omega_d*t + phi))));
%  format long g
%% 
% Now, we can qualitatively compare the numerical result against the theoretical 
% one as shown,

% figure, plot(tseries,func_entropy(tseries),'k','linewidth',1), hold on,
% plot(sys1.ode.t,sys1.ode.ds,'r--','linewidth',1); hold off, ax0 = gca;
% leg0 = legend('Theoretical','Numerical'); set(leg0,'interpreter','latex',...
%     'location','best','box','off','fontsize',12); 
% xlabel(ax0,strT,'Interpreter','LaTeX','fontsize',12); ylabel(ax0,strS,...
%     'Interpreter','LaTeX','fontsize',12);
% set(ax0,'ticklabelinterpreter','latex','fontsize',12,'linewidth',1);
%% 
% And, the quantitative comparison between the theoretical and numerical 
% form is obtained with: 
% 
% $$\Delta S_{gen} = \int_{0}^{\infty} \dot{S}_{gen}(\tau) d\tau$$
% 
% Thus,

totalEntropy_theoretical = integral(func_entropy,tseries(1),tseries(end));
%disp(totalEntropy_theoretical);
%disp(['$\Delta S_{gen} = ',num2str(totalEntropy_theoretical)*1e3,'\times 10^{-3}$~W/K'])
%% 
% and, 

totalEntropy_numerical = trapz(sys1.ode.ds)*diff(sys1.ode.t(1:2));
%disp(totalEntropy_theoretical);
%% 
% Finally, it is possible to determine the total entropy for all the previously 
% discussed approaches

for iig = 1 : ng    
    % Determine the total entropy from entropy generation rate results
    sys1.fde(iig).totalEntropy = trapz(sys1.fde(iig).ds)*diff(sys1.fde(iig).t(1:2));
end

% Store data in a table
TotalEntropy = table({'Traditional',sys1.fde.label}',...
   [totalEntropy_numerical,[sys1.fde.totalEntropy]]');
TotalEntropy.Properties.VariableNames = {'Model','Total_Entropy'};

%disp(TotalEntropy)
%% 
% Plot the Total Entropy variation
strTS = 'Total Entropy Gen., $$\Delta{S}_{gen}(t)$$~[J/K]';
strGamma = 'Fractional Order, $$\gamma$$';

Graph4 = Graphics('Sys1_totalentropy'); 
ax4 = axes(Graph4.objID,'NextPlot','add','Box','on');

gammas = gamma(end:-1:1);
totalEntropies = [sys1.fde(end:-1:1).totalEntropy];
[~,idx] = min(totalEntropies);

plot(gammas,totalEntropies,'k-','LineWidth',1,'DisplayName','Conformable');
plot(1,totalEntropy_numerical,'ob','LineWidth',1,'DisplayName','Traditional');
plot(gammas(idx),totalEntropies(idx),'*r','LineWidth',1,'DisplayName','Minimum');

leg4 = legend(ax4,'show'); 
set(leg4,'interpreter','latex','location','best','box','off','NumColumns',1); 
xlabel(ax4,strGamma,'Interpreter','LaTeX'); ylabel(ax4,strTS,'Interpreter','LaTeX');
setall(Graph4,1,[2 1],12,1);

%%
Graph5 = Graphics('Sys1_totalentropy'); 
ax5 = axes(Graph5.objID,'NextPlot','add','Box','on');

idx_ = numel(gamma) - idx + 1;

ss = [sys1.fde.ds]';
[Ti,Ga] = meshgrid(tseries,gamma);
surf(Ti,Ga,ss,'FaceAlpha',0.95), shading interp, hold on,
plot3(Ti(idx_,:),Ga(idx_,:),ss(idx_,:),'r','linewidth',1)
set(ax5,'ZScale','line'), view([35.5 29.2])

xlabel(ax5,strT,'Interpreter','LaTeX'); 
ylabel(ax5,strGamma,'Interpreter','LaTeX');
zlabel(ax5,strS,'Interpreter','LaTeX');

% set(gca,'ColorScale','log')
%ax5.Children(2).FaceAlpha = 0.95;
% zlim([1e-6 1e0])
setall(Graph5,1,[2 1],12,1);