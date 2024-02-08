%% System 1: Simple mechanical oscillator
% the body mass ($m~[\text{kg}]$),
m = 1; % kg

% the spring constant ($k~[\text{N/m}]$),
k = 1; % N/m

% the damping constant ($b~[\text{Ns/m}]$),
%b = 0; % N.s/m (undamped)
b = 1/4; % N.s/m (underdamped)
%b = 2; % N.s/m (critically damped)
% b = 4; % N.s/m (overdamped)

% the initial conditions ($x_0~[\text{m}]$ and $v_0~[\text{m/s}]$),
x_0 = 1; % m
v_0 = 0.0; % m/s

% and the ambient temperature ($T_0~[\text{K}]$),
T_0 = 25 + 273.15; % K

% Then, the natural frequency ($\omega_n=\sqrt{k/m}$) and the damping coefficient
% ($\zeta=b/2\sqrt{mk}$) are obtained with
omega_n = sqrt(k/m);
zeta    = b/(2*sqrt(k*m));

% Time series for the simulations is defined such as,
tseries = linspace(1e-6,20,500);

% and the $\gamma$ values are given by,
gamma = sort(unique([1,linspace(0.65,1.55,25)])); % 5, 25

% For plotting
fontsize  = 16;
aspect_2D = [3 2];
aspect_3D = [1 1];

% Assign label for files to save
sz = sprintf('%.3f',zeta);
sz(sz=='.') = 'd';
if zeta == 0
    condLabel = sprintf('_undamped_z%s',sz);
elseif zeta < 1
    condLabel = sprintf('_underdamped_z%s',sz);
elseif zeta == 1
    condLabel = sprintf('_criticallydamped_z%s',sz);
else
    condLabel = sprintf('_overdamped_z%s',sz);
end

%% Varying the fractional-derivative order ($\gamma$)

% Important code's features
ng = numel(gamma);
iTrad = find(gamma == 1); % this is the index for traditional derivative

% This curve is displayed as follows,
Colors = zeros(ng,3);
Colors(gamma ~= 1,:) = lines(ng - ~isempty(iTrad)); % Define colors

Graph1 = Graphics(['Sys1_displacement',condLabel]);
ax1 = axes(Graph1.objID,'NextPlot','add','Box','on');

Graph2 = Graphics(['Sys1_velocity',condLabel]);
ax2 = axes(Graph2.objID,'NextPlot','add','Box','on');

Graph3 = Graphics(['Sys1_entropy',condLabel]);
ax3 = axes(Graph3.objID,'NextPlot','add','Box','on');

% Plot curves from the FDE model
sys1 = struct();
for iig = 1 : ng
    g = gamma(iig); % Get the gamma value
    
    % Create the corresponding label
    if iig == iTrad
        sys1(iig).label = 'Traditional';
    else
        sys1(iig).label = sprintf('$$\\gamma = %.3f$$',g);
    end
    
    % Read the fractional-conformable model for a given gamma value
    sys1(iig).func = @(t,x) mechsys_01(t,x,[m,b,k,g]);
    
    % Solve this problem
    [sys1(iig).t,x_dummy] = ode45(sys1(iig).func, tseries, [x_0, v_0]);
    
    % Store the displacement and velocity
    sys1(iig).x = x_dummy(:,1);
    sys1(iig).v = x_dummy(:,2);
    
    % Determine the entropy generation rate
    sys1(iig).ds = 2*m/T_0*zeta*omega_n*sys1(iig).v.^2 ;
    
    % Plot FDE's displacement in Graph1
    plot(ax1,sys1(iig).t,sys1(iig).x,'LineWidth',1, ...
        'Color',Colors(iig,:),'LineStyle','-','DisplayName',sys1(iig).label);
    
    % Plot FDE's velocity in Graph2
    plot(ax2,sys1(iig).t,sys1(iig).v,'LineWidth',1, ...
        'Color',Colors(iig,:),'LineStyle','-','DisplayName',sys1(iig).label);
    
    % Plot FDE's velocity in Graph3
    plot(ax3,sys1(iig).t,sys1(iig).ds,'LineWidth',1, ...
        'Color',Colors(iig,:),'LineStyle','-','DisplayName',sys1(iig).label);
end

% Additional commands
strS = 'Entropy Gen., $$\dot{S}_{gen}(t)$$~[W/K]';
strV = 'Velocity, $$v(t)$$~[m/s]';
strX = 'Displacement, $$x(s)$$~[m]';
strT = 'Time, $$t$$~[s]';

leg1 = legend(ax1,'show');
set(leg1,'interpreter','latex','location','best','box','off','NumColumns',3);
xlabel(ax1,strT,'Interpreter','LaTeX'); ylabel(ax1,strX,'Interpreter','LaTeX');
setall(Graph1,1,aspect_2D,fontsize,2);

leg2 = legend(ax2,'show');
set(leg2,'interpreter','latex','location','best','box','off','NumColumns',3);
xlabel(ax2,strT,'Interpreter','LaTeX'); ylabel(ax2,strV,'Interpreter','LaTeX');
setall(Graph2,1,aspect_2D,fontsize,2);

leg3 = legend(ax3,'show');
set(leg3,'interpreter','latex','location','best','box','off','NumColumns',3);
xlabel(ax3,strT,'Interpreter','LaTeX'); ylabel(ax3,strS,'Interpreter','LaTeX');
setall(Graph3,1,aspect_2D,fontsize,2);

%% Total entropy for all found results
% Finally, it is possible to determine the total entropy for all the previously
% discussed approaches

for iig = 1 : ng
    % Determine the total entropy from entropy generation rate results
    sys1(iig).totalEntropy = trapz(sys1(iig).ds)*diff(sys1(iig).t(1:2));
end

% Store data in a table
TotalEntropy = table({sys1.label}',...
    [sys1.totalEntropy]');
TotalEntropy.Properties.VariableNames = {'Model','Total_Entropy'};

% Identify the minimal total entropy
totalEntropies = [sys1(:).totalEntropy];
[~,idmin] = min(totalEntropies);

% Plot the Total Entropy variation
strTS = 'Total Entropy Gen., $$\Delta{S}_{gen}(t)$$~[J/K]';
strGamma = 'Fractional Order, $$\gamma$$';

Graph4 = Graphics(['Sys1_totalentropy',condLabel]);
ax4 = axes(Graph4.objID,'NextPlot','add','Box','on');

plot(gamma,totalEntropies,'k-','LineWidth',1,'DisplayName','Conformable');
plot(gamma(idmin),totalEntropies(idmin),'*r','LineWidth',1,'DisplayName','Minimal');

leg4 = legend(ax4,'show');
set(leg4,'interpreter','latex','location','best','box','off','NumColumns',1);
xlabel(ax4,strGamma,'Interpreter','LaTeX'); ylabel(ax4,strTS,'Interpreter','LaTeX');
setall(Graph4,1,[3 1],12,1);

%%
Graph5 = Graphics(['Sys1_totalentropy_surf',condLabel]);
ax5 = axes(Graph5.objID,'NextPlot','add','Box','on');

TEntr = [sys1.ds]';
[Ti,Ga] = meshgrid(tseries,gamma);
surf(Ti,Ga,TEntr,'FaceAlpha',0.95), shading interp, hold on,
plot3(Ti(idmin,:),Ga(idmin,:),TEntr(idmin,:),'r','linewidth',1)
set(ax5,'ZScale','line'), view([35.5 29.2])

xlabel(ax5,strT,'Interpreter','LaTeX');
ylabel(ax5,strGamma,'Interpreter','LaTeX');
zlabel(ax5,strS,'Interpreter','LaTeX');
setall(Graph5,1,[2 1],12,1);
%print(Graph5.objID,Graph5.fileName,'-r360','-djpeg','-noui');

%% Plot x-v in 3D
Graph6 = Graphics(['Sys1_characteristic',condLabel]);
ax6 = axes(Graph6.objID,'NextPlot','add','Box','on');

X = [sys1(:).x];
V = [sys1(:).v];
plot3(Ga',X,V,'linewidth',1), 
set(ax6,'ZScale','line','XTick',gamma,'XGrid','on'), view([-32 28])

xlabel(ax6,strGamma,'Interpreter','LaTeX');
ylabel(ax6,strX,'Interpreter','LaTeX','Position',[-0.5 -0.18 -0.4]);
zlabel(ax6,strV,'Interpreter','LaTeX');
setall(Graph6,2,[1 1],12,1);

Graph7 = Graphics(['Sys1_characteristic_surf',condLabel]);
ax7 = axes(Graph7.objID,'NextPlot','add','Box','on');


X = [sys1(:).x];
V = [sys1(:).v];
surf(Ga',X,V,'FaceAlpha',0.9), shading interp
view([-44.8 28]), light

xlabel(ax7,strGamma,'Interpreter','LaTeX');
ylabel(ax7,strX,'Interpreter','LaTeX','Position',[-0.5 -0.5]);
zlabel(ax7,strV,'Interpreter','LaTeX');
setall(Graph7,2,[1 1],12,1);
%print(Graph7.objID,Graph7.fileName,'-r360','-djpeg','-noui');

%%