%% System 1: Simple mechanical oscillator
% the body mass ($m~[\text{kg}]$),
m = 1; % kg

% the spring constant ($k~[\text{N/m}]$),
k = 1; % N/m

% the damping constant ($b~[\text{Ns/m}]$),
%b = 0; % N.s/m (undamped)
%b = 1/4; % N.s/m (underdamped)
%b = 2; % N.s/m (critically damped)
b = 4; % N.s/m (overdamped)

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
tseries = linspace(1e-6,10,1000);

% and the $\gamma$ values are given by,
gamma = sort(unique([1,linspace(.65,1,6)])); % (5 | 25)

% For plotting
fontsize  = 14;
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
        sys1(iig).label = sprintf('$$\\gamma = %.2f$$',g);
    end
    
    % Read the fractional-conformable model for a given gamma value
    sys1(iig).func = @(t,x) mechsys_01(t,x,[m,b,k,g]);
    
    % Solve this problem
    [sys1(iig).t,x_dummy] = ode45(sys1(iig).func, tseries, [x_0, v_0]);
    
    % Store the displacement and velocity
    sys1(iig).x = x_dummy(:,1);
    sys1(iig).v = x_dummy(:,2);
    
    % Print a table
    sinf = stepinfo(sys1(iig).x,tseries);
    fprintf('z=%5.3f\tg=%5.3f\tt_s=%5.4f\tT_s=%5.4f\n',zeta,g,...
        4/zeta/omega_n,sinf.SettlingTime);
    
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
strS = '$$\dot{S}_{gen}(t)$$~[W/K]';
strV = '$$v(t)$$~[m/s]';
strX = '$$x(t)$$~[m]';
strT = '$$t$$~[s]';
strTS = '$$\Delta{S}_{gen}(t)$$~[J/K]';
strGamma = '$$\gamma$$';

leg1 = legend(ax1,'show');
set(leg1,'interpreter','latex','location','best','box','off','NumColumns',2);
xlabel(ax1,strT,'Interpreter','LaTeX'); ylabel(ax1,strX,'Interpreter','LaTeX');
ax1.XLim = [0 10];
setall(Graph1,2,aspect_2D,fontsize,1);

leg2 = legend(ax2,'show');
set(leg2,'interpreter','latex','location','best','box','off','NumColumns',2);
xlabel(ax2,strT,'Interpreter','LaTeX'); ylabel(ax2,strV,'Interpreter','LaTeX');
ax2.XLim = [0 10];
setall(Graph2,2,aspect_2D,fontsize,1);

leg3 = legend(ax3,'show');
set(leg3,'interpreter','latex','location','best','box','off','NumColumns',2);
xlabel(ax3,strT,'Interpreter','LaTeX'); ylabel(ax3,strS,'Interpreter','LaTeX');
ax3.XLim = [0 10];
setall(Graph3,2,aspect_2D,fontsize,1);


%% Plot x-v in 3D
Graph6 = Graphics(['Sys1_characteristic',condLabel]);
ax6 = axes(Graph6.objID,'NextPlot','add','Box','on');

TEntr = [sys1.ds]';
[Ti,Ga] = meshgrid(tseries,gamma);

X = [sys1(:).x];
V = [sys1(:).v];
plot3(Ga',X,V,'linewidth',1), 
set(ax6,'ZScale','line','XTick',gamma,'XGrid','on'), view([-25 25])
%ylim([-2 2]), zlim([-2 2])

xlabel(ax6,strGamma,'Interpreter','LaTeX');
ylabel(ax6,strX,'Interpreter','LaTeX');%,'Position',[0.57701 -0.084311 -2.3722]);
zlabel(ax6,strV,'Interpreter','LaTeX');
setall(Graph6,2,aspect_3D,fontsize,1);


%%