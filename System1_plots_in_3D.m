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
tseries = linspace(1e-6,60,500);

% and the $\gamma$ values are given by,
gamma = sort(unique([1,linspace(0.65,1,31)])); % (5 | 30)

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
    
end

% Additional commands
strS = '$$\dot{S}_{gen}(t)$$~[W/K]';
strV = '$$v(t)$$~[m/s]';
strX = '$$x(t)$$~[m]';
strT = '$$t$$~[s]';
strTS = '$$\Delta{S}_{gen}(t)$$~[J/K]';
strGamma = '$$\gamma$$';

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

%%
Graph5 = Graphics(['Sys1_totalentropy_surf',condLabel]);
ax5 = axes(Graph5.objID,'NextPlot','add','Box','on');

TEntr = [sys1.ds]';
[Ti,Ga] = meshgrid(tseries,gamma);
surf(Ti,Ga,TEntr,'FaceAlpha',0.95), colormap jet, shading interp, hold on,
plot3(Ti(idmin,:),Ga(idmin,:),TEntr(idmin,:),'r','linewidth',1)
set(ax5,'ZScale','line'), view([35.5 29.2])
ylim([min(gamma),max(gamma)]), 
xlim([0 10])
light('Position',[9,1,3e-3],'Style','local')

xlabel(ax5,strT,'Interpreter','LaTeX');
ylabel(ax5,strGamma,'Interpreter','LaTeX');
zlabel(ax5,strS,'Interpreter','LaTeX');
setall(Graph5,2,aspect_2D,fontsize,1);
%print(Graph5.objID,Graph5.fileName,'-r360','-djpeg','-noui');

%% Plot x-v in 3D

Graph7 = Graphics(['Sys1_characteristic_surf',condLabel]);
ax7 = axes(Graph7.objID,'NextPlot','add','Box','on');


X = [sys1(:).x];
V = [sys1(:).v];
% --------------------------------------- (m|g|c|y)
surf(Ga',X,V,'FaceAlpha',0.7,'FaceColor','m','EdgeColor','none','EdgeAlpha',0.3),% shading interp,
view([-44.8 28]), light('Position',[0.6,-2,0.8],'Style','local')

xlabel(ax7,strGamma,'Interpreter','LaTeX');
ylabel(ax7,strX,'Interpreter','LaTeX');%,'Position',[0.57701 -0.084311 -2.3722]);
zlabel(ax7,strV,'Interpreter','LaTeX');
setall(Graph7,2,aspect_3D,fontsize,1);
%print(Graph7.objID,Graph7.fileName,'-r360','-djpeg','-noui');

%%