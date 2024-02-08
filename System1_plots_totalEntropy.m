%% System 1: Simple mechanical oscillator
% the body mass ($m~[\text{kg}]$),
m = 1; % kg

% the spring constant ($k~[\text{N/m}]$),
k = 1; % N/m

% the damping constant ($b~[\text{Ns/m}]$),
%b = 0; % N.s/m (undamped)
%b = 1/4; % N.s/m (underdamped)
%b = 2; % N.s/m (critically damped)
% b = 4; % N.s/m (overdamped)
bseries = sort(unique([0,1/4,2,4,linspace(0,4,100),1./logspace(0,3,10)])); %linspace(0,4,30)
%bseries = [0,1/4,2,4];

% the initial conditions ($x_0~[\text{m}]$ and $v_0~[\text{m/s}]$),
x_0 = 1; % m
v_0 = 0.0; % m/s

% and the ambient temperature ($T_0~[\text{K}]$),
T_0 = 25 + 273.15; % K

% Then, the natural frequency ($\omega_n=\sqrt{k/m}$) and the damping coefficient
% ($\zeta=b/2\sqrt{mk}$) are obtained with
omega_n = sqrt(k/m);
zetaseries    = bseries/(2*sqrt(k*m));

% Time series for the simulations is defined such as,
tseries = linspace(1e-6,500,1000);

% and the $\gamma$ values are given by,
gamma = sort(unique([1,linspace(1e-3,1,91)])); % (6 | 31)

% For plotting
fontsize  = 14;
aspect_2D = [3 2];
aspect_3D = [3 2];

% Assign label for files to save
% sz = sprintf('%.3f',zeta);
% sz(sz=='.') = 'd';


%% Varying the fractional-derivative order ($\gamma$)

% Important code's features
ng = numel(gamma);
iTrad = find(gamma == 1); % this is the index for traditional derivative

[Ti,Ga] = meshgrid(tseries,gamma);

% This curve is displayed as follows
Colors = [
    0,0,0
    55,126,184; % azul
    228,26,28; % rojo
    77,175,74; % verde
    231,41,138; % fucsia
    152,78,163; % morado
    255,127,0; % naranja
    ]/255;
% Colors = lines(numel(zetaseries));
% Colors(gamma ~= 1,:) = (4 - ~isempty(iTrad)); % Define colors

% Additional commands
strS = '$$\dot{S}_{gen}(t)$$~[W/K]';
strV = '$$v(t)$$~[m/s]';
strX = '$$x(s)$$~[m]';
strT = '$$t$$~[s]';
strZ = '$$\zeta$$';
strTS = '$$\Delta{S}_{gen}$$~[J/K]';
strGamma = '$$\gamma$$';

% Plot curves from the FDE model
sys1 = struct();
totalEntropies = nan(size(gamma));
minTotalEntropies = nan(size(bseries));
idmin = nan(size(bseries));
qq = nan(numel(gamma),numel(bseries));
settlingtimes = qq;
peaks = qq;
zetaCell = cell(1,numel(bseries));
for iib = 1 : numel(bseries)
    % Read data from b and zeta series
    b = bseries(iib);
    zeta = zetaseries(iib);
 
    for iig = 1 : ng
        g = gamma(iig); % Get the gamma value
        
        % Create the corresponding label
        if iig == iTrad
            sys1(iig,iib).label = 'Traditional';
        else
            sys1(iig,iib).label = sprintf('$$\\gamma = %.2f$$',g);
        end
        
        % Read the fractional-conformable model for a given gamma value
        sys1(iig,iib).func = @(t,x) mechsys_01(t,x,[m,b,k,g]);
        
        % Solve this problem
        [sys1(iig,iib).t,x_dummy] = ode45(sys1(iig,iib).func, tseries, [x_0, v_0]);
        
        % Store the displacement and velocity
        sys1(iig,iib).x = x_dummy(:,1);
        sys1(iig,iib).v = x_dummy(:,2);
        
        sinf = stepinfo(sys1(iig,iib).x,tseries);
        settlingtimes(iig,iib) = sinf.SettlingTime;
        peaks(iig,iib) = sinf.Peak;
        
        % Determine the entropy generation rate
        sys1(iig,iib).ds = 2*m/T_0*zeta*omega_n*sys1(iig,iib).v.^2 ;
        
        % Determine the total entropy from entropy generation rate results
        sys1(iig,iib).totalEntropy = trapz(sys1(iig,iib).ds)*diff(sys1(iig,iib).t(1:2));
        qq(iig,iib) = sys1(iig,iib).totalEntropy;
    end
    
    % Identify the minimal total entropy
    totalEntropies = [sys1(:,iib).totalEntropy];
    [minTotalEntropies(iib),idmin(iib)] = min(totalEntropies);
    
    % Only plot specific cases
    plotFlag = 'on';
    switch b
        case 0
            stringDamped = 'Undamped'; %iiCol = 1;
        case 1/4
            stringDamped = 'Underdamped'; %iiCol = 1;
        case 2
            stringDamped = 'Critically damped'; %iiCol = 3;
        case 4
            stringDamped = 'Overdamped'; %iiCol = 4;
        otherwise
            stringDamped = 'xxx';
            plotFlag = 'off';
    end
    
    zetaCell{iib} = stringDamped;
    %         if zeta == 0
    %             stringDamped = 'Undamped'; iiCol = 1;
    %         elseif zeta < 1
    %             stringDamped = 'Underdamped'; iiCol = 2;
    %         elseif zeta == 1
    %             stringDamped = 'Critically damped'; iiCol = 3;
    %         else
    %             stringDamped = 'Overdamped'; iiCol = 4;
    %         end
    
    % Plot just the specific cases
    %     if plotFlag == true
    %         plot(gamma,totalEntropies,'Color',Colors(iib,:),...
    %             'LineWidth',1,'DisplayName',...
    %             sprintf('$$\\zeta = %.3f$$, %s,',zeta,stringDamped),...
    %             'HandleVisibility',plotFlag);
    %         plot(gamma(idmin(iib)),minTotalEntropies(iib),'*r',...
    %             'HandleVisibility','off');
    %     end
end
%%
% Plot information
Graph4 = Graphics(['Sys1_totalentropy_varying_zeta2']);
ax4 = axes(Graph4.objID,'NextPlot','add','Box','on');

% qq = reshape([sys1.totalEntropy],ng,[]);
[G,Z] = meshgrid(gamma,zetaseries);

s1 = surf(ax4,G,Z,qq'); colormap jet;
%cc = jet(100); colormap(cc(1:35,:));
s1.Parent.ColorScale = 'log';
s1.Parent.ZScale = 'log';
s1.EdgeColor = 'none';
s1.FaceAlpha = 1;
% contour(ax4,G,Z,(qq)'); %clabel(C,h1);

% light('Position',[.8 1.5 10]), view([130 24])
c1 = colorbar(); c1.FontSize = 14; c1.TickLabelInterpreter = 'latex'; view(2);


% set(ax4,'ColorScale','log','ZScale','log');

z_TotalEntropies = qq(:,zetaseries == 1)';
x_gamma = gamma;
y_zeta = zetaseries(zetaseries == 1)*ones(size(x_gamma));

plot3(x_gamma,y_zeta,z_TotalEntropies,'-','Color','r',... %[0.1 0.75 0.3]
    'LineWidth',1,'DisplayName','Minimal');
% ylim([0 max(totalEntropies)]);

% leg4 = legend(ax4,'show');
% set(leg4,'interpreter','latex','location','best','box','off','NumColumns',1);
xlabel(ax4,strGamma,'Interpreter','LaTeX');
ylabel(ax4,strZ,'Interpreter','LaTeX');
zlabel(ax4,strTS,'Interpreter','LaTeX');
setall(Graph4,3,[1 1.1],fontsize,1);
zlim([1e-4 1e4])
%print(Graph4.objID,Graph4.fileName,'-r360','-djpeg','-noui');

%% Plot x-v in 2D
for iig = 1 : ng
    
    g = gamma(iig); % Get the gamma value
    sg = sprintf('%.3f',g);
    sg(sg=='.') = 'd';
    
    % Plot for phase plane
    Graph6 = Graphics(['Sys1_PhasePlane_gamma',sg]);
    ax6 = axes(Graph6.objID,'NextPlot','add','Box','on');
    
    for iib = 1 : numel(bseries)
        % Read data from b and zeta series
        b = bseries(iib);
        zeta = zetaseries(iib);
        
        plotFlag = true;
        switch b
            case 0
                stringDamped = sprintf('$$\\zeta=%d$$',zeta);%'Undamped';
                iiCol = 1;
            case 1/4
                stringDamped = sprintf('$$\\zeta=%.3f$$',zeta);%'Underdamped';
                iiCol = 2;
            case 2
                stringDamped = sprintf('$$\\zeta=%d$$',zeta);%'Critically damped';
                iiCol = 3;
            case 4
                stringDamped = sprintf('$$\\zeta=%d$$',zeta);%'Overdamped';
                iiCol = 4;
            otherwise
                stringDamped = sprintf('$$\\zeta=%.3f$$',zeta);%'xxx';
                plotFlag = false;
        end
        
        if plotFlag == true
            X = sys1(iig,iib).x;
            V = sys1(iig,iib).v;
            plot(ax6,X,V,'LineWidth',1, 'Color',Colors(iiCol,:),...
                'LineStyle','-','DisplayName',stringDamped),
        end
    end
    xlim([-2 2]), ylim([-2 2])
    xlabel(ax6,strX,'Interpreter','LaTeX');
    ylabel(ax6,strV,'Interpreter','LaTeX');
    setall(Graph6,3,[1,1],16,1);
    if iig == ng
        leg6 = legend(ax6,'show'); axis square,
        set(leg6,'interpreter','latex','location','best','box','off',...
            'NumColumns',2,'FontSize',12);
    end
    
    
end

%% new figure about settling times
%Plot for phase plane
Graph1 = Graphics(['Sys1_SettlingTimes1']);
ax1 = axes(Graph1.objID,'NextPlot','add','Box','on','YScale','log');

plot(zetaseries,settlingtimes','LineWidth',1), ylim([0,150])
xlabel(ax1,strZ,'Interpreter','LaTeX');
leg1 = legend({sys1(:,1).label});

colormap(ax1,Colors);
set(leg1,'interpreter','latex','location','best','box','off',...
    'NumColumns',2,'FontSize',fontsize); 
ylabel(ax1,'Settling Time, $$t_s$$~[s]','Interpreter','LaTeX');
setall(Graph1,2,[3,2.5],16,1);

ts_min = @(g,z)  11.33 -1.474*g -7.967*z;
[ts_mins1,iz] = min(settlingtimes,[],2);
ts_mins2 = ts_min(gamma,zetaseries(iz));

%% Plot for phase plane
% Graph1 = Graphics(['Sys1_SettlingTimes2']);
% ax1 = axes(Graph1.objID,'NextPlot','add','Box','on','YScale','log');
% 
% plot(gamma,settlingtimes(:,2:end),'LineWidth',1),
% xlim([0.65,1])
% xlabel(ax1,strGamma,'Interpreter','LaTeX');
% leg1 = legend(zetaCell(2:end));
% 
% colormap(ax1,Colors);
% set(leg1,'interpreter','latex','location','best','box','off',...
%     'NumColumns',1,'FontSize',fontsize); 
% ylabel(ax1,'Settling Time, $$t_s$$~[s]','Interpreter','LaTeX');
% setall(Graph1,2,[3,2.5],16,1);