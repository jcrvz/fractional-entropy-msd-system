%% First simulation

% Set parameters default parameters

par.k   = 1; % spring constant [N/m]
par.m   = 1; % mass [kg]
par.tau = 2; % relaxation parameter

par.omega0  = sqrt(par.k/par.m); % natural frequency
par.zeta    = par.m/2/par.tau/par.omega0; % damping factor
par.gamma   = 0.9; % Fractional order

par.F0      = 0;     % External force [N]
par.omega   = 1.1*par.omega0; % External frequency

par.tf      = 20; % final time in s
par.nt      = 500; % number time samples

ic.x0       = 1; % initial displacement [m]
ic.v0       = 0.0; % initial velocity [m/s]

% Solve the problem for a simple case
[T,X] = SolveModelComf(par,ic);

% First test
rootname    = 'first_test';

f11 = Graphics([rootname,'-tx']);
plot(T,X(:,1),'b-',T,X(:,2),'r-','linewidth',1.5);
ylabel('Magnitude'), xlabel('$$t$$ [s]')
legend('$$x(t)$$ [m]','$$v(t)$$ [m/s]','location','northeast');
ylim([-1.5 1.5]);
setup(f11); setsize(f11,2,[9,4])
%save(f11);

f12 = Graphics([rootname,'-xx']);
plot(X(:,1), X(:,2),'b','linewidth',1.5);
ylabel('$$v(t)$$ [m/s]'), xlabel('$$x(t)$$ [m]')
%axis([-1 1 -1 1]*1.5)
setup(f12); setsize(f12,2,[1,1])
%save(f12,'Responses_xx');

%% Varying fractional order

ngammas = 9; % number of gamma values
gamma_vals = linspace(0.6,1,ngammas); % values for gamma

x_vals = zeros(par.nt,ngammas); % initialise displacement's data
v_vals = zeros(par.nt,ngammas); % initialise velocity's data
l_vals = cell(1,ngammas); % initialise cell for legends

for gg = 1 : ngammas
    
    % Set the gamma's value
    par.gamma = gamma_vals(gg);
    
    % Solve the problem
    [T,X] = SolveModelComf(par,ic);
    
    % Store new data
    x_vals(:,gg) = X(:,1);
    v_vals(:,gg) = X(:,2);
    l_vals{1,gg} = sprintf('$$\\gamma = %.2f$$',par.gamma);
end

% Plot displacements
% figure('Color','w','name','Displacements varying fractional order'),
% hp = plot(T,x_vals,'linewidth',1.5);
% set(hp,{'color'}, num2cell(jet(ngammas), 2));
% set(gca,'linewidth',1.5,'ticklabelinterpreter','latex','fontsize',20);
% xlabel('Time, $$t$$ [s]','interpreter','latex','fontsize',20);
% ylabel('$$x(t)$$ [m]','interpreter','latex','fontsize',20);
% legend(l_vals,'location','best','interpreter','latex','fontsize',20);

%colors = circshift(parula(ngammas),0,2);
colors = [  linspace(0,1,ngammas)',...
            zeros(ngammas,1),...
            linspace(1,0,ngammas)'];

% Plot Entropy generation rate
rootname = 'VarGamma';

f21 = Graphics([rootname,'-tx']);
hp = plot(T,x_vals,'linewidth',1.5);
set(hp,{'color'}, num2cell(colors, 2));
ylabel('$$x(t)$$ [m]','interpreter','latex');
xlabel('$$t$$ [s]','Interpreter','latex');
legend(l_vals,'interpreter','latex','location','best');
setup(f21); setsize(f21,1,[9,4]);

% Plot x-vs.-dx
f22 = Graphics([rootname,'-xx']);
% for n = 1 : ngammas
%     plot3(repmat(gamma_vals(n),1,par.nt),x_vals(:,n),v_vals(:,n),...
%         'linewidth',1.5,'color',colors(n,:)),
%     hold on,
% end
hp = plot(x_vals,v_vals,'linewidth',1.5);
set(hp,{'color'}, num2cell(colors, 2));
ylabel('$$v(t)$$ [m/s]'),
xlabel('$$x(t)$$ [m]','interpreter','latex');
legend(l_vals,'interpreter','latex','location','best');
%set(gca,'XTick',gamma_vals,'XGrid','on'); box on,
%xlim(minmax(gamma_vals) + max(0.01*minmax(gamma_vals)).*[-1 1])
%set(gca,'YDir','reverse');%view(140,40)
setup(f22); setsize(f22,2,[1,1]);
%save(f22);

% Plot x-vs.-dx-vs.-b
f23 = Graphics([rootname,'-ptx']);
for n = 1 : ngammas
    plot3(repmat(gamma_vals(n),1,par.nt),T,x_vals(:,n),'linewidth',1.5,...
        'color',colors(n,:)),
    hold on,
end
xlabel('$$\gamma$$','interpreter','latex');
zlabel('$$x(t)$$ [m]','interpreter','latex');
ylabel('$$t$$ [s]','Interpreter','latex');
set(gca,'XTick',gamma_vals,'XGrid','on'); box on,
xlim(minmax(gamma_vals) + max(0.01*minmax(gamma_vals)).*[-1 1])
set(gca,'YDir','reverse');%view(140,40)
setup(f23); setsize(f23,2,[1,1]);
%save(f24);

%

