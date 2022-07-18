% ---------------------------------------------------%
%             Orbit-Transfer Problem                  %
% ---------------------------------------------------%
% Solve the following optimal control problem:       %
% Minimize t_f                                       %
% subject to the differential equation constraints   %
%   dr/dt = v_r                                      %
%   d\theta/dt = v_\theta/r                          %
%   dv_r/dt = v_\theta^2/r-\mu/r^2+a u_1             %
%   dv_\theta/dt = -v_rv_\theta/r+a u_2/m            %
%   dm/dt = - a/v_e
% the equality path constraint                       %
%   u_1^2 + u_2^2 = 1                                %
% and the boundary conditions                        %
%   r(0) = 1                                         %
%   \theta(0) = 0                                    %
%   v_r(0) = 0                                       %
%   v_\theta(0) = sqrt(mu/r(0))                      %
%   m(0) = 1                                         %
%   r(t_f) = 1.5                                     %
%   v_r(t_f) = 0                                     %
%   v_\theta(t_f) = sqrt(mu/r(t_f)                   %
% -------------------------------------------------- %
% BEGIN: DO NOT ALTER THE FOLLOWING LINES OF CODE!!! %
% -------------------------------------------------- %
global igrid CONSTANTS psStuff nstates ncontrols npaths
% -------------------------------------------------- %
% END:   DO NOT ALTER THE FOLLOWING LINES OF CODE!!! %
% -------------------------------------------------- %

CONSTANTS.MU = 1;
CONSTANTS.m0 = 1;
CONSTANTS.mdot = 0.0749;
CONSTANTS.ve = 1.8758344;
nstates = 5;
ncontrols = 3;
npaths = 1;

% Bounds on State and Control
r0 = 1; theta0 = 0; vr0 = 0; vtheta0 = sqrt(CONSTANTS.MU/r0); m0 = 1;
rf = 1.5; vrf = 0;  vthetaf = sqrt(CONSTANTS.MU/rf);
rmin = 0.5; rmax = 2;
thetamin = 0; thetamax = 4*pi;
vrmin = -50; vrmax = 50;
vthetamin = -50; vthetamax = 50;
mmin = 0.5; mmax = m0;
u1min = -10; u1max = 10;
u2min = -10; u2max = 10;
Tmin = 0; Tmax = 0.1405;
t0min = 0; t0max = 0;
tfmin = 0; tfmax = 50;

%-----------------------------------------------------------------%
%      Compute Points, Weights, and Differentiation Matrix        %
%-----------------------------------------------------------------%
%-----------------------------------------------------------------%
% Choose Polynomial Degree and Number of Mesh Intervals           %
% numIntervals = 1 ===> p-method                                  %
% numIntervals > 1 ===> h-method                                  %
%-----------------------------------------------------------------%
N = 4;
numIntervals = 32;
%-----------------------------------------------------------------%
% DO NOT ALTER THE LINE OF CODE SHOWN BELOW!                      %
%-----------------------------------------------------------------%
meshPoints = linspace(-1,1,numIntervals+1).';  
polyDegrees = N*ones(numIntervals,1);
[tau,w,D] = lgrPS(meshPoints,polyDegrees);
psStuff.tau = tau; psStuff.w = w; psStuff.D = D; NLGR = length(w);
%-----------------------------------------------------------------%
% DO NOT ALTER THE LINES OF CODE SHOWN ABOVE!                     %
%-----------------------------------------------------------------%

% Set the bounds on the NLP variables.
zrmin = rmin*ones(length(tau),1);
zrmax = rmax*ones(length(tau),1);
zrmin(1) = r0; zrmax(1) = r0;
zrmin(end) = rf; zrmax(end) = rf;

zthetamin = thetamin*ones(length(tau),1);
zthetamax = thetamax*ones(length(tau),1);
zthetamin(1) = theta0; zthetamax(1) = theta0;

zvrmin = vrmin*ones(length(tau),1);
zvrmax = vrmax*ones(length(tau),1);
zvrmin(1) = vr0; zvrmax(1) = vr0;
zvrmin(end) = vrf; zvrmax(end) = vrf;

zvthetamin = vthetamin*ones(length(tau),1);
zvthetamax = vthetamax*ones(length(tau),1);
zvthetamin(1) = vtheta0; zvthetamax(1) = vtheta0;
zvthetamin(end) = vthetaf; zvthetamax(end) = vthetaf;

zmmin = mmin*ones(length(tau),1);
zmmax = mmax*ones(length(tau),1);
zmmin(1) = mmax; zmmax(1) = mmax;

zu1min = u1min*ones(length(tau)-1,1);
zu1max = u1max*ones(length(tau)-1,1);

zu2min = u2min*ones(length(tau)-1,1);
zu2max = u2max*ones(length(tau)-1,1);

zTmin = Tmin*ones(length(tau)-1,1);
zTmax = Tmax*ones(length(tau)-1,1);

zmin = [zrmin; zthetamin; zvrmin; zvthetamin; zmmin; zu1min; zu2min; zTmin; t0min; tfmin];
zmax = [zrmax; zthetamax; zvrmax; zvthetamax; zmmax; zu1max; zu2max; zTmax; t0max; tfmax];

% Set the bounds on the NLP constraints
% There are NSTATES sets of defect constraints.
defectMin = zeros(nstates*(length(tau)-1),1);
defectMax = zeros(nstates*(length(tau)-1),1);
% There is one path constraint
pathMin = ones(length(tau)-1,1); pathMax = ones(length(tau)-1,1);
% There is one nonlinear event constraint
bcMin = 0; bcMax = 0;
objMin = -inf; objMax = inf;
Fmin = [objMin; defectMin; pathMin];
Fmax = [objMax; defectMax; pathMax];

% Supply an initial guess
rguess = linspace(r0,rf,NLGR+1).';
thetaguess = linspace(theta0,theta0,NLGR+1).';
vrguess = linspace(vr0,vrf,NLGR+1).';
vthetaguess = linspace(vtheta0,vthetaf,NLGR+1).';
mguess = linspace(mmax,mmin,NLGR+1).';
u1guess = linspace(1,1,NLGR).';
u2guess = linspace(0,0,NLGR).';
Tguess = linspace(Tmax,Tmax,NLGR).';
t0guess = 0;
tfguess = 3.2481;
z0 = [rguess;thetaguess;vrguess;vthetaguess;mguess;u1guess;u2guess;Tguess;t0guess;tfguess];

%-----------------------------------------------------------------%
% Generate derivatives and sparsity pattern using Adigator        %
%-----------------------------------------------------------------%
% - Constraint Function Derivatives
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('OrbitTransferFun',{x});
S_jac  = output.JacobianStructure;
[iGfun,jGvar] = find(S_jac);

% - Objective Function Derivatives
xsize  = size(z0);
x      = adigatorCreateDerivInput(xsize,'z0');
output = adigatorGenJacFile('OrbitTransferObj',{x});
grd_structure = output.JacobianStructure;

%-----------------------------------------------------------------%
% set IPOPT callback functions
%-----------------------------------------------------------------%
funcs.objective   = @(Z)OrbitTransferObj(Z);
funcs.gradient    = @(Z)OrbitTransferGrd(Z);
funcs.constraints = @(Z)OrbitTransferCon(Z);
funcs.jacobian    = @(Z)OrbitTransferJac(Z);
funcs.jacobianstructure = @()OrbitTransferJacPat(S_jac);
options.ipopt.hessian_approximation = 'limited-memory';

%-----------------------------------------------------------------%
% Set IPOPT Options %
%-----------------------------------------------------------------%
options.ipopt.tol = 1e-8;
options.ipopt.linear_solver = 'ma57';
options.ipopt.max_iter = 2000;
options.ipopt.mu_strategy = 'adaptive';
options.ipopt.ma57_automatic_scaling = 'yes';
options.ipopt.print_user_options = 'yes';
options.ipopt.output_file = ['OrbitTransfer','IPOPTinfo.txt']; % print output file
options.ipopt.print_level = 5; % set print level default

options.lb = zmin; % Lower bound on the variables.
options.ub = zmax; % Upper bound on the variables.
options.cl = Fmin; % Lower bounds on the constraint functions.
options.cu = Fmax; % Upper bounds on the constraint functions.

%-----------------------------------------------------------------%
% Call IPOPT
%-----------------------------------------------------------------%
[z, info] = ipopt(z0,funcs,options);

%-----------------------------------------------------------------%
% extract lagrange multipliers from ipopt output, info
%-----------------------------------------------------------------%
Fmul = info.lambda;

% Extract the state and control from the decision vector z.
% Remember that the state is approximated at the LGR points
% plus the final point, while the control is only approximated 
% at only the LGR points.
r = z(1:NLGR+1);
theta = z(NLGR+2:2*(NLGR+1));
vr = z(2*(NLGR+1)+1:3*(NLGR+1));
vtheta = z(3*(NLGR+1)+1:4*(NLGR+1));
m = z(4*(NLGR+1)+1:5*(NLGR+1));
u1 = z(5*(NLGR+1)+1:5*(NLGR+1)+NLGR);
u2 = z(5*(NLGR+1)+NLGR+1:5*(NLGR+1)+2*NLGR);
T = z(5*(NLGR+1)+2*NLGR+1:5*(NLGR+1)+3*NLGR);
beta = 180/pi*atan2(u1,u2);
t0 = z(end-1);
tf = z(end);
t = (tf-t0)*(tau+1)/2+t0;
tLGR = t(1:end-1);
%-----------------------------------------------------------------%
% Extract the Lagrange multipliers corresponding                  %
% the defect constraints.                                         %
%-----------------------------------------------------------------%
multipliersDefects = Fmul(2:nstates*NLGR+1);
multipliersDefects = reshape(multipliersDefects,NLGR,nstates);
%-----------------------------------------------------------------%
% Compute the costates at the LGR points via transformation       %
%-----------------------------------------------------------------%
costateLGR = inv(diag(w))*multipliersDefects;
%-----------------------------------------------------------------%
% Compute the costate at the tau=+1 via transformation            %
%-----------------------------------------------------------------%
costateF = D(:,end).'*multipliersDefects;
%-----------------------------------------------------------------%
% Now assemble the costates into a single matrix                  %
%-----------------------------------------------------------------%
costate = [costateLGR; costateF];    
lamr = costate(:,1); lamtheta = costate(:,2);
lamvr = costate(:,3); lamvtheta = costate(:,4);
lamm = costate(:,5);

%-------------%
% Plot Results 
%-------------%
% figure(1);
% title('Plot of all the states');
% subplot(3,2,1);
% plot(t,r,'g-');
% xl = xlabel('$t$','Interpreter','LaTeX');
% yl = ylabel('$r(t)$','Interpreter','LaTeX');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontName','Times','FontSize',18);
% grid on;
% 
% subplot(3,2,2);
% plot(t,theta,'b-');
% xl = xlabel('$t$','Interpreter','LaTeX');
% yl = ylabel('$\theta(t)$','Interpreter','LaTeX');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontName','Times','FontSize',18);
% grid on;
% 
% subplot(3,2,3);
% plot(t,vr,'c-');
% xl = xlabel('$t$','Interpreter','LaTeX');
% yl = ylabel('$v_r(t)$','Interpreter','LaTeX');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontName','Times','FontSize',18);
% grid on;
% 
% subplot(3,2,4);
% plot(t,vtheta,'m-');
% xl = xlabel('$t$','Interpreter','LaTeX');
% yl = ylabel('$v_\theta(t)$','Interpreter','LaTeX');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontName','Times','FontSize',18);
% grid on;
% 
% subplot(3,2,5);
% plot(t,m,'k-');
% xl = xlabel('$t$','Interpreter','LaTeX');
% yl = ylabel('$v_\theta(t)$','Interpreter','LaTeX');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontName','Times','FontSize',18);
% grid on;

figure(1)
plot(t,r,'r',t,theta,'b',t,vr,'m',t,vtheta,'k', t, m,'g');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('states','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
title('Plot of States');
set(gca,'FontName','Times','FontSize',18);
grid on;

figure(2);
plot(tLGR,beta,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\beta(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
title('Throttle angle beta');
grid on;

figure(3);
subplot(3,1,1);
plot(tLGR,u1,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_1(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

subplot(3,1,2);
plot(tLGR,u2,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_2(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

subplot(3,1,3);
plot(tLGR,T,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$T(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

% figure(4);
% subplot(2,2,1);
% plot(t,lamr,'-o');
% xl = xlabel('$t$','Interpreter','LaTeX');
% yl = ylabel('$\lambda_r(t)$','Interpreter','LaTeX');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontName','Times','FontSize',18);
% grid on;
% 
% subplot(2,2,2);
% plot(t,lamtheta,'-o');
% xl = xlabel('$t$','Interpreter','LaTeX');
% yl = ylabel('$\lambda_\theta(t)$','Interpreter','LaTeX');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontName','Times','FontSize',18);
% grid on;
% 
% subplot(2,2,3);
% plot(t,lamvr,'-o');
% xl = xlabel('$t$','Interpreter','LaTeX');
% yl = ylabel('$\lambda_{v_r}(t)$','Interpreter','LaTeX');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontName','Times','FontSize',18);
% grid on;
% 
% subplot(2,2,4);
% plot(t,lamvtheta,'-o');
% xl = xlabel('$t$','Interpreter','LaTeX');
% yl = ylabel('$\lambda_{v_\theta}(t)$','Interpreter','LaTeX');
% set(xl,'FontSize',18);
% set(yl,'FontSize',18);
% set(gca,'FontName','Times','FontSize',18);
% grid on;

figure(4);
plot(t,lamr,'r',t,lamtheta,'b',t,lamvr,'m',t,lamvtheta,'k', t, lamm,'g');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('Costates','Interpreter','LaTeX');
title('Plot of CoStates');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
grid on;

