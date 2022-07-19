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
ncontrols = 2;
npaths = 0;

% Bounds on State and Control
r0 = 1; theta0 = 0; vr0 = 0; vtheta0 = sqrt(CONSTANTS.MU/r0); m0 = 1;
rf = 1.5; vrf = 0;  vthetaf = sqrt(CONSTANTS.MU/rf);
rmin = 0.5; rmax = 2;
thetamin = 0; thetamax = 4*pi;
vrmin = -30; vrmax = 50;
vthetamin = -30; vthetamax = 50;
mmin = 0.1; mmax = m0;
u1min = -3; u1max = 3;
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

% zu2min = u2min*ones(length(tau)-1,1);
% zu2max = u2max*ones(length(tau)-1,1);

zTmin = Tmin*ones(length(tau)-1,1);
zTmax = Tmax*ones(length(tau)-1,1);

zmin = [zrmin; zthetamin; zvrmin; zvthetamin; zmmin; zu1min; zTmin; t0min; tfmin];
zmax = [zrmax; zthetamax; zvrmax; zvthetamax; zmmax; zu1max; zTmax; t0max; tfmax];

% Set the bounds on the NLP constraints
% There are NSTATES sets of defect constraints.
defectMin = zeros(nstates*(length(tau)-1),1);
defectMax = zeros(nstates*(length(tau)-1),1);
% There is no path contraint in this problem
% pathMin = ones(length(tau)-1,1); pathMax = ones(length(tau)-1,1);
% There is one nonlinear event constraint
bcMin = 0; bcMax = 0;
objMin = tfmin; objMax = tfmax;
Fmin = [objMin; defectMin];
Fmax = [objMax; defectMax];

% Supply an initial guess
rguess = linspace(r0,rf,NLGR+1).';
thetaguess = linspace(theta0,theta0,NLGR+1).';
vrguess = linspace(vr0,vrf,NLGR+1).';
vthetaguess = linspace(vtheta0,vthetaf,NLGR+1).';
mguess = linspace(mmax,mmin,NLGR+1).';
u1guess = linspace(0,0,NLGR).';
u2guess = linspace(1,1,NLGR).';
Tguess = linspace(Tmax,Tmax,NLGR).';
t0guess = 0;
tfguess = 100;
z0 = [rguess;thetaguess;vrguess;vthetaguess;mguess;u1guess;Tguess;t0guess;tfguess];

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
T = z(5*(NLGR+1)+NLGR+1:5*(NLGR+1)+2*NLGR);
% T = z(5*(NLGR+1)+2*NLGR+1:5*(NLGR+1)+3*NLGR);
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
figure(1);
plot(t,r,'b', t, theta, 'g', t, vr, 'r', t, vtheta, 'c', t, m, 'k' );
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$r(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
title('Plot of States');
legend('r','theta','vr','vtheta','location','best');
set(gca,'FontName','Times','FontSize',18);
grid on;

figure(2);
subplot(1,2,1);
plot(tLGR,u1,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$u_1(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
title('Throttle angle beta');
grid on;

subplot(1,2,2);
plot(tLGR,T,'-o');
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$T(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
title('Control Input T');
grid on;

figure(3);
plot(t,lamr,'b', t, lamtheta, 'g', t, lamvr, 'r', t, lamvtheta, 'c', t, lamm, 'k' );
xl = xlabel('$t$','Interpreter','LaTeX');
yl = ylabel('$\lambda_r(t)$','Interpreter','LaTeX');
set(xl,'FontSize',18);
set(yl,'FontSize',18);
set(gca,'FontName','Times','FontSize',18);
legend('lamr','lamtheta','lamvr','lamvtheta','lamm','location','best');
title('Plot of CoStates');
grid on;

figure (4);
polarplot(theta,r);
title('Polar Plot of the location of the Spacecraft');

%OrbitTransferFun.m
function C = OrbitTransferFun(z)

%-----------------------------------------------------------------%
% Objective and constraint functions for the orbit-transfer       %
% problem.  This function is designed to be used with the NLP     %
% solver SNOPT.                                                   %
%-----------------------------------------------------------------%
%      DO NOT FOR ANY REASON ALTER THE LINE OF CODE BELOW!        %
global psStuff nstates ncontrols npaths CONSTANTS                 %
%      DO NOT FOR ANY REASON ALTER THE LINE OF CODE ABOVE!        %
%-----------------------------------------------------------------%

%-----------------------------------------------------------------%
%         Extract the constants used in the problem.              %
%-----------------------------------------------------------------%
MU = CONSTANTS.MU; mdot = CONSTANTS.mdot; ve = CONSTANTS.ve;

%-----------------------------------------------------------------%
% Radau pseudospectral method quantities required:                %
%   - Differentiation matrix (psStuff.D)                          %
%   - Legendre-Gauss-Radau weights (psStuff.w)                    %
%   - Legendre-Gauss-Radau points (psStuff.tau)                   %
%-----------------------------------------------------------------%
D = psStuff.D; tau = psStuff.tau; w = psStuff.w;

%-----------------------------------------------------------------%
% Decompose the NLP decision vector into pieces containing        %
%    - the state                                                  %
%    - the control                                                %
%    - the initial time                                           %
%    - the final time                                             %
%-----------------------------------------------------------------%
N = length(tau)-1;
stateIndices = 1:nstates*(N+1);
controlIndices = (nstates*(N+1)+1):(nstates*(N+1)+ncontrols*N);
t0Index = controlIndices(end)+1;
tfIndex = t0Index+1;
stateVector = z(stateIndices);
controlVector = z(controlIndices);
t0 = z(t0Index);
tf = z(tfIndex);
t = (tf-t0)*(tau+1)/2+t0;
tLGR = t(1:end-1);

%-----------------------------------------------------------------%
% Reshape the state and control parts of the NLP decision vector  %
% to matrices of sizes (N+1) by nstates and (N+1) by ncontrols,   %
% respectively.  The state is approximated at the N LGR points    %
% plus the final point.  Thus, each column of the state vector is %
% length N+1.  The LEFT-HAND SIDE of the defect constraints, D*X, %
% uses the state at all of the points (N LGR points plus final    %
% point).  The RIGHT-HAND SIDE of the defect constraints,         %
% (tf-t0)F/2, uses the state and control at only the LGR points.  %
% Thus, it is necessary to extract the state approximations at    %
% only the N LGR points.  Finally, in the Radau pseudospectral    %
% method, the control is approximated at only the N LGR points.   %
%-----------------------------------------------------------------%
statePlusEnd   = reshape(stateVector,N+1,nstates);
stateLGR = statePlusEnd(1:end-1,:);
control = reshape(controlVector,N,ncontrols);

%-----------------------------------------------------------------%
% Identify the components of the state column-wise from stateLGR. % 
%-----------------------------------------------------------------%
r = stateLGR(:,1);
theta = stateLGR(:,2);
vr = stateLGR(:,3);
vtheta = stateLGR(:,4);
m = stateLGR(:,5);
u1 = control(:,1);
T = control(:,2);


%-----------------------------------------------------------------%
% The quantity STATEF is the value of the state at the final      %
% time, tf, which corresponds to the state at $\tau=1$.           %
%-----------------------------------------------------------------%
stateF = statePlusEnd(end,:);
%-----------------------------------------------------------------%
% The orbit-raising problem contains one nonlinear boundary       %
% condition $\sqrt{mu/r(t_f)-v_\theta(t_f) = 0$.  Because $r(t)$  %
% and $v_\theta(t)$ are the first and fourth components of the    %
% state, it is necessary to extract stateF(1) and stateF(4) in    %
% order to compute this boundary condition function.              %
%-----------------------------------------------------------------%
rF = stateF(1);
vrF = stateF(3);
vthetaF = stateF(4);

% a = T./m;

%-----------------------------------------------------------------%
% Compute the right-hand side of the differential equations at    %
% the N LGR points.  Each component of the right-hand side is     %
% stored as a column vector of length N, that is each column has  %
% the form                                                        %
%                   [ f_i(x_1,u_1,t_1) ]                          %
%                   [ f_i(x_2,u_2,t_2) ]                          %
%                           .                                     %
%                           .                                     %
%                           .                                     %
%                   [ f_i(x_N,u_N,t_N) ]                          %
% where "i" is the right-hand side of the ith component of the    %
% vector field f.  It is noted that in MATLABB the calculation of %
% the right-hand side is vectorized.                              %
%-----------------------------------------------------------------%
rdot = vr;
thetadot = vtheta./r;
vrdot = vtheta.^2./r-MU./r.^2+T.*sin(u1)./m;
vthetadot = -vr.*vtheta./r+T.*cos(u1)./m;
m_dot = -T./ve;
diffeqRHS = [rdot, thetadot, vrdot, vthetadot, m_dot];

%-----------------------------------------------------------------%
% Compute the left-hand side of the defect constraints, recalling %
% that the left-hand side is computed using the state at the LGR  %
% points PLUS the final point.                                    %
%-----------------------------------------------------------------%
diffeqLHS = D*statePlusEnd;

%-----------------------------------------------------------------%
% Construct the defect constraints at the N LGR points.           %
% Remember that the right-hand side needs to be scaled by the     %
% factor (tf-t0)/2 because the rate of change of the state is     %
% being taken with respect to $\tau\in[-1,+1]$.  Thus, we have    %
% $dt/t\dau=(tf-t0)/2$.                                           %
%-----------------------------------------------------------------%
defects = diffeqLHS-(tf-t0)*diffeqRHS/2;

%-----------------------------------------------------------------%
% Construct the path constraints at the N LGR points.             %
%-----------------------------------------------------------------%
paths = 0;

%-----------------------------------------------------------------%
% Reshape the defect contraints into a column vector.             % 
%-----------------------------------------------------------------%
defects = reshape(defects,N*nstates,1);

%-----------------------------------------------------------------%
% Reshape the path contraints into a column vector.             % 
%-----------------------------------------------------------------%
% paths = reshape(paths,N*npaths,1);

%-----------------------------------------------------------------%
% Compute the nonlinear boundary condition.                       %
%-----------------------------------------------------------------%
bcs = [];

%-----------------------------------------------------------------%
% Construct the objective function plus constraint vector.        %
%-----------------------------------------------------------------%
C = [tf;defects];
end

%OrbitTransferObj.m
function obj = OrbitTransferObj(z)
% Computes the objective function of the problem

global psStuff nstates ncontrols CONSTANTS       

%-----------------------------------------------------------------%
%         Extract the constants used in the problem.              %
%-----------------------------------------------------------------%
MU = CONSTANTS.MU; mdot = CONSTANTS.mdot;  ve = CONSTANTS.ve;

%-----------------------------------------------------------------%
% Radau pseudospectral method quantities required:                %
%   - Differentiation matrix (psStuff.D)                          %
%   - Legendre-Gauss-Radau weights (psStuff.w)                    %
%   - Legendre-Gauss-Radau points (psStuff.tau)                   %
%-----------------------------------------------------------------%
D = psStuff.D; tau = psStuff.tau; w = psStuff.w;

%-----------------------------------------------------------------%
% Decompose the NLP decision vector into pieces containing        %
%    - the state                                                  %
%    - the control                                                %
%    - the initial time                                           %
%    - the final time                                             %
%-----------------------------------------------------------------%
N = length(tau)-1;
stateIndices = 1:nstates*(N+1);
controlIndices = (nstates*(N+1)+1):(nstates*(N+1)+ncontrols*N);
t0Index = controlIndices(end)+1;
tfIndex = t0Index+1;
stateVector = z(stateIndices);
controlVector = z(controlIndices);
t0 = z(t0Index);
tf = z(tfIndex);
t = (tf-t0)*(tau+1)/2+t0;
tLGR = t(1:end-1);

%-----------------------------------------------------------------%
% Reshape the state and control parts of the NLP decision vector  %
% to matrices of sizes (N+1) by nstates and (N+1) by ncontrols,   %
% respectively.  The state is approximated at the N LGR points    %
% plus the final point.  Thus, each column of the state vector is %
% length N+1.  The LEFT-HAND SIDE of the defect constraints, D*X, %
% uses the state at all of the points (N LGR points plus final    %
% point).  The RIGHT-HAND SIDE of the defect constraints,         %
% (tf-t0)F/2, uses the state and control at only the LGR points.  %
% Thus, it is necessary to extract the state approximations at    %
% only the N LGR points.  Finally, in the Radau pseudospectral    %
% method, the control is approximated at only the N LGR points.   %
%-----------------------------------------------------------------%
statePlusEnd   = reshape(stateVector,N+1,nstates);
stateLGR = statePlusEnd(1:end-1,:);
control = reshape(controlVector,N,ncontrols);

%-----------------------------------------------------------------%
% Identify the components of the state column-wise from stateLGR. % 
%-----------------------------------------------------------------%
r = stateLGR(:,1);
theta = stateLGR(:,2);
vr = stateLGR(:,3);
vtheta = stateLGR(:,4);
m = stateLGR(:,5);
u1 = control(:,1);
T = control(:,2);

%-----------------------------------------------------------------%
% The quantity STATEF is the value of the state at the final      %
% time, tf, which corresponds to the state at $\tau=1$.           %
%-----------------------------------------------------------------%
stateF = statePlusEnd(end,:);
%-----------------------------------------------------------------%
% (Edit) The orbit-transfer problem contains one nonlinear boundary%
% condition $\sqrt{mu/r(t_f)-v_\theta(t_f) = 0$.  Because $r(t)$  %
% and $v_\theta(t)$ are the first and fourth components of the    %
% state, it is necessary to extract stateF(1) and stateF(4) in    %
% order to compute this boundary condition function.              %
%-----------------------------------------------------------------%
tF = tf;

% Cost Function
J   = tF;
obj = J;

end

%OrbitTransferCon.m
function constraints = OrbitTransferCon(Z)
% computes the constraints

output      = OrbitTransferFun(Z);
constraints = output;

end

%OrbitTransferGrd.m
function grd = OrbitTransferGrd(Z)
% computes the gradient

output = OrbitTransferObj_Jac(Z);
grd    = output;

end

%OrbitTransferJac.m
function jac = OrbitTransferJac(Z)
% computes the jacobian

[jac,~] = OrbitTransferFun_Jac(Z);

end

%OrbitTransferJacPat.m
function jacpat = OrbitTransferJacPat(S_jac)
% computes the jacobian structure

jacpat = S_jac;

end



