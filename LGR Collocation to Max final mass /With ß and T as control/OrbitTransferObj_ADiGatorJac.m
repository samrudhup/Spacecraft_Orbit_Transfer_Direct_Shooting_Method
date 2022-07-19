% This code was generated using ADiGator version 1.4
% ©2010-2014 Matthew J. Weinstein and Anil V. Rao
% ADiGator may be obtained at https://sourceforge.net/projects/adigator/ 
% Contact: mweinstein@ufl.edu
% Bugs/suggestions may be reported to the sourceforge forums
%                    DISCLAIMER
% ADiGator is a general-purpose software distributed under the GNU General
% Public License version 3.0. While the software is distributed with the
% hope that it will be useful, both the software and generated code are
% provided 'AS IS' with NO WARRANTIES OF ANY KIND and no merchantability
% or fitness for any purpose or application.

function obj = OrbitTransferObj_ADiGatorJac(z)
global ADiGator_OrbitTransferObj_ADiGatorJac
if isempty(ADiGator_OrbitTransferObj_ADiGatorJac); ADiGator_LoadData(); end
Gator1Data = ADiGator_OrbitTransferObj_ADiGatorJac.OrbitTransferObj_ADiGatorJac.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: % Computes the objective function of the problem
global psStuff nstates ncontrols CONSTANTS 
%User Line: global
%User Line: %-----------------------------------------------------------------%
%User Line: %         Extract the constants used in the problem.              %
%User Line: %-----------------------------------------------------------------%
MU = CONSTANTS.MU;
%User Line: MU = CONSTANTS.MU;
mdot = CONSTANTS.mdot;
%User Line: mdot = CONSTANTS.mdot;
ve = CONSTANTS.ve;
%User Line: ve = CONSTANTS.ve;
%User Line: %-----------------------------------------------------------------%
%User Line: % Radau pseudospectral method quantities required:                %
%User Line: %   - Differentiation matrix (psStuff.D)                          %
%User Line: %   - Legendre-Gauss-Radau weights (psStuff.w)                    %
%User Line: %   - Legendre-Gauss-Radau points (psStuff.tau)                   %
%User Line: %-----------------------------------------------------------------%
D = psStuff.D;
%User Line: D = psStuff.D;
tau = psStuff.tau;
%User Line: tau = psStuff.tau;
w = psStuff.w;
%User Line: w = psStuff.w;
%User Line: %-----------------------------------------------------------------%
%User Line: % Decompose the NLP decision vector into pieces containing        %
%User Line: %    - the state                                                  %
%User Line: %    - the control                                                %
%User Line: %    - the initial time                                           %
%User Line: %    - the final time                                             %
%User Line: %-----------------------------------------------------------------%
cada1f1 = length(tau);
N.f = cada1f1 - 1;
%User Line: N = length(tau)-1;
cada1f1 = N.f + 1;
cada1f2 = nstates*cada1f1;
stateIndices.f = 1:cada1f2;
%User Line: stateIndices = 1:nstates*(N+1);
cada1f1 = N.f + 1;
cada1f2 = nstates*cada1f1;
cada1f3 = cada1f2 + 1;
cada1f4 = N.f + 1;
cada1f5 = nstates*cada1f4;
cada1f6 = ncontrols*N.f;
cada1f7 = cada1f5 + cada1f6;
controlIndices.f = cada1f3:cada1f7;
%User Line: controlIndices = (nstates*(N+1)+1):(nstates*(N+1)+ncontrols*N);
cada1f1 = length(controlIndices.f);
cada1f2 = controlIndices.f(cada1f1);
t0Index.f = cada1f2 + 1;
%User Line: t0Index = controlIndices(end)+1;
tfIndex.f = t0Index.f + 1;
%User Line: tfIndex = t0Index+1;
stateVector.dz0 = z.dz0(Gator1Data.Index1);
stateVector.f = z.f(stateIndices.f);
%User Line: stateVector = z(stateIndices);
controlVector.dz0 = z.dz0(Gator1Data.Index2);
controlVector.f = z.f(controlIndices.f);
%User Line: controlVector = z(controlIndices);
t0.dz0 = z.dz0(48);
t0.f = z.f(t0Index.f);
%User Line: t0 = z(t0Index);
tf.dz0 = z.dz0(49);
tf.f = z.f(tfIndex.f);
%User Line: tf = z(tfIndex);
cada1td1 = zeros(2,1);
cada1td1(2) = tf.dz0;
cada1td1(1) = cada1td1(1) + -t0.dz0;
cada1f1dz0 = cada1td1;
cada1f1 = tf.f - t0.f;
cada1f2 = tau + 1;
cada1tempdz0 = cada1f1dz0(Gator1Data.Index3);
cada1tf1 = cada1f2(Gator1Data.Index5);
cada1f3dz0 = cada1tf1(:).*cada1tempdz0(Gator1Data.Index4);
cada1f3 = cada1f1*cada1f2;
cada1f4dz0 = cada1f3dz0./2;
cada1f4 = cada1f3/2;
cada1tempdz0 = t0.dz0(Gator1Data.Index6);
cada1td1 = zeros(13,1);
cada1td1(Gator1Data.Index7) = cada1f4dz0;
cada1td1(Gator1Data.Index8) = cada1td1(Gator1Data.Index8) + cada1tempdz0;
t.dz0 = cada1td1;
t.f = cada1f4 + t0.f;
%User Line: t = (tf-t0)*(tau+1)/2+t0;
cada1f1 = length(t.f);
cada1f2 = cada1f1 - 1;
cada1f3 = 1:cada1f2;
tLGR.dz0 = t.dz0(Gator1Data.Index9);
tLGR.f = t.f(cada1f3);
%User Line: tLGR = t(1:end-1);
%User Line: %-----------------------------------------------------------------%
%User Line: % Reshape the state and control parts of the NLP decision vector  %
%User Line: % to matrices of sizes (N+1) by nstates and (N+1) by ncontrols,   %
%User Line: % respectively.  The state is approximated at the N LGR points    %
%User Line: % plus the final point.  Thus, each column of the state vector is %
%User Line: % length N+1.  The LEFT-HAND SIDE of the defect constraints, D*X, %
%User Line: % uses the state at all of the points (N LGR points plus final    %
%User Line: % point).  The RIGHT-HAND SIDE of the defect constraints,         %
%User Line: % (tf-t0)F/2, uses the state and control at only the LGR points.  %
%User Line: % Thus, it is necessary to extract the state approximations at    %
%User Line: % only the N LGR points.  Finally, in the Radau pseudospectral    %
%User Line: % method, the control is approximated at only the N LGR points.   %
%User Line: %-----------------------------------------------------------------%
cada1f1 = N.f + 1;
statePlusEnd.dz0 = stateVector.dz0;
statePlusEnd.f = reshape(stateVector.f,cada1f1,nstates);
%User Line: statePlusEnd   = reshape(stateVector,N+1,nstates);
cada1f1 = size(statePlusEnd.f,1);
cada1f2 = cada1f1 - 1;
cada1f3 = 1:cada1f2;
stateLGR.dz0 = statePlusEnd.dz0(Gator1Data.Index10);
stateLGR.f = statePlusEnd.f(cada1f3,:);
%User Line: stateLGR = statePlusEnd(1:end-1,:);
control.dz0 = controlVector.dz0;
control.f = reshape(controlVector.f,N.f,ncontrols);
%User Line: control = reshape(controlVector,N,ncontrols);
%User Line: %-----------------------------------------------------------------%
%User Line: % Identify the components of the state column-wise from stateLGR. %
%User Line: %-----------------------------------------------------------------%
r.dz0 = stateLGR.dz0(Gator1Data.Index11);
r.f = stateLGR.f(:,1);
%User Line: r = stateLGR(:,1);
theta.dz0 = stateLGR.dz0(Gator1Data.Index12);
theta.f = stateLGR.f(:,2);
%User Line: theta = stateLGR(:,2);
vr.dz0 = stateLGR.dz0(Gator1Data.Index13);
vr.f = stateLGR.f(:,3);
%User Line: vr = stateLGR(:,3);
vtheta.dz0 = stateLGR.dz0(Gator1Data.Index14);
vtheta.f = stateLGR.f(:,4);
%User Line: vtheta = stateLGR(:,4);
m.dz0 = stateLGR.dz0(Gator1Data.Index15);
m.f = stateLGR.f(:,5);
%User Line: m = stateLGR(:,5);
u1.dz0 = control.dz0(Gator1Data.Index16);
u1.f = control.f(:,1);
%User Line: u1 = control(:,1);
T.dz0 = control.dz0(Gator1Data.Index17);
T.f = control.f(:,2);
%User Line: T = control(:,2);
%User Line: % T = control(:,3);
%User Line: %-----------------------------------------------------------------%
%User Line: % The quantity STATEF is the value of the state at the final      %
%User Line: % time, tf, which corresponds to the state at $\tau=1$.           %
%User Line: %-----------------------------------------------------------------%
cada1f1 = size(statePlusEnd.f,1);
stateF.dz0 = statePlusEnd.dz0(Gator1Data.Index18);
stateF.f = statePlusEnd.f(cada1f1,:);
%User Line: stateF = statePlusEnd(end,:);
%User Line: %-----------------------------------------------------------------%
%User Line: % (Edit) The orbit-transfer problem contains one nonlinear boundary%
%User Line: % condition $\sqrt{mu/r(t_f)-v_\theta(t_f) = 0$.  Because $r(t)$  %
%User Line: % and $v_\theta(t)$ are the first and fourth components of the    %
%User Line: % state, it is necessary to extract stateF(1) and stateF(4) in    %
%User Line: % order to compute this boundary condition function.              %
%User Line: %-----------------------------------------------------------------%
tF.dz0 = tf.dz0; tF.f = tf.f;
%User Line: tF = tf;
%User Line: % Cost Function
cada1f1dz0 = stateF.dz0(5);
cada1f1 = stateF.f(5);
J.dz0 = -cada1f1dz0;
J.f = uminus(cada1f1);
%User Line: J   = -stateF(5);
obj.dz0 = J.dz0; obj.f = J.f;
%User Line: obj = J;
obj.dz0_size = 49;
obj.dz0_location = Gator1Data.Index19;
end


function ADiGator_LoadData()
global ADiGator_OrbitTransferObj_ADiGatorJac
ADiGator_OrbitTransferObj_ADiGatorJac = load('OrbitTransferObj_ADiGatorJac.mat');
return
end