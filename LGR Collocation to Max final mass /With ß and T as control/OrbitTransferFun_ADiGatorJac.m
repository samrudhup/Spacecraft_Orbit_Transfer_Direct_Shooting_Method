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

function C = OrbitTransferFun_ADiGatorJac(z)
global ADiGator_OrbitTransferFun_ADiGatorJac
if isempty(ADiGator_OrbitTransferFun_ADiGatorJac); ADiGator_LoadData(); end
Gator1Data = ADiGator_OrbitTransferFun_ADiGatorJac.OrbitTransferFun_ADiGatorJac.Gator1Data;
% ADiGator Start Derivative Computations
%User Line: %-----------------------------------------------------------------%
%User Line: % Objective and constraint functions for the orbit-transfer        %
%User Line: % problem.  This function is designed to be used with the NLP     %
%User Line: % solver SNOPT.                                                   %
%User Line: %-----------------------------------------------------------------%
%User Line: %      DO NOT FOR ANY REASON ALTER THE LINE OF CODE BELOW!        %
global psStuff nstates ncontrols npaths CONSTANTS 
%User Line: global
%User Line: %      DO NOT FOR ANY REASON ALTER THE LINE OF CODE ABOVE!        %
%User Line: %-----------------------------------------------------------------%
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
%User Line: % T =control(:,3);
%User Line: %-----------------------------------------------------------------%
%User Line: % The quantity STATEF is the value of the state at the final      %
%User Line: % time, tf, which corresponds to the state at $\tau=1$.           %
%User Line: %-----------------------------------------------------------------%
cada1f1 = size(statePlusEnd.f,1);
stateF.dz0 = statePlusEnd.dz0(Gator1Data.Index18);
stateF.f = statePlusEnd.f(cada1f1,:);
%User Line: stateF = statePlusEnd(end,:);
%User Line: %-----------------------------------------------------------------%
%User Line: % The orbit-raising problem contains one nonlinear boundary       %
%User Line: % condition $\sqrt{mu/r(t_f)-v_\theta(t_f) = 0$.  Because $r(t)$  %
%User Line: % and $v_\theta(t)$ are the first and fourth components of the    %
%User Line: % state, it is necessary to extract stateF(1) and stateF(4) in    %
%User Line: % order to compute this boundary condition function.              %
%User Line: %-----------------------------------------------------------------%
rF.dz0 = stateF.dz0(1);
rF.f = stateF.f(1);
%User Line: rF = stateF(1);
vrF.dz0 = stateF.dz0(3);
vrF.f = stateF.f(3);
%User Line: vrF = stateF(3);
vthetaF.dz0 = stateF.dz0(4);
vthetaF.f = stateF.f(4);
%User Line: vthetaF = stateF(4);
mF.dz0 = stateF.dz0(5);
mF.f = stateF.f(5);
%User Line: mF = stateF(5);
%User Line: % a = T./m;
%User Line: %-----------------------------------------------------------------%
%User Line: % Compute the right-hand side of the differential equations at    %
%User Line: % the N LGR points.  Each component of the right-hand side is     %
%User Line: % stored as a column vector of length N, that is each column has  %
%User Line: % the form                                                        %
%User Line: %                   [ f_i(x_1,u_1,t_1) ]                          %
%User Line: %                   [ f_i(x_2,u_2,t_2) ]                          %
%User Line: %                           .                                     %
%User Line: %                           .                                     %
%User Line: %                           .                                     %
%User Line: %                   [ f_i(x_N,u_N,t_N) ]                          %
%User Line: % where "i" is the right-hand side of the ith component of the    %
%User Line: % vector field f.  It is noted that in MATLABB the calculation of %
%User Line: % the right-hand side is vectorized.                              %
%User Line: %-----------------------------------------------------------------%
rdot.dz0 = vr.dz0; rdot.f = vr.f;
%User Line: rdot = vr;
cada1td1 = zeros(12,1);
cada1td1(Gator1Data.Index19) = vtheta.dz0./r.f(:);
cada1td1(Gator1Data.Index20) = cada1td1(Gator1Data.Index20) + -vtheta.f(:)./r.f(:).^2.*r.dz0;
thetadot.dz0 = cada1td1;
thetadot.f = vtheta.f./r.f;
%User Line: thetadot = vtheta./r;
cada1f1dz0 = 2.*vtheta.f(:).^(2-1).*vtheta.dz0;
cada1f1 = vtheta.f.^2;
cada1td1 = zeros(12,1);
cada1td1(Gator1Data.Index21) = cada1f1dz0./r.f(:);
cada1td1(Gator1Data.Index22) = cada1td1(Gator1Data.Index22) + -cada1f1(:)./r.f(:).^2.*r.dz0;
cada1f2dz0 = cada1td1;
cada1f2 = cada1f1./r.f;
cada1f3dz0 = 2.*r.f(:).^(2-1).*r.dz0;
cada1f3 = r.f.^2;
cada1f4dz0 = -MU./cada1f3(:).^2.*cada1f3dz0;
cada1f4 = MU./cada1f3;
cada1td1 = cada1f2dz0;
cada1td1(Gator1Data.Index23) = cada1td1(Gator1Data.Index23) + -cada1f4dz0;
cada1f5dz0 = cada1td1;
cada1f5 = cada1f2 - cada1f4;
cada1f6dz0 = cos(u1.f(:)).*u1.dz0;
cada1f6 = sin(u1.f);
cada1td1 = zeros(12,1);
cada1td1(Gator1Data.Index24) = cada1f6(:).*T.dz0;
cada1td1(Gator1Data.Index25) = cada1td1(Gator1Data.Index25) + T.f(:).*cada1f6dz0;
cada1f7dz0 = cada1td1;
cada1f7 = T.f.*cada1f6;
cada1tf1 = m.f(Gator1Data.Index26);
cada1td1 = zeros(18,1);
cada1td1(Gator1Data.Index27) = cada1f7dz0./cada1tf1(:);
cada1td1(Gator1Data.Index28) = cada1td1(Gator1Data.Index28) + -cada1f7(:)./m.f(:).^2.*m.dz0;
cada1f8dz0 = cada1td1;
cada1f8 = cada1f7./m.f;
cada1td1 = zeros(30,1);
cada1td1(Gator1Data.Index29) = cada1f5dz0;
cada1td1(Gator1Data.Index30) = cada1td1(Gator1Data.Index30) + cada1f8dz0;
vrdot.dz0 = cada1td1;
vrdot.f = cada1f5 + cada1f8;
%User Line: vrdot = vtheta.^2./r-MU./r.^2+T.*sin(u1)./m;
cada1f1dz0 = -vr.dz0;
cada1f1 = uminus(vr.f);
cada1td1 = zeros(12,1);
cada1td1(Gator1Data.Index31) = vtheta.f(:).*cada1f1dz0;
cada1td1(Gator1Data.Index32) = cada1td1(Gator1Data.Index32) + cada1f1(:).*vtheta.dz0;
cada1f2dz0 = cada1td1;
cada1f2 = cada1f1.*vtheta.f;
cada1tf1 = r.f(Gator1Data.Index33);
cada1td1 = zeros(18,1);
cada1td1(Gator1Data.Index34) = cada1f2dz0./cada1tf1(:);
cada1td1(Gator1Data.Index35) = cada1td1(Gator1Data.Index35) + -cada1f2(:)./r.f(:).^2.*r.dz0;
cada1f3dz0 = cada1td1;
cada1f3 = cada1f2./r.f;
cada1f4dz0 = -sin(u1.f(:)).*u1.dz0;
cada1f4 = cos(u1.f);
cada1td1 = zeros(12,1);
cada1td1(Gator1Data.Index36) = cada1f4(:).*T.dz0;
cada1td1(Gator1Data.Index37) = cada1td1(Gator1Data.Index37) + T.f(:).*cada1f4dz0;
cada1f5dz0 = cada1td1;
cada1f5 = T.f.*cada1f4;
cada1tf1 = m.f(Gator1Data.Index38);
cada1td1 = zeros(18,1);
cada1td1(Gator1Data.Index39) = cada1f5dz0./cada1tf1(:);
cada1td1(Gator1Data.Index40) = cada1td1(Gator1Data.Index40) + -cada1f5(:)./m.f(:).^2.*m.dz0;
cada1f6dz0 = cada1td1;
cada1f6 = cada1f5./m.f;
cada1td1 = zeros(36,1);
cada1td1(Gator1Data.Index41) = cada1f3dz0;
cada1td1(Gator1Data.Index42) = cada1td1(Gator1Data.Index42) + cada1f6dz0;
vthetadot.dz0 = cada1td1;
vthetadot.f = cada1f3 + cada1f6;
%User Line: vthetadot = -vr.*vtheta./r+T.*cos(u1)./m;
cada1f1dz0 = -T.dz0;
cada1f1 = uminus(T.f);
m_dot.dz0 = cada1f1dz0./ve;
m_dot.f = cada1f1/ve;
%User Line: m_dot = -T./ve;
cada1td1 = zeros(90,1);
cada1td1(Gator1Data.Index43) = rdot.dz0;
cada1td1(Gator1Data.Index44) = thetadot.dz0;
cada1td1(Gator1Data.Index45) = vrdot.dz0;
cada1td1(Gator1Data.Index46) = vthetadot.dz0;
cada1td1(Gator1Data.Index47) = m_dot.dz0;
diffeqRHS.dz0 = cada1td1;
diffeqRHS.f = [rdot.f thetadot.f vrdot.f vthetadot.f m_dot.f];
%User Line: diffeqRHS = [rdot, thetadot, vrdot, vthetadot, m_dot];
%User Line: %-----------------------------------------------------------------%
%User Line: % Compute the left-hand side of the defect constraints, recalling %
%User Line: % that the left-hand side is computed using the state at the LGR  %
%User Line: % points PLUS the final point.                                    %
%User Line: %-----------------------------------------------------------------%
cada1td1 = zeros(7,35);
cada1td1(Gator1Data.Index48) = statePlusEnd.dz0;
cada1td1 = D*cada1td1;
cada1td1 = cada1td1(:);
diffeqLHS.dz0 = cada1td1(Gator1Data.Index49);
diffeqLHS.f = D*statePlusEnd.f;
%User Line: diffeqLHS = D*statePlusEnd;
%User Line: %-----------------------------------------------------------------%
%User Line: % Construct the defect constraints at the N LGR points.           %
%User Line: % Remember that the right-hand side needs to be scaled by the     %
%User Line: % factor (tf-t0)/2 because the rate of change of the state is     %
%User Line: % being taken with respect to $\tau\in[-1,+1]$.  Thus, we have    %
%User Line: % $dt/t\dau=(tf-t0)/2$.                                           %
%User Line: %-----------------------------------------------------------------%
cada1td1 = zeros(2,1);
cada1td1(2) = tf.dz0;
cada1td1(1) = cada1td1(1) + -t0.dz0;
cada1f1dz0 = cada1td1;
cada1f1 = tf.f - t0.f;
cada1tempdz0 = cada1f1dz0(Gator1Data.Index50);
cada1tf1 = diffeqRHS.f(Gator1Data.Index51);
cada1td1 = zeros(150,1);
cada1td1(Gator1Data.Index52) = cada1tf1(:).*cada1tempdz0;
cada1td1(Gator1Data.Index53) = cada1td1(Gator1Data.Index53) + cada1f1.*diffeqRHS.dz0;
cada1f2dz0 = cada1td1;
cada1f2 = cada1f1*diffeqRHS.f;
cada1f3dz0 = cada1f2dz0./2;
cada1f3 = cada1f2/2;
cada1td1 = zeros(264,1);
cada1td1(Gator1Data.Index54) = diffeqLHS.dz0;
cada1td1(Gator1Data.Index55) = cada1td1(Gator1Data.Index55) + -cada1f3dz0;
defects.dz0 = cada1td1;
defects.f = diffeqLHS.f - cada1f3;
%User Line: defects = diffeqLHS-(tf-t0)*diffeqRHS/2;
%User Line: %-----------------------------------------------------------------%
%User Line: % Construct the path constraints at the N LGR points.             %
%User Line: %-----------------------------------------------------------------%
%User Line: % paths = [];
%User Line: %-----------------------------------------------------------------%
%User Line: % Reshape the defect contraints into a column vector.             %
%User Line: %-----------------------------------------------------------------%
cada1f1 = N.f*nstates;
defects.dz0 = defects.dz0;
defects.f = reshape(defects.f,cada1f1,1);
%User Line: defects = reshape(defects,N*nstates,1);
%User Line: %-----------------------------------------------------------------%
%User Line: % Reshape the path contraints into a column vector.             %
%User Line: %-----------------------------------------------------------------%
%User Line: % paths = reshape(paths,N*npaths,1);
%User Line: %-----------------------------------------------------------------%
%User Line: % Compute the nonlinear boundary condition.                       %
%User Line: %-----------------------------------------------------------------%
bcs.f =  [];
%User Line: bcs = [];
%User Line: %-----------------------------------------------------------------%
%User Line: % Construct the objective function plus constraint vector.        %
%User Line: %-----------------------------------------------------------------%
cada1f1dz0 = -mF.dz0;
cada1f1 = uminus(mF.f);
cada1td1 = zeros(265,1);
cada1td1(172) = cada1f1dz0;
cada1td1(Gator1Data.Index56) = defects.dz0;
C.dz0 = cada1td1;
C.f = [cada1f1;defects.f];
%User Line: C = [-mF;defects];
C.dz0_size = [31,49];
C.dz0_location = Gator1Data.Index57;
end


function ADiGator_LoadData()
global ADiGator_OrbitTransferFun_ADiGatorJac
ADiGator_OrbitTransferFun_ADiGatorJac = load('OrbitTransferFun_ADiGatorJac.mat');
return
end