clear all;
close all;
addpath('C:\Users\damsa\OneDrive\Desktop\casadi-matlabR2014b-v3.3.0');

import casadi.*

%%

delta_t = 0.1; % sampling time [s]
mass = 1;
N = 100; % prediction horizon

theta_max = pi/6; theta_min = -theta_max;
phi_max = pi/6;   phi_min = -phi_max;
T_max = 12;      T_min = 8;
psi_max = pi/4; psi_min = -psi_max;

x1 = SX.sym('x1'); x2 = SX.sym('x2'); 
x3 = SX.sym('x3'); x4 = SX.sym('x4'); 
x5 = SX.sym('x5'); x6 = SX.sym('x6'); 
states = [x1;
          x2;
          x3;
          x4;
          x5;
          x6]; 

n_states = length(states);

theta = SX.sym('theta'); phi = SX.sym('phi'); psi = SX.sym('psi'); T = SX.sym('T');

controls = [theta;
            phi;
            psi;
            T]; 
n_controls = length(controls);

% rhs = getDerivative(delta_t, mass, x1, x2, x3, x4, x5, x6, theta, phi, T);
rhs = getDerivative_fsm(delta_t, mass, x1, x2, x3, x4, x5, x6, theta, phi, psi, T);% system r.h.s

f = Function('f',{states, controls},{rhs}); % nonlinear mapping function f(x,u)
U = SX.sym('U', n_controls, N); % Decision variables (controls)
Param = SX.sym('P', n_states+ n_states*N);
% parameters (which include the initial and the reference state of the robot)

X = SX.sym('X',n_states,(N + 1));
% A Matrix that represents the states over the optimization problem.

Q = zeros(6,6); 
Q(1,1) = 1; Q(2,2) = 1; Q(3,3) = 1; 
Q(4,4) = 1; Q(5,5) = 1; Q(6,6) = 1; % weighing matrices (states)
R = zeros(4,4); R(1,1) = 1/theta_max; R(2,2) = 1/phi_max; R(3,3) = 1/psi_max;  R(4,4) = 1/(T_max+100); % weighing matrices (controls)

obj = 0; % Objective function
g = [];  % constraints vector

st  = X(:,1);% initial state 
g = [g; st - Param(1:6)]; % initial condition constraints

% compute solution symbolically
for k = 1:N
    st = X(:,k);  con = U(:,k);
    state_vec = st - Param(7+6*(k-1):12+6*(k-1));
    obj = obj  + state_vec'*Q*state_vec + con'*R*con; % calculate obj
    st_next = X(:,k+1);
    f_value = f(st, con);
    st_next_euler = st + (delta_t*f_value);
    g = [g;st_next-st_next_euler]; % compute constraints
end



% make the decision variable one column  vector
OPT_variables = [reshape(X,6*(N + 1),1);reshape(U,4*N,1)];

nlp_prob = struct('f', obj, 'x', OPT_variables, 'g', g, 'p', Param);

opts = struct;
opts.ipopt.max_iter = 2000;
opts.ipopt.print_level =0;%0,3
opts.print_time = 0;
opts.ipopt.acceptable_tol =1e-8;
opts.ipopt.acceptable_obj_change_tol = 1e-6;

solver = nlpsol('solver', 'ipopt', nlp_prob,opts);

args = struct;
args.lbg(1:6*(N+1)) = 0;
args.ubg(1:6*(N+1)) = 0;


args.lbx(1:6:6*(N+1),1) = -1; %state x lower bound
args.ubx(1:6:6*(N+1),1) = 40; %state x upper bound
args.lbx(2:6:6*(N+1),1) = -2; %state x_dot lower bound
args.ubx(2:6:6*(N+1),1) = 4; %state x_dot upper bound
args.lbx(3:6:6*(N+1),1) = -1; %state y lower bound
args.ubx(3:6:6*(N+1),1) = 40; %state y upper bound
args.lbx(4:6:6*(N+1),1) = -2; %state y_dot lower bound
args.ubx(4:6:6*(N+1),1) = 4; %state y_dot upper bound
args.lbx(5:6:6*(N+1),1) = -11; %state z lower bound
args.ubx(5:6:6*(N+1),1) = 20; %state z upper bound
args.lbx(6:6:6*(N+1),1) = -5; %state z_dot lower bound
args.ubx(6:6:6*(N+1),1) = 5; %state z_dot upper bound

args.lbx(6*(N+1)+1:4:6*(N+1)+4*N,1) = theta_min; %v lower bound
args.ubx(6*(N+1)+1:4:6*(N+1)+4*N,1) = theta_max; %v upper bound
args.lbx(6*(N+1)+2:4:6*(N+1)+4*N,1) = phi_min; %omega lower bound
args.ubx(6*(N+1)+2:4:6*(N+1)+4*N,1) = phi_max; %omega upper bound
args.lbx(6*(N+1)+3:4:6*(N+1)+4*N,1) = psi_min; %v lower bound
args.ubx(6*(N+1)+3:4:6*(N+1)+4*N,1) = psi_max; %v upper bound
args.lbx(6*(N+1)+4:4:6*(N+1)+4*N,1) = T_min; %v lower bound
args.ubx(6*(N+1)+4:4:6*(N+1)+4*N,1) = T_max; %v upper bound
%----------------------------------------------
% ALL OF THE ABOVE IS JUST A PROBLEM SET UP


% THE SIMULATION LOOP SHOULD START FROM HERE
%-------------------------------------------
t0 = 0;
x0 = [0.0; 0.0; 0.0; 0.0; 0; 0.0];    % initial condition.

xx(:,1) = x0; % xx contains the history of states
t(1) = t0;

u0 = zeros(N,4);        % two control inputs for each robot
X0 = repmat(x0,1,N+1)'; % initialization of the states decision variables

sim_tim = 15; % total sampling times

% Start MPC
mpciter = 0;
xx1 = [];
u_cl=[];

% the main simulaton loop... it works as long as the error is greater
% than 10^-6 and the number of mpc steps is less than its maximum
% value.
x_vel_drogue = 1;
y_vel_drogue = 1;
z_vel_drogue = 1;
init_drogue = [5; x_vel_drogue; 5; y_vel_drogue; -5; z_vel_drogue];
while(mpciter < sim_tim / delta_t) 
    tic
    current_time = mpciter*delta_t;
    args.p(1:6)  = x0;
    for k=1:N
        t_predict = current_time + (k-1)*delta_t;
        drogue_state = init_drogue + t_predict*drogue_traj(x_vel_drogue, y_vel_drogue, z_vel_drogue);
        args.p(7+6*(k-1):12+6*(k-1)) = drogue_state;
        drogue_state_traj(:,mpciter+k) = drogue_state;
    end
    
%     init_drogue = drogue_state;
    similarity = sqrt((drogue_state_traj(:,mpciter+1) - x0)'*(drogue_state_traj(:,mpciter+1) - x0));

    if similarity < 1e-2
        break;
    end
    drogue_states_evo(:,mpciter+1) = drogue_state_traj(:,mpciter+1);
%     disp("distance = ", similarity);
    args.x0  = [reshape(X0',6*(N+1),1);reshape(u0',4*N,1)];
    sol = solver('x0', args.x0, 'lbx', args.lbx, 'ubx', args.ubx,...
        'lbg', args.lbg, 'ubg', args.ubg,'p',args.p);
    u = reshape(full(sol.x(6*(N+1)+1:end))',4,N)'; % get controls only from the solution
%     xx1(:,1:6,mpciter+1)= reshape(full(sol.x(1:6*(N+1)))',6,N+1)'; % get solution TRAJECTORY
    u_cl= [u_cl ; u(1,:)];

    t(mpciter+1) = t0;
    % Apply the control and shift the solution
    [t0, x0, u0] = shift(delta_t, t0, x0, u,f);
    xx(:,mpciter+1) = x0;
    X0 = reshape(full(sol.x(1:6*(N+1)))',6,N+1)'; % get solution TRAJECTORY
    % Shift trajectory to initialize the next step
    X0 = [X0(2:end,:); X0(end,:)];
%     xs = [5 + delta_t*10; 0.0; 5; 0.0; -5; 0.0];
    mpciter;
    mpciter = mpciter + 1;
    toc
end


%%
figure(1)
subplot(611)
plot(t,xx(1,1:end),'b','linewidth',1.5); hold on;
plot(t,drogue_states_evo(1,1:end),'--r','linewidth',1.5);
ylabel('x (m)')
grid on
subplot(612)
plot(t,xx(2,1:end),'b','linewidth',1.5); hold on;
plot(t,drogue_states_evo(2,1:end),'--r','linewidth',1.5);
ylabel('xdot(m/s)')
grid on
subplot(613)
plot(t,xx(3,1:end),'b','linewidth',1.5); hold on;
plot(t,drogue_states_evo(3,1:end),'--r','linewidth',1.5);
ylabel('y(m)')
grid on
subplot(614)
plot(t,xx(4,1:end),'b','linewidth',1.5); hold on;
plot(t,drogue_states_evo(4,1:end),'--r','linewidth',1.5);
ylabel('ydot(m/s)')
grid on
subplot(615)
plot(t,xx(5,1:end),'b','linewidth',1.5); hold on;
plot(t,drogue_states_evo(5,1:end),'--r','linewidth',1.5);
ylabel('z (m)')
grid on
subplot(616)
plot(t,xx(6,1:end),'b','linewidth',1.5); hold on;
plot(t,drogue_states_evo(6,1:end),'--r','linewidth',1.5);
xlabel('time (seconds)')
ylabel('zdot(m/s)')
grid on

figure(2)
subplot(411)
plot(t,u_cl(:, 1),'k','linewidth',1.5);
ylabel('theta_c (rad)')
grid on
subplot(412)
plot(t,u_cl(:, 2),'k','linewidth',1.5); 
ylabel('phi_c (rad)')
grid on
subplot(413)
plot(t,u_cl(:, 3),'k','linewidth',1.5); 
ylabel('psi_c (rad)')
grid on
subplot(414)
plot(t,u_cl(:, 4),'k','linewidth',1.5); 
xlabel('time (seconds)')
ylabel('Thr (N)')
grid on

figure(2)
subplot(411)
plot(t,u_cl(:, 1),'k','linewidth',1.5);
ylabel('theta_c (rad)')
grid on
subplot(412)
plot(t,u_cl(:, 2),'k','linewidth',1.5); 
ylabel('phi_c (rad)')
grid on
subplot(413)
plot(t,u_cl(:, 3),'k','linewidth',1.5); 
ylabel('psi_c (rad)')
grid on
subplot(414)
plot(t,u_cl(:, 4),'k','linewidth',1.5); 
xlabel('time (seconds)')
ylabel('Thr (N)')
grid on

function [t0, x0, u0] = shift(delta_t, t0, x0, u,f)
    % add noise to the control actions before applying it
    con_cov = diag([0.005 deg2rad(2)]).^2;
    con = u(1,:)'; 
    st = x0;
    
    f_value = f(st,con);   
    st = st + (delta_t*f_value);
    
    x0 = full(st);
    t0 = t0 + delta_t;
    
    u0 = [u(2:size(u,1),:);u(size(u,1),:)]; % shift the control action 
end

function rhs = drogue_traj(x2, x4, x6) 
    r1 = x2;
    r2 = 0;
    r3 = x4;
    r4 = 0;
    r5 = x6;
    r6 = 0;
    rhs = [r1; r2; r3; r4; r5; r6];
end




