%% Clean up the workspace and console

clear; clc; close all;

%% Load the variables which appear in the example

load('variables.mat');

%% Form the realzation of the plant

G = ss([-0.5 * eye(6) - 0.1 * F_b -0.1 * F_b; zeros(6) eye(6)],...
       [zeros(6); eye(6)], [eye(6) eye(6)], zeros(6));

%% Check the minimality of the plant's realization

nc_modes_G = tzero(G.a, G.b, [], []); % Empty, so the realization is controllable
no_modes_G = tzero(G.a, [], G.c, []); % Empty, so the realization is observable

%% Form the realization of the centralized controller

K_ctr = ss([A11_K A12_K; A21_K A22_K], [B1_K; B2_K],...
           [eye(6) zeros(6)], zeros(6));
K_ctr_0 = -K_ctr.c*(K_ctr.a\K_ctr.b) + K_ctr.d; % Full matrix, no sparsity pattern

%% Validate the stability of the closed-loop system with centralized control laws

CL_sys_ctr = feedback(G, K_ctr, 1);
if max(real(eig(CL_sys_ctr.a))) >= 0
    msg_ctr_stab = 'The centralized control laws are not stabilizing.';
else
    msg_ctr_stab = 'The centralized control laws are stabilizing.';
end
disp(msg_ctr_stab);

%% Form the realization of the stable SRTR pair

SRTR = ss(A22_K + K * A12_K,...
       [-(A22_K * K + K * A12_K * K - K * A11_K - A21_K), K * B1_K + B2_K],...
          A12_K, [A11_K - A12_K * K B1_K]);

%% Check the satisfaction of the model-matching conditions

row_SRTR = cell(1,6);
cond_thm = zeros(6);
model_matching_fail = 0;

for i = 1:6

    row_SRTR{i} = ss2ss(SRTR(i,:), Q_i{i});
    null_cols = 1 + mod((i:i+3),6);
    cond_thm(i,1) = norm(row_SRTR{i}.d(null_cols)) <  1e-9;
    cond_thm(i,2) = norm(row_SRTR{i}.d(6 + null_cols)) <  1e-9;
    cond_thm(i,3) = norm(row_SRTR{i}.b(6,null_cols)) <  1e-9;
    cond_thm(i,4) = norm(row_SRTR{i}.b(6,6 + null_cols)) <  1e-9;
    cond_thm(i,5) = norm(row_SRTR{i}.a(6,1:5)) <  1e-9;
    cond_thm(i,6) = real(row_SRTR{i}.a(6,6)) <  0;

    if prod(cond_thm(i,:)) < 1
        model_matching_fail = 1;
        break
    end

end

if model_matching_fail > 0
    msg_procedure = 'The procedure has not been succesfully carried out.';
else
    msg_procedure = 'The procedure has been succesfully carried out.';
end
disp(msg_procedure);

%% Form minimal realizations for each distributed subcontroller

subctl_dist = cell(1,6);
K_dist = ss([],[],[],[]); % Realization of the global, agregated controller
for i = 1:6
    row_SRTR{i} = ss(row_SRTR{i}.a(6,6),round(row_SRTR{i}.b(6,:),9),...
                     row_SRTR{i}.c(:,6),round(row_SRTR{i}.d(:,:),9));
    subctl_dist{i} = ss(tf(1,[1 0]) * row_SRTR{i}, 'min');
    K_dist = ss(blkdiag(K_dist.a, subctl_dist{i}.a),...
                vertcat(K_dist.b, subctl_dist{i}.b),...
                blkdiag(K_dist.c, subctl_dist{i}.c),...
                vertcat(K_dist.d, subctl_dist{i}.d));
end

nc_modes_K_dist = tzero(K_dist.a, K_dist.b, [], []); % Empty, so the realization is controllable
no_modes_K_dist = tzero(K_dist.a, [], K_dist.c, []); % Empty, so the realization is observable

%% Validate the stability of the closed-loop system with distributed control laws

CL_sys_dist = feedback([eye(6); G], K_dist, 1);
if max(real(eig(CL_sys_dist.a))) >= 0
    msg_dist_stab = 'The distributed control laws are not stabilizing.';
else
    msg_dist_stab = 'The distributed control laws are stabilizing.';
end
disp(msg_dist_stab);

%% Compute the nonzero entries of the SRTR pair

W_local = cell(1,6);
V_local = cell(1,6);
for i = 1:6
    W_local{i} = tf(row_SRTR{i}(i));
    V_local{i} = tf(row_SRTR{i}(i+6));
end

W_prev  = cell(1,6);
V_prev  = cell(1,6);

W_prev{1} = tf(row_SRTR{1}(6));
V_prev{1} = tf(row_SRTR{1}(12));
for i = 2:6
    W_prev{i} = tf(row_SRTR{i}(i-1));
    V_prev{i} = tf(row_SRTR{i}(i+5));
end