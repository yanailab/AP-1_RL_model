% Simulation for parameter scan loops.
% N,M,seed determined from calling script
% no plots (optional)

% AP-1 model - Asymmetric, 3-state model. May 2026
%   present=0, x=0 : absent
%   present=1, x=0 : present but inactive
%   present=1, x=1 : present and active
%
% Stress still depends only on activity pattern x.
% Belief update keeps the same structure, but now distinguishes:
%   active   -> +1
%   inactive -> -1
%   absent   ->  0
%
% Discrete modification variable m (previously continuous "belief")
% Ladder of states [-M:M] with stochastic transitions, rates proportional
% to -Delta S * sgn(activity). only modified for present dimers.

%close all; clear all; clc

%% ----------------------------
% parameters
% ----------------------------

% physical / environmental
%N = 20;     % number of dimers / slots in the cell

%M = 10;     % number of modification states in ladder

T = 20000;  % number of time steps

% flip dynamics
 alpha = 3;    % belief strength
 p0     = 0.05;
 
% occupancy dynamics
p_create = 0.005;  % stress-scaled absent -> inactive prefactor
p_remove = 0.05; % inactive -> absent

% stress sensing
S_delta  = 1; % window size for cellular stress estimate

% belief update value
delta_belief = 2/M;
eta          = delta_belief;

% bookkeeping
record_every = 50;
 %seed = 137;
 %seed  = 142;
 %rng(seed);

%% ----------------------------
% initialization
% ----------------------------

% species identities
%s = repelem(1:N);

% state and belief
present = zeros(N,1); % dimer occupancy
x       = zeros(N,1); % activity
x_new   = zeros(N,1);
m       = zeros(N,1); % stability) v 

% optimal configuration
f     = 0.5;
n1    = floor(f*N);
x_opt = [ones(n1,1); zeros(N-n1,1)];

% recency clock
R = zeros(N,1);

%% ----------------------------
% helper functions
% ----------------------------
clamp           = @(z,lo,hi) max(min(z,hi),lo);
%stress         = @(x,xopt) sqrt(sum(x ~= x_opt));
%stress          = @(x,xopt) sqrt(sum((x-xopt).^2));
% stress functions that don't change scale with N
%stress          = @(x,xopt) mean(x ~= x_opt);
stress          = @(x,xopt) mean((x-xopt).^2);

% % Display flip kinetics
% %-------------------------
% M_axis = [-M:1:M];
% rate   = clamp(p0*exp(-alpha.*M_axis),0,1);
% figure;plot(M_axis,rate,'o-','linewidth',1.5); 
% hold on;
% plot(M_axis,1./(1+exp(-alpha*M_axis)),'sk-','linewidth',1.5);

%% ----------------------------
% storage
% ----------------------------
keep_idx   = 1:record_every:(T+1);
keep_steps = keep_idx - 1;
n_keep     = numel(keep_idx);

X_keep   = nan(n_keep, N);
P_keep   = nan(n_keep, N);
B_keep   = nan(n_keep, N);
S_keep   = nan(n_keep, 1);
S_c_keep = nan(n_keep, 1);
S_hist   = [];

%% ----------------------------
% record initial state (t = 0)
% ----------------------------
S_prev   = stress(x, x_opt);
S_hist   = [S_hist; S_prev];
S_c_prev = mean(S_hist(max(end-S_delta+1,1):end));

k = 1;
X_keep(k,:)   = x';
P_keep(k,:)   = m';
B_keep(k,:)   = m';
S_keep(k)     = S_prev;
S_c_keep(k)   = S_c_prev;


for t = 1:T

    %% A) State Transitions - creation and removal

    present_new = present;   % refrain from flipping back and forth in the same timestep
    x_new       = x;
    
    % stress-dependent creation: production noise decreases as stress is relieved
    %p_create_eff = min(1, p_create * S_c_prev);
    p_create_eff = min(1, p_create * S_c_prev.^2/(S_c_prev.^2+0.5));
    %p_create_eff = min(1, p_create *exp(10*(S_c_prev-0.1))./(1+exp(10*(S_c_prev-0.1))));

    % absent -> present 
    create_idx              = find((present == 0) & (rand(N,1) < p_create_eff));   %indices
    present_new(create_idx) = 1;
    x_new(create_idx)       = 0;        % created with 0 activity
    m(create_idx)           = 0;        % created with 0 stability


    % inactive -> absent  transitions depend on belief
    p_remove_eff            = p_remove ./ (1 + exp(10 * (m./M - (-0.5)))); 
    remove_idx              = find((present == 1)  & (rand(N,1) < p_remove_eff));
    present_new(remove_idx) = 0;
    x_new(remove_idx)       = 0;        % degraded with 0 activity
    m(remove_idx)           = 0;        % degraded with 0 stability


    %% B) activity flips among PRESENT dimers only
    % inactive -> active with fixed probability 
    % active   -> inactive with m-dependent probability

    p_i          = zeros(size(x_new));

    inactive_idx = find((present_new == 1) & (x_new == 0));   
    active_idx   = find((present_new == 1) & (x_new == 1));

    p_i(inactive_idx) = p0;
    %p_i(inactive_idx) = p0 .* exp(alpha .* m(inactive_idx)./M);
    
    %p_i(active_idx)  = p0;
    p_i(active_idx)   = p0 .* exp(-alpha .* m(active_idx)./M);

    p_i   = min(p_i, 1);
    flips = (rand(size(x_new)) < p_i);

    to_activate        = inactive_idx(flips(inactive_idx));
    x_new(to_activate) = 1;

    to_inactivate        = active_idx(flips(active_idx));
    x_new(to_inactivate) = 0;

    % hard safety constraint: absent slots cannot be active
    x_new(present_new == 0) = 0;

    % debug check
    if any(x_new(present_new == 0) ~= 0)
        error('Inconsistent state: active entry in an absent slot');
    end

    %% C) stress
    % stress depends only on activity pattern
    S_new  = stress(x_new, x_opt);
    S_hist = [S_hist; S_new];

    idx0    = max(numel(S_hist)-S_delta+1,1);
    S_c_new = mean(S_hist(idx0:end));

    dS = S_c_new - S_c_prev;

    %% D) belief update
    % 3-state sign:
    %   active   -> +1
    %   inactive -> -1
    %   absent   ->  0

    state_sign = zeros(N,1);
    state_sign((present_new == 1) & (x_new == 0)) = -1;
    state_sign((present_new == 1) & (x_new == 1)) = +1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- discrete ladder update ---


dir = -sign(dS) .* state_sign;   % desired step direction: -1, 0, or +1
                                 % decrease in stress + activity - increase
                                 % 

% Fixed jump probability, only sign feedback. Too crude
jump = rand(N,1) < eta;
m = m + dir .* jump;
m = max(min(m, M), -M);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Simple deterministic update
    %belief = clamp(belief - delta_belief * sign(dS) .* state_sign, -1, 1);

    %% E) advance
    present  = present_new;
    x        = x_new .* present_new;   % extra hard enforcement
    S_prev   = S_new;
    S_c_prev = S_c_new;

    if ismember(t+1, keep_idx)
        k = k + 1;
        X_keep(k,:)   = x';
        P_keep(k,:)   = present';
        B_keep(k,:)   = m';
        S_keep(k)     = S_prev;
        %S_c_keep(k)   = S_c_prev;
        S_c_keep(k)   = dS;
    end
end

%% ----------------------------
% package results
% ----------------------------
res.params = struct( ...
    'N',N,'T',T, ...
    'alpha',alpha,'p0',p0, ...
    'p_create',p_create,'p_remove',p_remove, ...
    'S_delta',S_delta, ...
    'delta_belief',delta_belief, ...
    'record_every',record_every);
%    'record_every',record_every,'seed',seed);    % Record also seed


res.x_opt = x_opt;
res.X = X_keep;
res.present = P_keep;
res.belief = B_keep;
res.stress = S_keep;
res.stress_cellular = S_c_keep;
res.keep_steps = keep_steps;

%% ----------------------------
% plots
% ----------------------------
time = keep_steps;

% % stress
% %%%%%%%%
% figure; hold on;
% plot(time, res.stress./f, 'LineWidth',1.5,'Color',[0.4 0.4 0.5]);hold on;
% %plot(time, res.stress_cellular, 'LineWidth',1.5);
% xlabel('Time step'); ylabel('stress S');
% box on; ylim([0 1.1]);
% ax = gca;                     
% ax.FontSize = 14; 

% % % belief
% %%%%%%%%%%
%  figure; hold on;
% % Plot all smoothed traj
% b=B_keep(:,1:n1);
% Bp = smoothdata(b, 1, 'movmean', 3);
% b=B_keep(:,n1+1:N);
% Bm = smoothdata(b, 1, 'movmean', 3);
% plot(time,Bp, 'LineWidth',1,'Color',[0 0.4470 0.7410]); hold on
% plot(time,Bm, 'LineWidth',1,'Color',[0.85 0.325 0.098]);
% % Plot average over +- groups
 mBp = mean(B_keep(:,1:n1)');
 mBm = mean(B_keep(:,n1+1:N)');
% plot(time,mBp,'k-','linewidth',2); hold on
% plot(time,mBm, 'k-','linewidth',2);
% ylim([-M M].*1.2); box on;
% xlabel('Time step'); ylabel('Stability');
% ax = gca;                     
% ax.FontSize = 14; 
% % % Plot bare traj
% % plot(time,B_keep(:,1:N/2), 'LineWidth',1,'Color',[0 0.4470 0.7410]); hold on
% % plot(time,B_keep(:,N/2+1:N), 'LineWidth',1,'Color',[0.85 0.325 0.098]);
% % Plot average +- trajectories

% % Fraction Dimers
% %%%%%%%%%%%%%%%%%
% % fraction present
 mP = mean(P_keep,2);
% figure;
% plot(time, mP, 'LineWidth',1.5);
% ylim([0 1.1]);
% xlabel('time step'); ylabel('fraction present');
% %title('fraction present (all dimers)');
% %mean active among present only
 frac_active_given_present = sum(X_keep,2) ./ max(sum(P_keep,2),1);
% hold on;
% plot(time, frac_active_given_present, 'LineWidth',1.5,'Color','k');
% ylim([0 1.1]);
% xlabel('time step'); ylabel('fraction active among present');
% legend('present','active|present')
% %title('activity conditional on presence');




