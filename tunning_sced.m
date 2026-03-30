format short
clc;clear;close;

cfg = struct();

% Core
cfg.n = 300;
cfg.p = 10;
cfg.N = 200;
cfg.seed = 1;

% Scenario
cfg.rho = 0.25;                 % correlation used in this run
cfg.tau = 1.60;                 % censoring upper bound

% True Beta
cfg.Beta = 0.5*[-1.0;1.0;0;0;0;-1.0;zeros(cfg.p-6,1)];

% Tuning
cfg.theta = (0.575:-0.002:0.535)*1.25e4;   % or [] to use defaults by case
cfg.delta = 0.067;                           % or vector, or [] for default
cfg.Q = 1;                                   % number of initial-value replicates
cfg.tol = 1e-5;
cfg.tkini = 4;   % step size in ini_beta
cfg.tkmain = 46; % step size hlselo_lqa


out = sced_lqa(cfg);