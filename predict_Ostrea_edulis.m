

function [prdData, info] = predict_Ostrea_edulis(par, data, auxData)
  
  % unpack par, data, auxData
  cPar = parscomp_st(par); vars_pull(par); 
  vars_pull(cPar);  vars_pull(data);  vars_pull(auxData);
  
  if f_Hutc>1 || f_Hutc<0.2 ||...
     f_LN>1 ||f_LN<0.2 ||...
     f_Robe0>1 ||f_Robe0<0.2 ||...
     f_Laba>1 ||f_Laba<0.2 ||...
     f_Robe1>0.8 ||f_Robe1<0.2 ||...
     f_Robe2>0.8 || f_Robe2<0.2 ||...
     f_Saur>0.9 ||f_Saur<0.2 ||...
     f_Wils>0.9 ||f_Wils<0.2 ||...
     ERrobe<0 ||...
     L_0Saur<0 ||...
     del_Mu<0 ||...
     del_Mb<0 ||...
     L_0Saur <0 ||...
     E_Hb > E_Hr || E_Hr > E_Hs || E_Hs > E_Hj || E_Hj > E_Hp
        info = 0; prdData = {}; return;
  end 
    
  % compute temperature correction factors
  pars_T = [T_A, T_L, T_H, T_AL, T_AH];  
  %
  TC_am = tempcorr(temp.am,T_ref, pars_T);
  TC_ab = tempcorr(temp.ab,T_ref, pars_T);
  TC_ar = tempcorr(temp.ar, T_ref, pars_T);
  TC_Robe_15 = tempcorr(temp.tL_Robe_15,T_ref,pars_T);
  TC_Robe_20 = tempcorr(temp.tL_Robe_20,T_ref,pars_T);
  TC_Robe_25 = tempcorr(temp.tL_Robe_25,T_ref,pars_T);
  TC_Robe_30 = tempcorr(temp.tL_Robe_30,T_ref,pars_T);
  TC_Laba = tempcorr(temp.tL_Laba,T_ref,pars_T);
  TC_Saur = tempcorr(temp.tL_Saur,T_ref,pars_T);
%   TC_Mann_12 = tempcorr(temp.tWd_Mann_12,T_ref,pars_T);
%   TC_Mann_15 = tempcorr(temp.tWd_Mann_15,T_ref,pars_T);
%   TC_Mann_18 = tempcorr(temp.tWd_Mann_18,T_ref,pars_T);
 %TC_Mann_21 = tempcorr(temp.tWd_Mann_21,T_ref,pars_T);
  TC_LN = tempcorr(temp.LN, T_ref,pars_T);
  TC_Hutc = tempcorr(C2K(data.WdJO(:,1)), T_ref, pars_T);
  TC_Walne = tempcorr(temp.LN, T_ref, pars_T);
  
  %% zero-variate data

  % life cycle
  pars_ts = [g k l_T v_Hb v_Hs v_Hj v_Hp];
  [tau_s, tau_j, tau_p, tau_b, l_s, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_ts(pars_ts, f);  
  if info ~= 1 % numerical procedure failed
    info = 0; prdData = []; return
  end
  pars_tr = [g k l_T v_Hb v_Hr]; % v_Hr replaces v_Hp in call to get_tp
  [tau_r, ~, l_r] = get_tp(pars_tr, f); % scaled age, length at release
 
  % initial
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  E_0 = p_Am * initial_scaled_reserve(f, pars_UE0); % J, initial energy in egg
  Ww_0 = E_0 * w_E/ mu_E/ d_E;        % g, initial wet weight
  
  % birth
  L_b = L_m * l_b;                   % cm, structural length at birth at f
  Lw_b = L_b/ del_Mb;                 % cm, shell length at birth
  aT_b = tau_b / k_M / TC_ab;         % d,  age at birth at f and T
  
  % release
  L_r = L_m * l_r;                   % cm, structural length at release at f
  Lw_r = L_r/ del_Mb;                 % cm, shell length at settlement at f
  Wd_r = L_r^3 * (1 + f * ome) * d_V; % g, dry weight at settlement
  aT_r  = tau_r/ k_M / TC_ar;         % d,  age at birth at f and T
  
  % start acceleration
  L_s = L_m * l_s;                  % cm, structural length at start acceleration
  Lw_s = L_s/ del_Mb;               % cm, shell length at start acceleration
  Wd_s = L_s^3 * (1 + f * ome) * d_V; % g, dry weight at settlement

  % metamorphosis(= end acceleration)
  s_M = l_j/ l_s;                     % -, acceleration factor

  % puberty 
  L_p = L_m * l_p;                  % cm, structural length at puberty at f
  Ww_p = L_p^3 * (1 + f * ome);     % g, wet weight at puberty 
  
  % ultimate
  [tau_s, tau_j, tau_p, tau_b, l_s, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_ts(pars_ts, f_Saur); 
  L_i  = L_m * l_i;                 % cm, ultimate structural length at f
  Wd_i = L_i^3 * d_V * (1 + f_Saur * w); % g,  ultimate dry weight
  
  % reproduction  
  pars_R = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hs; U_Hj; U_Hp]; % compose parameter vector at T

  % life span
    [tau_s, tau_j, tau_p, tau_b, l_s, l_j, l_p, l_b, l_i, rho_j, rho_B, info] = get_ts(pars_ts, f);
  pars_tm = [g; l_T; h_a/ k_M^2; s_G];  % compose parameter vector at T_ref
  t_m = get_tm_s(pars_tm, f, l_b);      % -, scaled mean life span at T_ref
  aT_m = t_m/ k_M/ TC_am;               % d, mean life span at T
   
  % max clearance rate
  ECR_max   = (1 / d_V / (1 + f * ome)).^(2/3) * F_m * 0.9 * 1 * s_M/ 24; % l/h, clearance rate

  % pack to output 
  prdData.am    = aT_m;
  prdData.ab    = aT_b;
  prdData.ar    = aT_r;

  prdData.Lb    = Lw_b;
  prdData.Lr    = Lw_r;
  prdData.Ls    = Lw_s;
  
 % prdData.Ww0   = Ww_0;
  prdData.Wdr   = Wd_r;
  prdData.Wds   = Wd_s;
  prdData.Wwp   = Ww_p;
  prdData.Wdi   = Wd_i;
  prdData.CR_max =ECR_max;
  
  %% uni-variate data  
  
% LARVAL GROWTH 
  % tL from Robe2017
  [tau_s, tau_j, tau_p, tau_b, l_s, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_ts(pars_ts, f_Robe0);
  [tau_r, ~, l_r] = get_tp(pars_tr, f_Robe0); % scaled age, length at release
  
%   %15°C       
  kT_M = k_M * TC_Robe_15; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M; 
  tT_s = (tau_s - tau_r)/ kT_M; L_r = 0.0180 *del_Mb; L_s = L_m * l_s;  L_ib = f_Robe0 * L_m; L_i = L_m * l_i;  
  L_bs = L_ib - (L_ib - L_r) * exp( - rT_B * tL_Robe_15(tL_Robe_15(:,1) < tT_s,1));  
  L_sj = L_s * exp((tL_Robe_15(tL_Robe_15(:,1) >= tT_s) - tT_s) * rT_j/ 3);
  EtL_Robe_15 = [L_bs; L_sj]./ del_Mb;   % cm, shell length
  %20°C
  kT_M = k_M * TC_Robe_20; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M; 
  tT_s = (tau_s - tau_r)/ kT_M; L_r = 0.0180 *del_Mb; L_s = L_m * l_s;  L_ib = f_Robe0 * L_m; L_i = L_m * l_i;  
  L_bs = L_ib - (L_ib - L_r) * exp( - rT_B * tL_Robe_20(tL_Robe_20(:,1) < tT_s,1));  
  L_sj = L_s * exp((tL_Robe_20(tL_Robe_20(:,1) >= tT_s,1) - tT_s) * rT_j/ 3);
  EtL_Robe_20 = [L_bs; L_sj]./ del_Mb;   % cm, shell length
%   %25°C
  kT_M = k_M * TC_Robe_25; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M; 
  tT_s = (tau_s - tau_r)/ kT_M; L_r = 0.0180 *del_Mb; L_s = L_m * l_s;  L_ib = f_Robe0 * L_m; L_i = L_m * l_i;  
  L_bs = L_ib - (L_ib - L_r) * exp( - rT_B * tL_Robe_25(tL_Robe_25(:,1) < tT_s,1));  % cm, struc length
  L_sj = L_s * exp((tL_Robe_25(tL_Robe_25(:,1) >= tT_s,1) - tT_s) * rT_j/ 3);
  EtL_Robe_25 = [L_bs; L_sj]./ del_Mb;   % cm, shell length
%   %30 °C  
  kT_M = k_M * TC_Robe_30; rT_j = rho_j * kT_M; rT_B = rho_B * kT_M; 
  tT_s = (tau_s - tau_r)/ kT_M; L_r = 0.0180 *del_Mb; L_s = L_m * l_s;  L_ib = f_Robe0 * L_m; L_i = L_m * l_i;  
  L_bs = L_ib - (L_ib - L_r) * exp( - rT_B * tL_Robe_30((tL_Robe_30(:,1) < tT_s),1));  % cm, struc length
  L_sj = L_s * exp((tL_Robe_30(tL_Robe_30(:,1) >= tT_s,1) - tT_s) * rT_j/ 3);
  EtL_Robe_30 = [L_bs; L_sj]./ del_Mb;   % cm, shell length

  
% POST-LARVAL GROWTH

  % tWd from Laba1999
  [tau_s, tau_j, tau_p, tau_b, l_s, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_ts(pars_ts, f_Laba);
  [tau_r, ~, l_r] = get_tp(pars_tr, f_Laba); % scaled age, length at release
  L_r = l_r * L_m; L_s = l_s * L_m; L_j = l_j * L_m; 
  eLRH0 = [f_Laba, L_r, 0, E_Hr];
  [t eLRH]=ode45(@dget_eLRH, tWd_Laba(:,1), eLRH0, [], f_Laba, TC_Laba, v, L_s, L_j, L_T, L_m, g, kap, E_m, k_J, E_Hp);
    EtL_Laba = eLRH(:,2)./del_Mb;
    EtWd_Laba = eLRH(:,2).^3 * (1 + f_Laba * ome) * d_V; % g, tissue dry weight
    
% ADULT GROWTH
  % tL data from Saurel 2003
  [tau_s, tau_j, tau_p, tau_b, l_s, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_ts(pars_ts, f_Saur);
  kT_M = k_M * TC_Saur; rT_B = rho_B * kT_M; L_i = L_m * l_i;   
  L  = L_i - (L_i - L_0Saur) * exp( - rT_B * tL_Saur(:,1)); % cm, struc length 
  EtL_Saur = L./ del_Mu;   % cm, shell length   
  
  % Mann1975
  [tau_s, tau_j, tau_p, tau_b, l_s, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_ts(pars_ts, 0.5);
  TC12 = tempcorr(C2K(12),T_ref,pars_T); 
%   TC15 = tempcorr(C2K(15),T_ref,pars_T); 
%   TC18 = tempcorr(C2K(18),T_ref,pars_T); 
%   TC21 = tempcorr(C2K(21),T_ref,pars_T); 

%   %12°C
  [tau_s, tau_j, tau_p, tau_b, l_s, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_ts(pars_ts, f_Mann);
 L_i = L_m * l_i; s_M = l_j/ l_s; r_B = rho_B * k_M;    
   options = odeset('AbsTol',1e-8, 'RelTol',1e-8);   
  LR_0 = [L_p 0];
  t_R = spawn.tWd_Mann_12+t_Mann;   % d, day at spawning
  tt = tWd_Mann_12 + t_Mann;        % d, since puberty
  t = [0; tt(tt(:,1) < t_R,1); t_R]; 
  %Untill reproduction
  [t_0R, E_R0R] = ode45(@dget_LR, t, LR_0, options, [], E_Hp, TC12, f_Mann, r_B, L_i, k_J, E_m, v *s_M, kap); % find state variables at start experiment
  E_R0R(1,:) = []; E_R0R(end,:) = [];
  %during brooding
  t = [t_R; tt(tt(:,1) >= t_R,1)];
  LR_0 = E_R0R(end,:);
  LR_0(2) = 0;
  [t_Ri, E_RRi] = ode45(@dget_LR, t, LR_0, options, [], E_Hp, TC12, f_Mann, r_B, L_i, k_J, E_m, v *s_M, kap); % find state variables at start experiment
  E_RRi(1,:) = [];
  LR = [E_R0R; E_RRi];
  EtWd_Mann_12 = LR(:,1).^3 * (1 + f_Mann * ome) * d_V + LR(:,2)/ mu_E * w_E; % g, dry weight
  % cluch weight
  pars_UE0 = [V_Hb; g; k_J; k_M; v]; % compose parameter vector
  E_0 = p_Am * initial_scaled_reserve(f, pars_UE0); % J, initial energy in egg
  RR = E_R0R(end,2)/E_0;    % reproduction rate
  EWd_clutch12 = RR * (Wd_r*3);    % weight loss during reproduction

% 
% %   %15°C
%   [tau_s, tau_j, tau_p, tau_b, l_s, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_ts(pars_ts, f_Mann);
%  L_i = L_m * l_i; s_M = l_j/ l_s; r_B = rho_B * k_M;    
%    options = odeset('AbsTol',1e-8, 'RelTol',1e-8);   
%   LR_0 = [L_p 0];
%   t_R = spawn.tWd_Mann_15+t_Mann;   % d, day at spawning
%   tt = tWd_Mann_15 + t_Mann;        % d, since puberty
%   t = [0; tt(tt(:,1) < t_R,1); t_R]; 
%   %Untill reproduction
%   [t_0R, E_R0R] = ode45(@dget_LR, t, LR_0, options, [], E_Hp, TC12, f_Mann, r_B, L_i, k_J, E_m, v *s_M, kap); % find state variables at start experiment
%   E_R0R(1,:) = []; E_R0R(end,:) = [];
%   %during brooding
%   t = [t_R; tt(tt(:,1) >= t_R,1)];
%   LR_0 = E_R0R(end,:);
%   LR_0(2) = 0;
%   [t_Ri, E_RRi] = ode45(@dget_LR, t, LR_0, options, [], E_Hp, TC15, f_Mann, r_B, L_i, k_J, E_m, v *s_M, kap); % find state variables at start experiment
%   E_RRi(1,:) = [];
%   LR = [E_R0R; E_RRi];
%   EtWd_Mann_15 = LR(:,1).^3 * (1 + f_Mann * ome) * d_V + LR(:,2)/ mu_E * w_E; % g, dry weight
%   
% %   %18°C
%   [tau_s, tau_j, tau_p, tau_b, l_s, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_ts(pars_ts, f_Mann);
%  L_i = L_m * l_i; s_M = l_j/ l_s; r_B = rho_B * k_M;    
%    options = odeset('AbsTol',1e-8, 'RelTol',1e-8);   
%   LR_0 = [L_p 0];
%   t_R = spawn.tWd_Mann_18+t_Mann;   % d, day at spawning
%   tt = tWd_Mann_18 + t_Mann;        % d, since puberty
%   t = [0; tt(tt(:,1) < t_R,1); t_R]; 
%   
%   %Untill reproduction
%   [t_0R, E_R0R] = ode45(@dget_LR, t, LR_0, options, [], E_Hp, TC12, f_Mann, r_B, L_i, k_J, E_m, v *s_M, kap); % find state variables at start experiment
%   E_R0R(1,:) = []; E_R0R(end,:) = [];
%   %during brooding
%   t = [t_R; tt(tt(:,1) >= t_R,1)];
%   LR_0 = E_R0R(end,:);
%   LR_0(2) = 0;
%   [t_Ri, E_RRi] = ode45(@dget_LR, t, LR_0, options, [], E_Hp, TC18, f_Mann, r_B, L_i, k_J, E_m, v *s_M, kap); % find state variables at start experiment
%   E_RRi(1,:) = [];
%   LR = [E_R0R; E_RRi];
%   EtWd_Mann_18 = LR(:,1).^3 * (1 + f_Mann * ome) * d_V + LR(:,2)/ mu_E * w_E; % g, dry weight
%    %21°C
%   [tau_s, tau_j, tau_p, tau_b, l_s, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_ts(pars_ts, f_Mann);
%  L_i = L_m * l_i; s_M = l_j/ l_s; r_B = rho_B * k_M;    
%    options = odeset('AbsTol',1e-8, 'RelTol',1e-8);   
%   LR_0 = [L_p 0];
%   t_R = spawn.tWd_Mann_21+t_Mann;   % d, day at spawning
%   tt = tWd_Mann_21 + t_Mann;        % d, since puberty
%   t = [0; tt(tt(:,1) < t_R,1); t_R]; 
%   %Untill reproduction
%   [t_0R, E_R0R] = ode45(@dget_LR, t, LR_0, options, [], E_Hp, TC12, f_Mann, r_B, L_i, k_J, E_m, v *s_M, kap); % find state variables at start experiment
%   E_R0R(1,:) = []; E_R0R(end,:) = [];
%   %during brooding
%   t = [t_R; tt(tt(:,1) >= t_R,1)];
%   LR_0 = E_R0R(end,:);
%   LR_0(2) = 0;
%   [t_Ri, E_RRi] = ode45(@dget_LR, t, LR_0, options, [], E_Hp, TC21, f_Mann, r_B, L_i, k_J, E_m, v *s_M, kap); % find state variables at start experiment
%   E_RRi(1,:) = [];
%   LR = [E_R0R; E_RRi];
%   EtWd_Mann_21 = LR(:,1).^3 * (1 + f_Mann * ome) * d_V + LR(:,2)/ mu_E * w_E; % g, dry weight

  % Wils1987      
    L_0Wils = (0.01/ d_V / (1 + f_Wils * ome) ).^(1/3);
  tT = temp.tWd_Wils;
  [tau_s, tau_j, tau_p, tau_b, l_s, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_ts(pars_ts, f_Wils);
  L_i = L_m * l_i; s_M = l_j/ l_s; r_B = rho_B * k_M;    
  LHR_0 = [L_0Wils; E_H0; 0];  % state variables at 0
  options = odeset('Events', @puberty, 'AbsTol',1e-8, 'RelTol',1e-8);    
  [t, LHR] = ode45(@dget_LHR, tWd_Wils(:,1), LHR_0, options, E_Hp, tT, T_ref, pars_T, f_Wils, r_B, L_i, k_J, E_m, v * s_M, kap); % find state variables at start experiment
  L = LHR(:,1); E_R = LHR(:,3);
%   del_M = find_del_M(L, del_c, del_Mu, del_Mb);
%   EtL_Wils = L./ del_M;
  EtWd_Wils = L.^3 .* (1 + f_Wils * ome) * d_V + E_R * w_E/ mu_E;   % g, tissue dry weight
        
  
  % Robe1991
  %Robe1
  L_0Robe = (0.08966/ d_V / (1 + f_Robe0 * ome) ).^(1/3);
  tT = temp.tWd_Robe1;
  [tau_s, tau_j, tau_p, tau_b, l_s, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_ts(pars_ts, f_Robe1);
  L_i = L_m * l_i; s_M = l_j/ l_s; r_B = rho_B * k_M;    
  LHR_0 = [L_0Robe; E_Hp; ERrobe];  % state variables at 0
  options = odeset('AbsTol',1e-7, 'RelTol',1e-7);
  [t, LHR] = ode45(@dget_LHR, tWd_Robe1(:,1), LHR_0, options, E_Hp, tT, T_ref, pars_T, f_Robe1, r_B, L_i, k_J, E_m, v * s_M, kap); % find state variables at start experiment
  L = LHR(:,1); E_R = LHR(:,3);
  EtWd_Robe1 = L.^3 .* (1 + f_Robe1 * ome) * d_V + E_R * w_E/ mu_E;   % g, tissue dry weight
  
  %Robe2
  tT = temp.tWd_Robe2;
  [tau_s, tau_j, tau_p, tau_b, l_s, l_j, l_p, l_b, l_i, rho_j, rho_B] = get_ts(pars_ts, f_Robe2);
  L_i = L_m * l_i; s_M = l_j/ l_s; r_B = rho_B * k_M;    
  LHR_0 = [L_0Robe; E_Hp; ERrobe];  % state variables at 0
  options = odeset('AbsTol',1e-7, 'RelTol',1e-7);
  [t, LHR] = ode45(@dget_LHR, tWd_Robe2(:,1), LHR_0, options, E_Hp, tT, T_ref, pars_T, f_Robe2, r_B, L_i, k_J, E_m, v * s_M, kap); % find state variables at start experiment
  L = LHR(:,1); E_R = LHR(:,3);
  EtWd_Robe2 = L.^3 .* (1 + f_Robe2 * ome) * d_V + E_R * w_E/ mu_E;   % g, tissue dry weight
 
% OTHER DATA
% L-N data Walne
    L = (LN(:,1)/ d_V / (1 + f_LN * ome) ).^(1/3);
    ELN = 365 .*   TC_Walne.* reprod_rate_s(L, f_LN, pars_R);

% Wd-JO data
  p_ref     = p_Am * L_m^2;                                   % J/d,   max assimilation power at max size and T_ref
  pars_p    = [kap; kap_R; g; k_J; k_M; L_T; v; U_Hb; U_Hs; U_Hj; U_Hp];
  O2M       = (- n_M \ n_O)';                             % -,     matrix that converts organic to mineral fluxes. O2M is prepared for post-multiplication (eq. 4.35)
  %
  L         = (1/ d_V / (1 + f_Hutc * ome) ).^(1/3);  % cm,    structural length
  pACSJGRD  = p_ref * scaled_power_s(L, f_Hutc, pars_p, l_b, l_s, l_j, l_p);
  pADG      = pACSJGRD(:,[1 7 5]);  %pADG(:,1) = 0;       %        starved animals (for 3 days), no assimilation
  J_M       = pADG * eta_O' * O2M;                        % mol/d  mineral fluxes
  EJO       = - TC_Hutc * J_M(:,3) * 32e3 / 24;           % mg02/h,  oxygen consumption rate
  EJO       = EJO/13.3 * 10 *10;                              % uL / h, oxygen consumption 


 % pack to output
%Larval growth
  prdData.tL_Robe_15 = EtL_Robe_15;
  prdData.tL_Robe_20 = EtL_Robe_20;
  prdData.tL_Robe_25 = EtL_Robe_25;
  prdData.tL_Robe_30 = EtL_Robe_30;
  
%post-larval growth
  prdData.tL_Laba    = EtL_Laba;
  prdData.tWd_Laba   = EtWd_Laba;

%adult growth
  prdData.tL_Saur = EtL_Saur;
   prdData.tWd_Mann_12 = EtWd_Mann_12;
%   prdData.tWd_Mann_15 = EtWd_Mann_15;
%   prdData.tWd_Mann_18 = EtWd_Mann_18;
%  prdData.tWd_Mann_21 = EtWd_Mann_21;
  
  prdData.Wd_clutch12 = EWd_clutch12;
  %prdData.tL_Wils  = EtL_Wils;
  prdData.tWd_Wils = EtWd_Wils;
  prdData.tWd_Robe1 = EtWd_Robe1;
  prdData.tWd_Robe2 = EtWd_Robe2; 
  
%Other data
   prdData.LN = ELN;
  prdData.WdJO = EJO;

end

function dE_R = dget_E_R(t, E_R, f, L_0, L_i, r_B, v, E_m, kap, p_J)
  L = L_i - (L_i - L_0) * exp(- r_B * t); % struc length at t
  r = 3 * r_B * (L_i/ L - 1); % 1/d, spec growth rate
  p_C = f * E_m * L^3 * (v/ L - r); % J/d, reserve mobilisation rate
  dE_R = (1-kap) * p_C - p_J; % allocation to reprod buffer
end

% event find L
function [value,isterminal,direction] = find_L(t, LR, L_stop, varargin)
  value = LR(1) - L_stop;  % trigger 
  isterminal = 0;  % proceed at event
  direction  = [];         % get all the zeros
end

function dLR = dget_LR(t, LR, L_stop, E_Hp, TC, f, r_B, L_i, k_J, E_m, v, kap)
   L = LR(1);  E_R = LR(2); % unpack variables
 
   rT_B = r_B * TC; kT_J = k_J * TC; vT = TC * v;
   rT = 3 * rT_B * (L_i/ L - 1); % 1/d, spec growth rate
   pT_C = f * E_m * L^3 * (vT/ L - rT); % J/d, reserve mobilisation rate
   dL = L * rT/ 3; % cm/d, change in structural length L
   dE_R = (1 - kap) * pT_C - E_Hp * kT_J; % J/d, change in reproduction buffer
    
  dLR = [dL; dE_R]; % pack output
end

function dLHR = dget_LHR(t, LHR, E_Hp, tT, T_ref, pars_T, f, r_B, L_i, k_J, E_m, v, kap)
   L = LHR(1);  E_H = LHR(2); E_R = LHR(3); % unpack variables
 
   T = spline1(t,tT);
   TC = tempcorr(T, T_ref, pars_T); 

   rT_B = r_B * TC; kT_J = k_J * TC; vT = TC * v;
   rT = 3 * rT_B * (L_i/ L - 1); % 1/d, spec growth rate
   pT_C = f * E_m * L^3 * (vT/ L - rT); % J/d, reserve mobilisation rate
   dL = L * rT/ 3; % cm/d, change in structural length L
   
  if E_H < E_Hp
      dE_H = (1 - kap) * pT_C - E_H * kT_J;  % J/d, change in maturity
      dE_R = 0;
  else
      dE_H = 0;
      dE_R = (1 - kap) * pT_C - E_Hp * kT_J; % J/d, change in reproduction buffer
  end
    
  dLHR = [dL; dE_H ;dE_R]; % pack output
end

% event puberty
function [value,isterminal,direction] = puberty(t, LHR, E_Hp, varargin)
  value = LHR(2) - E_Hp;  % trigger 
  isterminal = 0;  % proceed at event
  direction  = [];         % get all the zeros
end


% full DEB function
function deLRH = dget_eLRH(t, eLRHS, f, Tcorr, v, L_s, L_j, L_T, L_m, g, kap, E_m, k_J, E_Hp)
    % created at 2020/02/19 by Brecht Stechele
      
% find initial conditions
  ee = eLRHS(1); L = eLRHS(2); ER = eLRHS(3); EH = eLRHS(4);
  
% acceleration factor
  s_M = max(1,min(L/L_s,L_j/L_s));  
% Calculate parameters to be used in differential equations  
  vT = v * Tcorr * s_M; 
  kTJ = k_J * Tcorr ; 
  L_m_j = L_m * s_M;
  
% Calculate specific growth rate and mobilisation
  r = vT * (ee/L - (1 + L_T/ L)/ L_m_j)/ (ee + g); % 1/d, spec growth rate
  p_C = (vT/ L - r) * ee * E_m * L^3;              % J/d, mobilisation

% Differential equations  DEB
  de = (f - ee) .* vT/ L; % 1/d, change in scaled reserve density e
  dL = L .* r/ 3; % cm/d, change in structural length L
  if EH < E_Hp
      dEH = (1 - kap) * p_C - EH * kTJ; % J/d, change in maturity buffer 
      dER = 0;
  else
      dEH = 0;
      dER = (1 - kap) * p_C - E_Hp * kTJ; % J/d, change in reproduction buffer
  end  
  
  %Pack to output
  deLRH = [de; dL; dER; dEH]; 
end

function [ del_M] = find_del_M(L, del_c, del_Mu, del_Mb)
del_M = del_Mb * exp(-del_c * L) + del_Mu;
end

