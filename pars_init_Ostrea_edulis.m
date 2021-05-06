function [par, metaPar, txtPar] = pars_init_Ostrea_edulis(metaData)

metaPar.model = 'asj'; 

%% reference parameter (not to be changed) 
par.T_ref = 293.15;   free.T_ref = 0;   units.T_ref = 'K';        label.T_ref = 'Reference temperature'; 

%% core primary parameters 
par.z = 1.692;           free.z    = 1;   units.z = '-';            label.z = 'zoom factor'; 
par.F_m = 14.13;         free.F_m   = 1;   units.F_m = 'l/d.cm^2';   label.F_m = '{F_m}, max spec searching rate'; 
par.kap_X = 0.8;        free.kap_X = 0;   units.kap_X = '-';        label.kap_X = 'digestion efficiency of food to reserve'; 
par.kap_P = 0.1;        free.kap_P = 0;   units.kap_P = '-';        label.kap_P = 'faecation efficiency of food to faeces'; 
par.v = 0.04018;          free.v     = 1;   units.v = 'cm/d';         label.v = 'energy conductance'; 
par.kap = 0.9583;          free.kap   = 1;   units.kap = '-';          label.kap = 'allocation fraction to soma'; 
par.kap_R = 0.95;       free.kap_R = 0;   units.kap_R = '-';        label.kap_R = 'reproduction efficiency'; 
par.p_M = 7.193;        free.p_M   = 1;   units.p_M = 'J/d.cm^3';   label.p_M = '[p_M], vol-spec somatic maint'; 
par.p_T = 0;            free.p_T   = 0;   units.p_T = 'J/d.cm^2';   label.p_T = '{p_T}, surf-spec somatic maint'; 
par.k_J = 0.002;        free.k_J   = 0;   units.k_J = '1/d';        label.k_J = 'maturity maint rate coefficient'; 
par.E_G = 2374;         free.E_G   = 0;   units.E_G = 'J/cm^3';     label.E_G = '[E_G], spec cost for structure'; 
par.E_Hb = 4.985e-5;   free.E_Hb  = 1;   units.E_Hb = 'J';         label.E_Hb = 'maturity at birth'; 
par.E_Hr = 0.0002493;   free.E_Hr  = 1;     units.E_Hr = 'J';        label.E_Hr = 'maturity at release'; 
par.E_Hs = 0.002436;   free.E_Hs  = 1;   units.E_Hs = 'J';         label.E_Hs = 'maturity at settlement'; 
par.E_Hj = 0.3849;      free.E_Hj  = 1;   units.E_Hj = 'J';         label.E_Hj = 'maturity at metam'; 
par.E_Hp = 18.11;       free.E_Hp  = 1;   units.E_Hp = 'J';         label.E_Hp = 'maturity at puberty'; 
par.h_a = 5.135e-10;    free.h_a   = 1;   units.h_a = '1/d^2';      label.h_a = 'Weibull aging acceleration'; 
par.s_G = 0.0001;       free.s_G   = 0;   units.s_G = '-';          label.s_G = 'Gompertz stress coefficient'; 

%% Temperature parameters
par.T_A = 5000;         free.T_A   = 0;   units.T_A = 'K';          label.T_A = 'Arrhenius temperature'; 
par.T_L = 285.8;        free.T_L   = 0;   units.T_L = 'K';          label.T_L = 'Low boundary of optimal temperature';
par.T_AL= 23070;        free.T_AL  = 0;   units.T_AL = 'K';          label.T_AL= 'Arrhenius temperautre of lowest boundary';
par.T_H = 302.9;        free.T_H   = 0;   units.T_H = 'K';          label.T_H = 'Low boundary of optimal temperature';
par.T_AH= 55610;       free.T_AH  = 0;   units.T_AH = 'K';          label.T_AH= 'Arrhenius temperautre of lowest boundary';

%% Conversion parameters
par.del_Mu = 0.4534;    free.del_Mu = 1;  units.del_Mu = '-';       label.del_Mu = 'shape coefficient at birth'; 
par.del_Mb = 0.7747;    free.del_Mb = 1;  units.del_Mb = '-';       label.del_Mb = 'ultimate shape coefficient';

%% Auxilary parameters related to datasets
par.E_H0 = 50;          free.E_H0    = 0;     units.E_H0 = 'J';        label.E_H0 = 'initial maturity in tL_Wils and tWd_Wils data'; 
par.ERrobe = 35.41;     free.ERrobe  = 0;    units.ERrobe ='J';       label.ERrobe = 'initial matturity in tWd_Robe1 and tWd_Robe2';
par.L_0Saur = 2.715;    free.L_0Saur = 1;   units.L_0Saur = 'cm';    label.L_0Saur = 'initial physical length in tL_Saur data'; 
par.t_Mann = 128.6;    free.t_Mann = 1;   units.t_Mann = 'd';    label.t_Mann = 'initial time before experiment Mann'; 

par.f = 1;              free.f       = 0;     units.f = '-';           label.f = 'scaled function respons, univariate'; 
par.f_Hutc = 0.2;    free.f_Hutc  = 1;    units.f_Hutc = '-';      label.f_Hutc = 'scaled functional response for WdJO data'; 
par.f_LN = 0.2021;      free.f_LN    = 1;     units.f_LN = '-';        label.f_LN = 'scaled functional response for LN data'; 
par.f_Mann = 1.544;       free.f_Mann  = 1;    units.f_Mann = '-';      label.f_Mann = 'scaled functional response for tWd_Mann data'; 
par.f_Robe0 = 0.5797;   free.f_Robe0 = 1;   units.f_Robe0 = '-';     label.f_Robe0 = 'scaled functional response for tL_Robe data'; 
par.f_Laba = 1;    free.f_Laba  = 1;   units.f_Laba = '-';     label.f_Laba = 'scaled functional response for tL_Robe data'; 
par.f_Robe1 = 0.4731;   free.f_Robe1 = 1;   units.f_Robe1 = '-';     label.f_Robe1 = 'scaled functional response for tWd_Robe1 data'; 
par.f_Robe2 = 0.3853;   free.f_Robe2 = 1;   units.f_Robe2 = '-';     label.f_Robe2 = 'scaled functional response for tWd_Robe2 data'; 
par.f_Saur = 0.60;    free.f_Saur  = 0;    units.f_Saur = '-';      label.f_Saur = 'scaled functional response for tL_adult1 data'; 
par.f_Wils = 0.9;    free.f_Wils  = 1;    units.f_Wils = '_';      label.f_Wils = 'scaled functional response for tL_Wils and tWd_Wils data'; 

%% set chemical parameters from Kooy2010 
[par, units, label, free] = addchem(par, units, label, free, metaData.phylum, metaData.class); 

%% Pack output: 
txtPar.units = units; txtPar.label = label; par.free = free; 
