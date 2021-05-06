function [data, auxData, metaData, txtData, weights] = mydata_Ostrea_edulis 

%% set metaData
metaData.phylum     = 'Mollusca'; 
metaData.class      = 'Bivalvia'; 
metaData.order      = 'Ostreoida'; 
metaData.family     = 'Ostreidae';
metaData.species    = 'Ostrea_edulis'; 
metaData.species_en = 'European flat oyster'; 
metaData.ecoCode.climate = {'MC'};
metaData.ecoCode.ecozone = {'MANW','MAm', 'MAn'};
metaData.ecoCode.habitat = {'jiMcb', 'jiMr', 'jiMi'};
metaData.ecoCode.embryo  = {'Mp'};
metaData.ecoCode.migrate = {};
metaData.ecoCode.food    = {'biPp'};
metaData.ecoCode.gender  = {'Hsb'};
metaData.ecoCode.reprod  = {'O'};
metaData.T_typical  = C2K(20); % K. body temp
metaData.data_0     = {'Lw0'; 'Lb'; 'am'; 'Ls'; 'Lm'; 'Wwp'};
metaData.data_1     = {'t-L', 't-Wd', 'L-N', 'Wd-JO', 'Wd-CR'}; 

metaData.COMPLETE = 3.5; % using criteria of LikaKear2011

metaData.author   = {'Bram Elbers'};                             
metaData.date_subm = [2013 08 20];                           
metaData.email    = {'Bram.Elbers@gmail.com'};                 
metaData.address  = {'VU University amsterdam. 1081 HV. The Netherlands'};     

metaData.author_mod_1   = {'Bas Kooijman'};    
metaData.date_mod_1     = [2013 12 22];              
metaData.email_mod_1    = {'bas.kooijman@vu.nl'};            
metaData.address_mod_1  = {'VU University Amsterdam'};   

metaData.author_mod_2   = {'Brecht Stechele'};    
metaData.date_mod_2     = [2020 03 25];              
metaData.email_mod_2    = {'brecht.stechele@ugent.be'};            
metaData.address_mod_2  = {'Ghent University, Coupure Links 653, 9000 Gent, Belgium'};  

metaData.curator     = {'Bas Kooijman'};
metaData.email_cur   = {'bas.kooijman@vu.nl'}; 
metaData.date_acc    = [2020 03 26]; 

%% set data
% zero-variate data
data.am = 12777;    units.am = 'd';    label.am = 'life span';                           bibkey.am = 'Saur2003'; 
    temp.am = C2K(15);  units.temp.am = 'K'; label.temp.am = 'temperature';
data.ab   = 1.5;      units.ab      = 'd'; label.ab      = 'age at birth';         bibkey.ab = {'Gals1964','Carr2009'};
  temp.ab = C2K(15);  units.temp.ab = 'K'; label.temp.ab = 'temperature';
data.ar   = 10;       units.ar      = 'd'; label.ar      = 'age at release'; bibkey.aj = {'Gals1964','Carr2009'};
  temp.ar = C2K(15);  units.temp.ar = 'K'; label.temp.ar = 'temperature';
    %comment.ar = 'include brooding time of 10days';
    
data.Lb  = 0.010;      units.Lb  = 'cm';   label.Lb  = 'length at birth';               bibkey.Lb  = 'Buro1985';
data.Lr = 0.0183;       units.Lr  = 'cm';   label.Lr  = 'length at release';             bibkey.Lr  = 'Robe2017';  
data.Ls = 0.0300;       units.Ls = 'cm';    label.Ls = 'length at eyespot';              bibkey.Ls  = 'Rodr2019';   

% data.Ww0 = 1.44e-6;      units.Ww0  = 'g';   label.Ww0  = 'initial wet weight';          bibkey.Ww0 = 'Bayn2017'; 
%   comment.Ww0 = 'based on egg diameter of 0.014 cm: pi/6*0.014^3;';
data.Wdr = 1.7e-7*1.2;      units.Wdr  = 'g';   label.Wdr  = 'weight at release';            bibkey.Wdr = 'Laba1999'; 
    %comment.Wdr = 'times 1.2 to include inorganic content.'
data.Wds = 2.34e-6*1.2;      units.Wds  = 'g';   label.Wds  = 'weight at settlement';            bibkey.Wds = 'Laba1999'; 
    %comment.Wdr = 'times 1.2 to include inorganic content.'

data.Wwp = 0.3864/2;      units.Wwp = 'g';    label.Wwp = 'wet weight at puberty';         bibkey.Wwp = 'Cole';  
    %comment.Wwp = 'might invest in puberty before riping of the oocytes'
data.Wdi = 15.886;      units.Wdi = 'g';    label.Wdi = 'ultimate dry weight';           bibkey.Wdi = 'Saur2003';
data.CR_max = 12.6;     units.CR_max = 'L/h/g'; label.CR_max='maximum clearance rate';   bibkey.CR_max = 'Niel2017';   
 data.Wd_clutch12 = 0.150;   units.Wd_clutch12 = 'g'; label.Wd_clutch12 = 'weight loss during reproductin'; bibkey.Wd_clutch12 = 'Mann1975';
  %2.34
% uni-variate data 

% LARVAL GROWTH 
% tL data at 15°C
data.tL_Robe_15 = [1	2	4	6	8	10	11	12	13	15;
    0.018325	0.019178	0.020131	0.020834	0.022988	0.023991	0.024617	0.024594	0.02482	0.024898]';
units.tL_Robe_15 = {'d', 'cm'}; label.tL_Robe_15 = {'time since release', 'length'};  
temp.tL_Robe_15    = C2K(15);  units.temp.tL_Robe_15 = 'K'; label.temp.tL_Robe_15 = 'temperature';
bibkey.tL_Robe_15 = 'Robe2017';
% tL data at 20°C
data.tL_Robe_20 = [ 1	2	4	6	7	8	9	10; 
    0.018325	0.020353	0.023131	0.024809	0.026261	0.026713	0.027964	0.028316]';
units.tL_Robe_20 = {'d', 'cm'}; label.tL_Robe_20 = {'time since release', 'length'};  
temp.tL_Robe_20    = C2K(20);  units.temp.tL_Robe_20 = 'K'; label.temp.tL_Robe_20 = 'temperature';
bibkey.tL_Robe_20 = 'Robe2017';
% % tL data at 25°C
data.tL_Robe_25 = [ 1	2	4	6	7	8; 
    0.018325	0.020853	0.025081	0.029509	0.031286	0.032313]';
units.tL_Robe_25 = {'d', 'cm'}; label.tL_Robe_25 = {'time since release', 'length'};  
temp.tL_Robe_25    = C2K(25);  units.temp.tL_Robe_25 = 'K'; label.temp.tL_Robe_25 = 'temperature';
bibkey.tL_Robe_25 = 'Robe2017';
% tL data at 30°C
data.tL_Robe_30 = [ 1	2	4	6	7	8; 
    0.018325	0.021228	0.026581	0.029983	0.031684	0.032763]';
units.tL_Robe_30 = {'d', 'cm'}; label.tL_Robe_30 = {'time since release', 'length'};  
temp.tL_Robe_30    = C2K(30);  units.temp.tL_Robe_30 = 'K'; label.temp.tL_Robe_30 = 'temperature';
bibkey.tL_Robe_30 = 'Robe2017';

% LARVAE - to - POST-LARVAL GROWTH
% tL data at 16°C
data.tL_Laba = [0 7 11 17 18 20 23 27; 0.01736 0.02160 0.02490 0.0303 0.0345 0.0526 0.0727 0.0986 ]';
units.tL_Laba   = {'d', 'cm'};  label.tL_Laba = {'time since release', 'length'};  
temp.tL_Laba   = C2K(16);  units.temp.tL_Laba = 'K'; label.temp.tL_Laba = 'temperature';
bibkey.tL_Laba = 'Laba1999';

% tWd data at 16°C
data.tWd_Laba = [0 7 11 17 18 20 23 27; 1.7e-7 5.7e-7 1.15e-6 2.34e-6 2.26e-6 2.04e-6 4.65e-6 8.09e-6]'; 
units.tWd_Laba   = {'d', 'g'};  label.tWd_Laba = {'time since release', 'tissue dry weight'};  
temp.tWd_Laba   = C2K(16);  units.temp.tWd_Laba = 'K'; label.temp.tWd_Laba = 'temperature';
bibkey.tWd_Laba = 'Laba1999';


% ADULT GROWTH
% tL data
data.tL_Saur = [1	366	731	1096	1461	1826	2191	2556	2921	3286	3651	4016	4381	5111	5476	5841	6206	6571	7301	7666;
   5.77	6.65	8.02	7.84	8.14	8.77	9.14	9.38	9.46	10.29	10.01	9.63	10.47	11.09	10.81	11.26	11.49	11.5	11.6	12.19]';
units.tL_Saur = {'d', 'cm'}; label.tL_Saur = {'time since start experiment', 'length'};  
temp.tL_Saur    = C2K(12.4);  units.temp.tL_Saur = 'K';  label.temp.tL_Saur = 'temperature';
bibkey.tL_Saur = 'Saur2003';

% tWd data
Mann1975 = [...
1       95      95      95      95
21      171     140     189     211
35      201     322     321 	325
49      241 	327     419 	447
63      304     352     448     478
77      334     435 	412     456
91      347     430     442     375
105     442     521     529     591
119     404     525     650     838
133     444 	566     617     756];

% % %t- dW
data.tWd_Mann_12 = Mann1975(:,[1 2]);
data.tWd_Mann_12(:,2) = data.tWd_Mann_12(:,2)./1000;
units.tWd_Mann_12 = {'d', 'g'}; label.tWd_Mann_12 = {'time since start experiment', 'tissue dry weight'};
temp.tWd_Mann_12 = C2K(12);  units.temp.tWd_Mann_12 = 'K'; label.temp.tWd_Mann_12 = 'temperature';
bibkey.tWd_Mann_12 = 'Mann1975';
spawn.tWd_Mann_12 = 105.1; units.spawn.tWd_Mann_12 = 'd'; label.spawn.tWd_Mann_12 = 'day of spawning';

% data.tWd_Mann_15 = Mann1975(:,[1 3]);
% data.tWd_Mann_15(:,2) = data.tWd_Mann_15(:,2)./1000;
% units.tWd_Mann_15 = {'d', 'g'}; label.tWd_Mann_15 = {'time since start experiment', 'tissue dry weight'};
% temp.tWd_Mann_15 = C2K(15);  units.temp.tWd_Mann_15 = 'K'; label.temp.tWd_Mann_15 = 'temperature';
% bibkey.tWd_Mann_15 = 'Mann1975';
% spawn.tWd_Mann_15 = 77.1; units.spawn.tWd_Mann_15 = 'd'; label.spawn.tWd_Mann_15 = 'day of spawning';
% 
% data.tWd_Mann_18 = Mann1975(:,[1 4]);
% data.tWd_Mann_18(:,2) = data.tWd_Mann_18(:,2)./1000;
% units.tWd_Mann_18 = {'d', 'g'}; label.tWd_Mann_18 = {'time since start experiment', 'tissue dry weight'};
% temp.tWd_Mann_18 = C2K(18);  units.temp.tWd_Mann_18 = 'K'; label.temp.tWd_Mann_18 = 'temperature';
% bibkey.tWd_Mann_18 = 'Mann1975';
% spawn.tWd_Mann_18 = 49.1; units.spawn.tWd_Mann_18 = 'd'; label.spawn.tWd_Mann_18 = 'day of spawning';
% 
% data.tWd_Mann_21 = Mann1975(:,[1 5]);
% data.tWd_Mann_21(:,2) = data.tWd_Mann_21(:,2)./1000;
% units.tWd_Mann_21 = {'d', 'g'}; label.tWd_Mann_21 = {'time since start experiment', 'tissue dry weight'};
% temp.tWd_Mann_21 = C2K(21);  units.temp.tWd_Mann_21 = 'K'; label.temp.tWd_Mann_21 = 'temperature';
% bibkey.tWd_Mann_21 = 'Mann1975';
% spawn.tWd_Mann_21 = 65.1; units.spawn.tWd_Mann_21 = 'd'; label.spawn.tWd_Mann_21 = 'day of spawning';
% 
%tL and tWd data at varying temp and transition juv -> adult
% data.tL_Wils = [1	36	58	91	123	172	281	335	382	464	528; 
%     1.03	2.04	2.91	3.49	4	3.74	3.75	3.91	4.87	6.18	6.38]';
% units.tL_Wils = {'d', 'cm'};  label.tL_Wils = {'time since start experiment', 'length'};
% bibkey.tL_Wils = {'Wils1987'};
% temp.tL_Wils = [1 40 56	91	158	195	293	321	343	355	376	386	405	422	443	450	458	504	528; 
%         10.01 16.47 16.52 15.94 11.14 8.83 8.65 12.02	11.02	12.45	12.72	13.37	15.81	15.91	16.77	16.13	9.9	8.61	7.82]'; % Temperature data (°C)
% temp.tL_Wils(:,2) = C2K(temp.tL_Wils(:,2)); units.temp.tL_Wils = {'d','K'}; label.temp.tL_Wils = 'temperature';
% % 
data.tWd_Wils = [1	36	58	91	123	172	281	335	382	464	528; 
    0.01	0.02	0.02	0.04	0.06	0.05	0.05	0.1	0.21	0.28	0.64]';
units.tWd_Wils = {'d', 'g'};  label.tWd_Wils = {'time since start experiment', 'tissue dry weight'};
bibkey.tWd_Wils = {'Wils1987'};
temp.tWd_Wils = [1 40 56	91	158	195	293	321	343	355	376	386	405	422	443	450	458	504	528; 
        10.01 16.47 16.52 15.94 11.14 8.83 8.65 12.02	11.02	12.45	12.72	13.37	15.81	15.91	16.77	16.13	9.9	8.61	7.82]'; % Temperature data (°C)
temp.tWd_Wils(:,2) = C2K(temp.tWd_Wils(:,2)); units.temp.tWd_Wils = {'d','K'}; label.temp.tWd_Wils = 'temperature';

% t- Wd
Robert1991=[...
65      0.08966 	0.09310
106     0.18786    	0.17926
192     0.41011     0.24143
277     0.35007     0.31737
379     0.48974     0.31418
465     0.40905     0.37119
583     0.73465     0.43861
638     0.61086     0.47833
];
data.tWd_Robe1 = Robert1991(:,[1 2]);
units.tWd_Robe1 = {'d', 'g'}; label.tWd_Robe1 = {'time since start experiment', 'tissue dry weight'};
temp.tWd_Robe1 = [1	65	100	130	160	190	250	280	330	360	380	420	450	480	510	540	570	610	640	670;
13.93	16.93	17.74	18.5	15.48	13.95	11.28	12.21	12.19	18.75	18.01	19.73	20.71	19.18	19.06	14.99	13.28	13.12	6.58	7.07]'; % Temperature data (°C)
temp.tWd_Robe1(:,2) = C2K(temp.tWd_Robe1(:,2)); units.temp.tWd_Robe1 = {'d','K'}; label.temp.tWd_Robe1 = 'temperature';
bibkey.tWd_Robe1 = 'Robe1991';
%
data.tWd_Robe2 = Robert1991(:,[1 3]);
units.tWd_Robe2 = {'d', 'g'}; label.tWd_Robe2 = {'time since start experiment', 'tissue dry weight'};
temp.tWd_Robe2 = [1	65	100	130	160	190	250	280	330	360	380	420	450	480	510	540	570	610	640	670;
15.93	19.11	20.45	19.98	15.5	13.04	10.05	10.05	11.2	16.84	19.79	20.55	23.11	20.97	18.29	16.9	11.41	2.71	4.02	8.85]'; % Temperature data (°C)
temp.tWd_Robe2(:,2) = C2K(temp.tWd_Robe2(:,2)); units.temp.tWd_Robe2 = {'d','K'}; label.temp.tWd_Robe2 = 'temperature';
bibkey.tWd_Robe2 = 'Robe1991';


% % OTHER DATA
% % LN data

data.LN = [0.13	0.29	0.30	0.38	0.39	0.41	0.43	0.38	0.62	0.52	0.41	0.42	0.52	0.61	0.69	0.74	0.82	0.67	0.64	0.53	0.51	0.48	0.85	0.77	0.78	1.34	1.45	1.23	1.28	1.09	0.94	1.04	1.10	1.33	1.48	1.44	1.34	1.63	1.64	1.68	1.69	1.66	1.58	1.62	2.10	2.27	2.30	1.98	2.01	1.75	1.81	1.75	1.75	1.90	1.91	2.06	2.05	2.29	2.25	2.40	2.38	2.64	2.62	2.62	2.69	2.74	2.14	1.95	1.74	1.34	1.62	1.49	1.30	1.40	2.24	2.11	1.82	1.66	1.47	1.28	1.14	1.01	1.06	1.00	2.26	2.35	2.80;
0.12	0.22	0.24	0.22	0.25	0.27	0.29	0.29	0.30	0.32	0.33	0.35	0.35	0.38	0.42	0.44	0.45	0.47	0.51	0.50	0.45	0.41	0.55	0.59	0.62	0.74	0.73	0.76	0.81	0.79	0.86	0.98	1.02	0.90	0.93	0.96	1.02	1.03	1.06	0.99	0.96	0.77	0.84	0.69	0.96	0.99	1.04	1.06	1.10	1.07	1.24	1.26	1.30	1.32	1.36	1.36	1.26	1.29	1.41	1.27	1.30	1.25	1.50	1.59	1.61	1.54	1.62	1.68	1.41	1.44	1.25	1.26	1.23	1.14	0.82	0.78	0.73	0.55	0.36	0.48	0.47	0.54	0.63	0.73	2.07	2.10	2.26]';
data.LN(:,2)=data.LN(:,2)*1000000;
units.LN = {'g', '#'};  label.LN = {'dry weight tissue', 'fecundity'};
temp.LN = C2K(20);  units.temp.LN = 'K';  label.temp.LN = 'temperature';
bibkey.LN = 'Walne1964';
  
% Wd-JO and Wd-clearance rate data
data.WdJO = [...
5,  0.367
10, 0.483
15, 0.974
20, 1.422
25, 2.056];
units.WdJO = {'K', 'uL/h.g'};  label.WdJO = {'temperature', 'O_2 consumption/gDW'};
temp.WdJO = C2K(10);  units.temp.WdJO = 'K'; label.temp.WdJO = 'temperature';
bibkey.WdJO = 'Hutc1992';

%% set weights for all real data
 weights = setweights(data, []);  
 weights.Ls = weights.Ls*2;
 weights.ab = weights.ab;
 weights.ar = weights.ar*2;
 weights.Wd_clutch12 = weights.Wd_clutch12*2;
 weights.tWd_Mann_12 = 0.3 * weights.tWd_Mann_12;

 
%% set pseudodata and respective weights
 [data, units, label, weights] = addpseudodata(data, units, label, weights);
 weights.psd.k_J = 0; weights.psd.k = 0.1;
 data.psd.k = 0.6; units.psd.k  = '-'; label.psd.k  = 'maintenance ratio';  

%% pack auxData and txtData for output
 auxData.temp = temp;
 auxData.spawn = spawn;
 txtData.units = units;
 txtData.label = label;
 txtData.bibkey = bibkey;
% txtData.comment = comment;

%% Group plots
 set1 = {'tL_Robe_30','tL_Robe_25','tL_Robe_20','tL_Robe_15'}; comment1 = {'pre-settlement growth (cm) at 30, 25, 20, 15 C'};
 set2 = {'tWd_Robe1', 'tWd_Robe1'}; comment2 = {'Tissue dry weight at two locations in Arcachon, France'};
 %set3 = {'tWd_Mann_12', 'tWd_Mann_15', 'tWd_Mann_18', 'tWd_Mann_21'}; comment3 = {'Tissue dry weight at 21, 18, 15, 12 C'};
 metaData.grp.sets = { set1, set2}; 
 metaData.grp.comment = {comment1, comment2}; 

%% Discussion points
 D1 = 'mod_2: larval release 6 to 8 d after hatch and feed within mothers  mantle cavity: https://doi.org/10.2307/1542477, https://doi.org/10.2307/1542477';
 D2 = 'mod_2: settlement is 7 d after release: https://doi.org/10.1007/BF02114680';     
 D3 = 'mod_2: model asj is used, instead of  abj'; 
 D4 = 'mod_2: max length and max age are much longer: https://doi.org/10.1007/BF02114680';
 D5 = 'mod_2: Wwi = 4.98 g is referred to https://doi.org/10.1006/jmsc.1993.1052, but this data is not mentioned in the reference';
 D6 = 'mod_2: Ri is at 9 cm, not at Li';
 D7 = 'mod_2: culture conditions of tL data are not given';
 metaData.discussion = struct('D1', D1, 'D2', D2, 'D3', D3, 'D4', D4, 'D5', D5, 'D6', D6, 'D7', D7);

%% Links
metaData.links.id_CoL = '3063d7e4904e854f23e2d5ac9861a140'; % Cat of Life
metaData.links.id_EoL = '449502'; % Ency of Life
metaData.links.id_Wiki = 'Ostrea_edulis'; % Wikipedia
metaData.links.id_Taxo = '39291'; % Taxonomicon
metaData.links.id_WoRMS = '140658'; % WoRMS
metaData.links.id_molluscabase = '140658'; % MolluscaBase

%% References
bibkey = 'Wiki'; type = 'Misc'; bib = ...
'howpublished = {\url{http://en.wikipedia.org/wiki/Ostrea_edulis}}';
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Kooy2010'; type = 'Book'; bib = [ ...  % used in setting of chemical parameters and pseudodata
'author = {Kooijman. S. A. L. M.}, ' ...
'year = {2010}, ' ...
'title  = {Dynamic Energy Budget theory for metabolic organisation}, ' ...
'publisher = {Cambridge Univ. Press. Cambridge}, ' ...
'pages = {Table 4.2 (page 150). 8.1 (page 300)}, ' ...
'howpublished = {\url{http://www.bio.vu.nl/thb/research/bib/Kooy2010.html}}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Saur2003'; type = 'techreport'; bib = [ ... 
'author = {Saurel, C. and C. A. Richardson}, ' ...
'year = {2003}, ' ...
'title  = {Age and growth analysis of native oyster bed \emph{Ostrea edulis} in {W}ales}, ' ...
'institution = {CCW Contract Science Report}, ' ...
'volume = {549}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Cole1941'; type = 'article'; bib = [ ... 
'author = {Cole, H. A.}, ' ...
'year = {1941}, ' ...
'title  = {The fecundity of \emph{Ostrea edulis}}, ' ...
'journal = {Journal of the Marine Biological Association of the United Kingdom}, ' ...o
'volume = {25}, ' ...
'doi = {10.1017/S0025315400054710}, '...
'pages = {243--260}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Robe2017'; type = 'article'; bib = [ ... 
'author = {Robert, R. and J. Vignier and B. Petton}, ' ...
'year = {2017}, ' ...
'title  = {Influence of feeding regime and temperature on development and settlement of oyster \emph{Ostrea edulis}  ({L}innaeus. 1758) larvae}, '...
'journal = {Aquaculture Research}, ' ...
'volume = {1}, ' ...
'dio = {10.1111/are.13297},' ...
'pages = {1-18}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Laba1999';type = 'Article'; bib = [ ... 
'author = {Labarta, U.}, ' ... 
'year = {1999}, ' ...
'title = {Energy. biochemical substrates and growth in the larval development. metamorphosis and postlarvae of \emph{Ostrea edulis}}, ' ...
'journal = {Journal of Experimental Marine Bioogy and Ecology}, ' ...
'volume = {238}, ' ...
'number = {2}, '...
'doi = {10.1016/S0022-0981(98)00171-3}, '...
'pages = {225-242}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Mann1975';type = 'Article'; bib = [ ... 
'author = {Mann, P. R.}, ' ... 
'year = {1975}, ' ...
'title = {Growth and biochemical composition in \emph{Ostrea edulis} and \emph{Crassostrea gigas}}, ' ...
'journal = {Nineth European Marine Biology Symposium; Aberden University press}, ' ...
'doi = {10.1017/S0025315400046208}, '...
'pages = {587-607}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Buro1985';type = 'Article'; bib = [ ... 
'author = {Buroker, N. E.}, ' ... 
'year = {1985}, ' ...
'title = {Evolutionary patterns in the family {O}streidae: Larviparity vs. oviparity}, ' ...
'journal = {Journal of Experimental Marine Biology and Ecology}, ' ...
'volume = {90}, '...;
'pages = {233-247}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Haur1998'; type = 'Article'; bib = [ ... 
'author = {Buxton, C.D. and R.C. Newell and J.G. Field}, ' ... 
'year = {1998}, ' ...
'title = {Response Surface analysis of the combined effects of Exposure and Acclimation Temperatures on Filtration, Oxygen consumption and Scope for Growth in the Oyster \emph{Ostrea edulis}: determination of allometric coefficients}, ' ...
'journal = {Marine Ecology}, ' ...
'volume = {6}, '...;
'pages = {73-82}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Wils1987'; type = 'Article'; bib = [ ... 
'author = {Wilson, J. H.}, ' ... 
'year = {1987}, ' ...
'title = {Environmental parameters controling growth of \emph{Ostrea edulis}  and \emph{Pecten maximums} in suspended culture}, ' ...
'journal = {Aquaculture}, ' ...
'volume = {64}, '...;
'pages = {119-131}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Bayn2017'; type = 'Book'; bib = [ ... 
'author = {Bayne, B. L.}, ' ... 
'year = {2017}, ' ...
'title = {Biology of Oysters}, ' ...
'publisher = {Academic Press}, ' ...
'pages = {860}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Hutc1992'; type = 'article'; bib = [ ... 
'author = {Hutchinson, S. and L. E. Hawkins}, ' ... 
'year = {1992}, ' ...
'title = {Quantification of the physiological response of the {E}uropean flat oyster, \emph{Ostrea edulis} to temperature and salinity}, ' ...
'journal = {Journal of Molluscan Studies}, ' ...
'pages = {215-226}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Rodr2019'; type = 'article'; bib = [ ... 
'author = {Rodriguez-Perez, A. and J. Mark and D. H. David and H. B. Theodore and L. F. Miller and W. G. Sanderson}, ' ... 
'year = {2019}, ' ...
'title = {Conservation and restoration of a keystone species: {U}nderstanding the settlement preferences of the {E}uropean oyster \emph{Ostrea edulis}}, ' ...
'journal = {Marine pollution bulletin}, ' ...
'volume = {138},' ...
'pages = {312--321}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Elbe2013'; type = 'Misc'; bib =[...
'author = {Elbers, B.}, '...
'year = {2013}, ' ...
'note = {Own data}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Robe1991';type = 'Article'; bib = [ ... 
'author = {Robert, R. and M. Borel and Y. Pichot and G. Trut}. ' ... 
'year = {1991}. ' ...
'title = {Growth and mortality of the {E}uropean oyster \emph{Ostrea edulis} in the {B}ay of {A}rcachon ({F}rance)}. ' ...
'journal = {Aquatic living resources}. ' ...
'number = {4}. '...
'doi = {10.1051/alr:1991028}. '...
'pages = {265--274}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Pogo2011';type = 'Article'; bib = [ ... 
'author = {Rogoda, B. and B.H. Buck and W. Hagen}. ' ... 
'year = {2011}. ' ...
'title = {Growth performance and condition of oysters (\emph{Crassostrea gigas} and \emph{Ostrea edulis}) farmed in an offshore environment ({N}orth {S}ea, {G}ermany}. ' ...
'journal = {Aquaculture}. ' ...
'number = {319}. '...;
'pages = {484-492}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];
%
bibkey = 'Niel2017';type = 'Article'; bib = [ ... 
'author = {Nielsen, M. and B.W. Hansen and B. Vismann}. ' ... 
'year = {2017}. ' ...
'title = {Feeding traits of the european flat oyster, \emph{Ostrea edulis} and the invasive Pacific oyster, \emph{Crassostrea gigas}}. ' ...
'journal = {Aquaculture}. ' ...
'number = {319}. '...;
'pages = {484-492}'];
metaData.biblist.(bibkey) = ['''@', type, '{', bibkey, ', ' bib, '}'';'];

