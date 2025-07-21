%% Matlab script for the data analysis for technosolraingardens
% Last edited on 7/16/25 by DMW
% The first part of the script calculates variables, treatment means, and performs statistical analysis on data collected from runs 2 through 6 as presented in manuscript tables
% The second part of the script calculates treatment means within each run for each pollutant and creates the time series figures for each treatment for runs 1 through 6 as presented in manuscript figures 


%% Part I
% Import data files
filename = 'Infiltration_Tests_clean.xlsx';
sheetname = 'Matlab2to6';
dataInf = readtable(filename,'Sheet', sheetname);

%read in variables
treatment= table2array(dataInf(2:76,"Treatment"));
treatments = categorical(treatment); % change from cell to categorical
treatments_run = treatment(1:15,:);
treatments_6runs = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6]'; % update to runs_treatment
reps_run = [1 1 1 1 1 2 2 2 2 2 3 3 3 3 3]';
column = table2array(dataInf(2:76,"Column"));
run = table2array(dataInf(2:76,"Run"));
volume = table2array(dataInf(2:76,"VolumeApplied")); % ml
ROstart = table2array(dataInf(2:76,"RunoffStart")); % hr:min
ROend = table2array(dataInf(2:76,"RunoffEnd")); % hr:min
Drainage = table2array(dataInf(2:76,"DrainageStart")); % hr:min
Vol15min = table2array(dataInf(2:76,"Vol15min")); % ml
Vol1hr = table2array(dataInf(2:76,"Vol1hr")); % ml
Vol24hr = table2array(dataInf(2:76,"Vol24hr")); % ml
Vol48hr = table2array(dataInf(2:76,"Vol48hr")); % ml
Vol72hr = table2array(dataInf(2:76,"Vol72hr")); % ml
Vol96hr = table2array(dataInf(2:76,"Vol96hr")); % ml

filename = 'runs_1to5';
runs_onetofive_table= readtable(filename);
runs_onetofive = table2array(runs_onetofive_table);

% extract Time points (in hours) *ChatGPT suggested code
TimePT15min = repmat(0.25,75,1);
TimePT1hr = ones(75, 1); 
TimePT24hr = repmat(24, 75, 1);

% Concentration (mg/l) of pollutants in synthetic stormwater, including log
Total_NO3_N = 2.53; % mg L-1

Total_NH4_N = 4.61; % mg L-1

Total_P = 0.64; % mg L-1

Total_Cu = 0.08; % mg L-1

Total_Zn = 0.48; % mg L-1

Total_Pb = 0.02; % mg L-1

% load in pollutant results
% NO3-N
filename='Nitrate_Analysis_clean.xlsx';
sheetname = 'Matlab2to6';
dataNO3= readtable(filename,'Sheet', sheetname);
NO3absorbance = table2array(dataNO3(:,6));
NO3slope = table2array(dataNO3(:,7));
NO3intercept = table2array(dataNO3(:,8));

% NH4-N
filename='Ammonium_Analysis565_clean.xlsx';
sheetname = 'Matlab2to6';
dataNH4= readtable(filename,'Sheet', sheetname);
NH4absorbance = table2array(dataNH4(:,6));
NH4slope = table2array(dataNH4(:,7));
NH4intercept = table2array(dataNH4(:,8));

% ICP – P, Cu, Zn, Pb
filename = 'ICP_Results_clean.xlsx';
sheetname = 'Matlab2to6';
opts = detectImportOptions(filename, 'Sheet',sheetname,'TextType','string');
opts = setvartype (opts,{'P','P_SD','Pb','Pb_SD'},'string');
dataICP= readtable(filename, opts);

% extract treatments and runs from ICP data 
treatment_ICP = dataICP{2:76,"Treatment"};
treatments_unique = unique(treatment_ICP);
%treatments_ICP = categorical(treatment_ICP); % change from cell to categorical
run_ICP = table2array(dataICP(2:76,"Run"));
runs_unique = unique(run_ICP);
%runs_ICP = categorical(run_ICP);

% Particle size distribution (mm) of waste materials 
filename = 'SwM_MiscResults_clean';
sheetname = 'Matlab_PSD';
dataPSD= readtable(filename,'Sheet', sheetname);

% Tempe Cell --> -33 kPa
filename = 'SwM_MiscResults_clean';
sheetname = 'Matlab_TempCell';
data33kpa= readtable(filename,'Sheet', sheetname);

treatments_lab = data33kpa{:,1};
treatments_labcat = categorical(treatments_lab);
treatments_unique = unique(treatments_lab);
reps_lab = data33kpa{:,2};
reps_unique = unique(reps_lab);

% Pressure plate --> -1500 kPa
filename = 'SwM_MiscResults_clean';
sheetname = 'PressurePlate';
data1500kpa= readtable(filename,'Sheet', sheetname);

% Ksat 
filename = 'SwM_MiscResults_clean';
sheetname = 'Ksat';
dataKsat= readtable(filename,'Sheet', sheetname);

ksat_treatments = dataKsat{:,1};
ksat_treatmentscat = categorical(ksat_treatments);
ksat_rep = dataKsat{:,2};
ksat_run = dataKsat{:,3};

% effluent pH 
filename = 'SwM_MiscResults_clean';
sheetname = 'Matlab_pH2to6';
datapH= readtable(filename,'Sheet', sheetname);

% Import equipment weights
filename = 'average_weights_clean';
equip_weights= readtable(filename);

core_length = table2array(equip_weights(1:10,7)); % mm
length_core =  (mean(core_length))/10; % cm
core_diameter = table2array(equip_weights(1:10,6)); % mm
diameter_core = (mean(core_diameter))/10; % cm
radius_core = diameter_core/2; % cm 
area_core = pi*(radius_core^2); % cm^2 
tins = table2array(equip_weights(1:23,2));
weight_tin = mean(tins);  

% Calculate variables 
% Time for Drainage to start (minutes)
T2D = (Drainage - ROstart).*1440; % minutes
T2D_Table = table(treatment,run,T2D,'VariableNames',{'Treatment','Run','T2D'});

% Tables for T2D individual treatments
filter_idx = (treatments == '1');
T2D_clay =T2D(filter_idx);

filter_idx = (treatments == '2');
T2D_glass = T2D(filter_idx);

filter_idx = (treatments == '3');
T2D_sand = T2D(filter_idx);

filter_idx = (treatments == '4');
T2D_shale = T2D(filter_idx);

filter_idx = (treatments == '5');
T2D_shell = T2D(filter_idx);

% Mean T2D
Means_T2D_All = groupsummary(T2D_Table,'Treatment','mean','T2D');
Means_T2D_Runs = groupsummary(T2D_Table,["Treatment","Run"],"mean","T2D");

% Std dev T2D
std_T2Dclay = std(T2D_clay);
std_T2Dglass = std(T2D_glass);
std_T2Dsand = std(T2D_sand);
std_T2Dshale = std(T2D_shale);
std_T2Dshell = std(T2D_shell);

% infiltration rate (cm hr-1)
ROstart_hrmin = datetime(ROstart,'ConvertFrom','datenum','Format','HH:mm');
ROend_hrmin = datetime(ROend,'ConvertFrom','datenum','Format','HH:mm');
ROtime = ROend_hrmin - ROstart_hrmin;
ROtime_min = minutes(ROtime);
IR = volume./ROtime_min; % ml/min
IR_cmmin = IR*0.01; % cm min-1 (convert ml=1 to cm3=1; convert cm3=1 to cm=0.01)
IR_cmhr = IR_cmmin*60; % cm hr-1 (convert min=60 to hour=1)
IR_Table = table(treatment,run,IR_cmhr,'VariableNames',{'Treatment','Run','InfiltrationRate'}); % cm hr-1

% create Tables for individual treatments 
filter_idx = (treatments == '1');
IR_clay = IR_cmhr(filter_idx);

filter_idx = (treatments == '2');
IR_glass = IR_cmhr(filter_idx);

filter_idx = (treatments == '3');
IR_sand = IR_cmhr(filter_idx);

filter_idx = (treatments == '4');
IR_shale = IR_cmhr(filter_idx);

filter_idx = (treatments == '5');
IR_shell = IR_cmhr(filter_idx);

% Means infiltration rate – treatments
Means_IR_All = groupsummary(IR_Table,'Treatment','mean','InfiltrationRate');
Means_IR_Runs = groupsummary(IR_Table,["Treatment","Run"],"mean","InfiltrationRate");

% Std Dev IR - treatments
std_IRclay = std(IR_clay);
std_IRglass =std(IR_glass); % REPLACE/REMOVE NaN VALUES 
std_IRsand = std(IR_sand);
std_IRshale = std(IR_shale);
std_IRshell = std(IR_shell);

% drainage volume @ 15 min (ml)
DV15min_Table = table(treatment,TimePT15min,run,Vol15min,'VariableNames',{'Treatment','TimePt','Run','Vol15min'});
DV15min = table2array(DV15min_Table(:,4));

% Tables DV15min treatments
filter_idx = (treatments == '1');
DV15min_clay = DV15min(filter_idx);

filter_idx = (treatments == '2');
DV15min_glass = DV15min(filter_idx);

filter_idx = (treatments == '3');
DV15min_sand = DV15min(filter_idx);

filter_idx = (treatments == '4');
DV15min_shale = DV15min(filter_idx);

filter_idx = (treatments == '5');
DV15min_shell = DV15min(filter_idx);

% Means Drainage volume @ 15 minutes
Means_DV15min_All = groupsummary(DV15min_Table,'Treatment','mean','Vol15min');
Means_DV15min_All.TimePoint = repmat(0.25, height(Means_DV15min_All),1); % add timepoint column
Means_DV15min_All.Properties.VariableNames{'mean_Vol15min'} = 'Volume'; % change name of column 

% Std Dev DV15min
std_DV15minclay = std(DV15min_clay);
std_DV15minlass = std(DV15min_glass);
std_DV15minsand = std(DV15min_sand);
std_DV15minshale = std(DV15min_shale);
std_DV15minshell = std(DV15min_shell);

% Drainage volume @ 1 hr (ml)
DV1hr_Table = table(treatment,TimePT1hr,run,Vol1hr,'VariableNames',{'Treatment','TimePt','Run','Vol1hr'});
DV1hr = table2array(DV1hr_Table(:,4));

% Tables DV1hr individual treatments
filter_idx = (treatments == '1');
DV1hr_clay = DV1hr(filter_idx);

filter_idx = (treatments == '2');
DV1hr_glass = DV1hr(filter_idx);

filter_idx = (treatments == '3');
DV1hr_sand = DV1hr(filter_idx);

filter_idx = (treatments == '4');
DV1hr_shale = DV1hr(filter_idx);

filter_idx = (treatments == '5');
DV1hr_shell = DV1hr(filter_idx);

% Means Drainage volume @ 1 hour
Means_DV1hr_All = groupsummary(DV1hr_Table,'Treatment','mean','Vol1hr');
Means_DV1hr_All.TimePoint = repmat(0.25, height(Means_DV1hr_All),1); % add timepoint column
Means_DV1hr_All.Properties.VariableNames{'mean_Vol1hr'} = 'Volume'; % change name of column

% Std Dev DV1hr
std_DV1hrclay = std(DV1hr_clay);
std_DV1hrglass = std(DV1hr_glass);
std_DV1hrsand = std(DV1hr_sand);
std_DV1hrshale = std(DV1hr_shale);
std_DV1hrshell = std(DV1hr_shell);

% Drainage volume @ 24 hr (ml)
DV24hr_Table = table(treatment,TimePT24hr,run,Vol24hr,'VariableNames',{'Treatment','TimePt','Run','Vol24hr'});
DV24hr = table2array(DV24hr_Table(:,4));

% Tables DV24hr individual treatments
filter_idx = (treatments == '1');
DV24hr_clay = DV24hr(filter_idx);

filter_idx = (treatments == '2');
DV24hr_glass = DV24hr(filter_idx);

filter_idx = (treatments == '3');
DV24hr_sand = DV24hr(filter_idx);

filter_idx = (treatments == '4');
DV24hr_shale = DV24hr(filter_idx);

filter_idx = (treatments == '5');
DV24hr_shell = DV24hr(filter_idx);

% Means Drainage volume @ 24 hours
Means_DV24hr_All = groupsummary(DV24hr_Table,'Treatment','mean','Vol24hr');
Means_DV24hr_All.TimePoint = repmat(0.25, height(Means_DV24hr_All),1); % add timepoint column
Means_DV24hr_All.Properties.VariableNames{'mean_Vol24hr'} = 'Volume'; % change name of column 

% Std Dev DV24hr
std_DV24hrclay = std(DV24hr_clay);
std_DV24hrglass = std(DV24hr_glass);
std_DV24hrsand = std(DV24hr_sand);
std_DV24hrshale = std(DV24hr_shale);
std_DV24hrshell = std(DV24hr_shell);

% Total Volume Drained (ml)
TVD = Vol15min + Vol1hr + Vol24hr;
TVD_Table = table(treatment,run,TVD,'VariableNames',{'Treatment','Run','TotalVolDrained'});

% Tables total volume drained individual treatments
filter_idx = (treatments == '1');
TVD_clay = TVD(filter_idx);

filter_idx = (treatments == '2');
TVD_glass = TVD(filter_idx);

filter_idx = (treatments == '3');
TVD_sand = TVD(filter_idx);

filter_idx = (treatments == '4');
TVD_shale = TVD(filter_idx);

filter_idx = (treatments == '5');
TVD_shell = TVD(filter_idx);

% means water retained
Means_TVD_All = groupsummary(TVD_Table,'Treatment','mean','TotalVolDrained');
Means_TVD_Runs = groupsummary(TVD_Table,["Treatment","Run"],"mean","TotalVolDrained");

% water retained in column (ml)
WR = volume - Vol15min - Vol1hr - Vol24hr;
WR_Table = table(treatment,run,WR,'VariableNames',{'Treatment','Run','WR'});

% Tables water retained individual treatments
filter_idx = (treatments == '1');
WR_clay = WR(filter_idx);

filter_idx = (treatments == '2');
WR_glass = WR(filter_idx);

filter_idx = (treatments == '3');
WR_sand = WR(filter_idx);

filter_idx = (treatments == '4');
WR_shale = WR(filter_idx);

filter_idx = (treatments == '5');
WR_shell = WR(filter_idx);

% means water retained
Means_WR_All = groupsummary(WR_Table,'Treatment','mean','WR');
Means_WR_Runs = groupsummary(WR_Table,["Treatment","Run"],"mean","WR");

% Std Dev WR 
std_WRclay = std(WR_clay);
std_WRglass = std(WR_glass);
std_WRsand = std(WR_sand);
std_WRshale = std(WR_shale);
std_WRshell = std(WR_shell);

% NO3-N concentration
% extract NO3-N concentration from dataNO3 table 
NO3 = (NO3absorbance - NO3intercept)./NO3slope;
NO3_concentration = []; 

for i = 1:3:length(NO3)
    means = mean(NO3(i:i+2)); 
    NO3_concentration = [NO3_concentration, means];
end 

NO3_vertical = NO3_concentration'; % convert horizontal vector to vertical vector
NO3_Table = table(treatment,run,NO3_vertical,'VariableNames',{'Treatment','Run','NO3'});

% Tables NO3 individual treatments
filter_idx = (treatments == '1');
NO3_clay = NO3_vertical(filter_idx);

filter_idx = (treatments == '2');
NO3_glass = NO3_vertical(filter_idx);

filter_idx = (treatments == '3');
NO3_sand = NO3_vertical(filter_idx);

filter_idx = (treatments == '4');
NO3_shale = NO3_vertical(filter_idx);

filter_idx = (treatments == '5');
NO3_shell = NO3_vertical(filter_idx);

% Means NO3-N – treatments
Means_NO3_All = groupsummary(NO3_Table,'Treatment','mean','NO3');
Means_NO3_All.Properties.VariableNames{'mean_NO3'} = 'NO3';
Means_NO3_Runs = groupsummary(NO3_Table,["Treatment","Run"],"mean","NO3");
Means_NO3_Runs.Properties.VariableNames{'mean_NO3'} = 'NO3';

% Std Dev NO3 - treatments
std_NO3clay = std(NO3_clay);
std_NO3glass = std(NO3_glass);
std_NO3sand = std(NO3_sand);
std_NO3shale = std(NO3_shale);
std_NO3shell = std(NO3_shell);

% NH4-N concentration
NH4 = (NH4absorbance - NH4intercept)./NH4slope;

NH4_concentration = []; 

for i = 1:2:length(NH4)
    means = mean(NH4(i:i+1)); 
    NH4_concentration = [NH4_concentration, means];
end 

NH4_vertical = NH4_concentration'; % convert horizontal vector to vertical vector
NH4_Table = table(treatment,run,NH4_vertical,'VariableNames',{'Treatment','Run','NH4'});

% Tables NH4 individual treatments 
filter_idx = (treatments == '1');
NH4_clay = NH4_vertical(filter_idx);

filter_idx = (treatments == '2');
NH4_glass = NH4_vertical(filter_idx);

filter_idx = (treatments == '3');
NH4_sand = NH4_vertical(filter_idx);

filter_idx = (treatments == '4');
NH4_shale = NH4_vertical(filter_idx);

filter_idx = (treatments == '5');
NH4_shell = NH4_vertical(filter_idx);

% Means NH4-N – treatments
Means_NH4_All = groupsummary(NH4_Table,'Treatment','mean','NH4');
Means_NH4_All.Properties.VariableNames{'mean_NH4'} = 'NH4';
Means_NH4_Runs = groupsummary(NH4_Table,["Treatment","Run"],"mean","NH4");
Means_NH4_Runs.Properties.VariableNames{'mean_NH4'} = 'NH4';

% Std Dev NH4 - treatments
std_NH4clay = std(NH4_clay);
std_NH4glass = std(NH4_glass);
std_NH4sand = std(NH4_sand);
std_NH4shale = std(NH4_shale);
std_NH4shell = std(NH4_shell);

% detect outlier in glass treatment run 3
NH4_glassR3 = NH4_glass(4:6);
glassoutlierRun3 = isoutlier(NH4_glassR3);


% P concentration
% Extract P data from dataICP
P = table2array(dataICP(2:76,"P")); % units: ng/ml

% recalculate all values under detection limit (i.e. ">") to 65% DL
for i = 1:length(P)
    if contains(P(i),"<")
        val = str2double(erase(P(i),"<"));
        P(i) = string(val*0.65);
    end    
end

%convert ng ml-1 to mg l-1
P_ng = str2double(P);
P_mg = (P_ng_double/1000000)*1000; % units: mg/l
P_Table = table(treatment_ICP,run_ICP,P_mg,'VariableNames',{'Treatment','Run','P'});

% Tables P individual treatments
filter_idx = (treatments == '1');
P_clay = P_mg(filter_idx);

filter_idx = (treatments == '2');
P_glass = P_mg(filter_idx);

filter_idx = (treatments == '3');
P_sand = P_mg(filter_idx);

filter_idx = (treatments == '4');
P_shale = P_mg(filter_idx);

filter_idx = (treatments == '5');
P_shell = P_mg(filter_idx);

% Means P – treatments 
Means_P_All = groupsummary(P_Table,'Treatment','mean','P');
Means_P_All.Properties.VariableNames{'mean_P'} = 'P';
Means_P_Runs = groupsummary(P_Table,["Treatment","Run"],"mean","P");
Means_P_Runs.Properties.VariableNames{'mean_P'} = 'P';

% Std Dev P - treatments
std_Pclay = std(P_clay);
std_Pglass = std(P_glass);
std_Psand = std(P_sand);
std_Pshale = std(P_shale);
std_Pshell = std(P_shell);

% Cu concentration
% Extract Cu data from dataICP
Cu_ng = table2array(dataICP(2:76,"Cu"));
%convert ng ml-1 to mg l-1
Cu_mg = (Cu_ng/1000000)*1000; % units: mg/l
Cu_Table = table(treatment_ICP,run_ICP,Cu_mg,'VariableNames',{'Treatment','Run','Cu'});

% Tables Cu individual treatments
filter_idx = (treatments == '1');
Cu_clay = Cu_mg(filter_idx);

filter_idx = (treatments == '2');
Cu_glass = Cu_mg(filter_idx);

filter_idx = (treatments == '3');
Cu_sand = Cu_mg(filter_idx);

filter_idx = (treatments == '4');
Cu_shale = Cu_mg(filter_idx);

filter_idx = (treatments == '5');
Cu_shell = Cu_mg(filter_idx);

% Means Cu – treatments 
Means_Cu_All = groupsummary(Cu_Table,'Treatment','mean','Cu');
Means_Cu_All.Properties.VariableNames{'mean_Cu'} = 'Cu';
Means_Cu_Runs = groupsummary(Cu_Table,["Treatment","Run"],"mean","Cu");
Means_Cu_Runs.Properties.VariableNames{'mean_Cu'} = 'Cu';

% Std Dev Cu – treatments 
std_Cuclay = std(Cu_clay);
std_Cuglass = std(Cu_glass);
std_Cusand = std(Cu_sand);
std_Cushale = std(Cu_shale);
std_Cushell = std(Cu_shell);

% Zn concentration
% Extract Zn data from dataICP
Zn = table2array(dataICP(2:76,"Zn"));

%convert ng ml-1 to mg l-1
Zn = str2double(Zn_ng);
%convert ng ml-1 to mg l-1
Zn_mg = (Zn_ng/1000000)*1000; % units: mg/l
Zn_Table = table(treatment_ICP,run_ICP,Zn_mg,'VariableNames',{'Treatment','Run','Zn'});

% Tables Zn individual treatments
filter_idx = (treatments == '1');
Zn_clay = Zn_mg(filter_idx);

filter_idx = (treatments == '2');
Zn_glass = Zn_mg(filter_idx);

filter_idx = (treatments == '3');
Zn_sand = Zn_mg(filter_idx);

filter_idx = (treatments == '4');
Zn_shale = Zn_mg(filter_idx);

filter_idx = (treatments == '5');
Zn_shell = Zn_mg(filter_idx);

% Means Zn – treatments 
Means_Zn_All = groupsummary(Zn_Table,'Treatment','mean','Zn');
Means_Zn_All.Properties.VariableNames{'mean_Zn'} = 'Zn';
Means_Zn_Runs = groupsummary(Zn_Table,["Treatment","Run"],"mean","Zn");
Means_Zn_Runs.Properties.VariableNames{'mean_Zn'} = 'Zn';

% Std Dev Zn – treatments 
std_Znclay = std(Zn_clay);
std_Znglass = std(Zn_glass);
std_Znsand = std(Zn_sand);
std_Znshale = std(Zn_shale);
std_Znshell = std(Zn_shell);

% Pb concentration
% extract Pb data from dataICP
Pb = table2array(dataICP(2:76,"Pb"));

% recalculate all values under detection limit (i.e. ">") to 65% DL
for i = 1:length(Pb)
    if contains(Pb(i),"<")
        val = str2double(erase(Pb(i),"<"));
        Pb(i) = string(val*0.65);
    end    
end


%convert ng ml-1 to mg l-1
Pb_ng = str2double(Pb);
Pb_mg = (Pb_ng/1000000)*1000; % units: mg/l
Pb_Table = table(treatment_ICP,run_ICP,Pb_mg,'VariableNames',{'Treatment','Run','Pb'});

% Tables Pb individual treatments
filter_idx = (treatments == '1');
Pb_clay = Pb_mg(filter_idx);

filter_idx = (treatments == '2');
Pb_glass = Pb_mg(filter_idx);

filter_idx = (treatments == '3');
Pb_sand = Pb_mg(filter_idx);

filter_idx = (treatments == '4');
Pb_shale = Pb_mg(filter_idx);

filter_idx = (treatments == '5');
Pb_shell = Pb_mg(filter_idx);

 % Means Pb – treatments 
Means_Pb_All = groupsummary(Pb_Table,'Treatment','mean','Pb');
Means_Pb_All.Properties.VariableNames{'mean_Pb'} = 'Pb';
Means_Pb_Runs = groupsummary(Pb_Table,["Treatment","Run"],"mean","Pb");
Means_Pb_Runs.Properties.VariableNames{'mean_Pb'} = 'Pb';

% Std Dev Pb – treatments 
std_Pbclay = std(Pb_clay);
std_Pbglass = std(Pb_glass);
std_Pbsand = std(Pb_sand);
std_Pbshale = std(Pb_shale);
std_Pbshell = std(Pb_shell);

% pollutant removal rates (%)
% (Total_X - Measured_X)/Total_X * 100 = percent change 
removal_NO3 = ((Total_NO3_N - NO3_Table.NO3)/Total_NO3_N)*100;

removal_NH4 = ((Total_NH4_N - NH4_Table.NH4)/Total_NH4_N)*100;

removal_P = ((Total_P - P_Table.P)/Total_P)*100;

removal_Cu = ((Total_Cu - Cu_Table.Cu)/Total_Cu)*100;

removal_Zn = ((Total_Zn - Zn_Table.Zn)/Total_Zn)*100;

removal_Pb = ((Total_Pb - Pb_Table.Pb)/Total_Pb)*100;

% Particle Size Distribution of Waste Materials 
% extract data by waste material
PSD_Shale = dataPSD(1:2,:);
PSD_Glass = dataPSD(3:4,:);
PSD_Shell = dataPSD(5:6,:);

%average reps 1 and 2
PSD_ShaleAve = mean(PSD_Shale(:,3:8));
PSD_GlassAve = mean(PSD_Glass(:,3:8));
PSD_ShellAve = mean(PSD_Shell(:,3:8));

% extract total sieve weight and change table to array
PSD_ShaleTotal = table2array(PSD_ShaleAve(:,1));
PSD_GlassTotal = table2array(PSD_GlassAve(:,1));
PSD_ShellTotal = table2array(PSD_ShellAve(:,1));

% extract data for each size fraction and change table to array
PSD_Shale1 = table2array(PSD_ShaleAve(:,2));
PSD_Shale1to2 = table2array(PSD_ShaleAve(:,3));
PSD_Shale2to4 = table2array(PSD_ShaleAve(:,4));
PSD_Shale4to8 = table2array(PSD_ShaleAve(:,5));
PSD_Shale8 = table2array(PSD_ShaleAve(:,6));

PSD_Glass1 = table2array(PSD_GlassAve(:,2));
PSD_Glass1to2 = table2array(PSD_GlassAve(:,3));
PSD_Glass2to4 = table2array(PSD_GlassAve(:,4));
PSD_Glass4to8 = table2array(PSD_GlassAve(:,5));
PSD_Glass8 = table2array(PSD_GlassAve(:,6));

PSD_Shell1 = table2array(PSD_ShellAve(:,2));
PSD_Shell1to2 = table2array(PSD_ShellAve(:,3));
PSD_Shell2to4 = table2array(PSD_ShellAve(:,4));
PSD_Shell4to8 = table2array(PSD_ShellAve(:,5));
PSD_Shell8 = table2array(PSD_ShellAve(:,6));

% Calculate proportion of each size fraction based on total sieve weight
PSD_ShaleFrac1 = PSD_Shale1/PSD_ShaleTotal;
PSD_ShaleFrac1to2 = PSD_Shale1to2/PSD_ShaleTotal;
PSD_ShaleFrac2to4 = PSD_Shale2to4/PSD_ShaleTotal;
PSD_ShaleFrac4to8 = PSD_Shale4to8/PSD_ShaleTotal;
PSD_ShaleFrac8 = PSD_Shale8/PSD_ShaleTotal;

PSD_GlassFrac1 = PSD_Glass1/PSD_GlassTotal;
PSD_GlassFrac1to2 = PSD_Glass1to2/PSD_GlassTotal;
PSD_GlassFrac2to4 = PSD_Glass2to4/PSD_GlassTotal;
PSD_GlassFrac4to8 = PSD_Glass4to8/PSD_GlassTotal;
PSD_GlassFrac8 = PSD_Glass8/PSD_GlassTotal;

PSD_ShellFrac1 = PSD_Shell1/PSD_ShellTotal;
PSD_ShellFrac1to2 = PSD_Shell1to2/PSD_ShellTotal;
PSD_ShellFrac2to4 = PSD_Shell2to4/PSD_ShellTotal;
PSD_ShellFrac4to8 = PSD_Shell4to8/PSD_ShellTotal;
PSD_ShellFrac8 = PSD_Shell8/PSD_ShellTotal;

% combine proportions into one table 
PSD_ShaleFracTotal = [PSD_ShaleFrac1,PSD_ShaleFrac1to2,PSD_ShaleFrac2to4,PSD_ShaleFrac4to8,PSD_ShaleFrac8];
PSD_GlassFracTotal = [PSD_GlassFrac1,PSD_GlassFrac1to2,PSD_GlassFrac2to4,PSD_GlassFrac4to8,PSD_GlassFrac8];
PSD_ShellFracTotal = [PSD_ShellFrac1,PSD_ShellFrac1to2,PSD_ShellFrac2to4,PSD_ShellFrac4to8,PSD_ShellFrac8];

% VWC @ -33 kPa
VWC33_wet =  table2array(data33kpa(:,11));
VWC33_wetsoil = VWC33_wet - weight_tin;
VWC33_dry = table2array(data33kpa(:,12));
VWC33_drysoil = VWC33_dry - weight_tin;
VWC33_volume = pi * (radius_core^2) * length_core; % cm^3
VWC_bd = VWC33_drysoil/VWC33_volume; % g cm-3
VWC33 = ((VWC33_wetsoil-VWC33_drysoil)/VWC33_drysoil)*VWC_bd; % decimal 

VWC33_Table = table(treatments_lab,reps_lab,VWC33,'VariableNames',{'Treatment','Rep','VWC33'});

% Tables VWC @ 33 individual treatments
VWC33_clay = table2array(VWC33_Table(1:3,3));
VWC33_glass = table2array(VWC33_Table(4:6,3));
VWC33_sand = table2array(VWC33_Table(7:9,3));
VWC33_shale = table2array(VWC33_Table(10:12,3));
VWC33_shell = table2array(VWC33_Table(13:15,3));

% Means VWC @ -33 kPa
Means_VWC33_All = groupsummary(VWC33_Table,'Treatment','mean','VWC33');
Means_VWC33_Reps = groupsummary(VWC33_Table,["Treatment","Rep"],"mean","VWC33");

% Std Dev VWC33
std_VWC33clay = std(VWC33_clay);
std_VWC33glass = std(VWC33_glass);
std_VWC33sand = std(VWC33_sand);
std_VWC33shale = std(VWC33_shale);
std_VWC33shell = std(VWC33_shell);

% VWC @ -1500 kPa
VWC1500_wet = table2array(data1500kpa(:,5));
VWC1500_dry = table2array(data1500kpa(:,6));
VWC1500_TinWt = table2array(data1500kpa(:,4));
VWC1500_wetsoil = VWC1500_wet - VWC1500_TinWt;
VWC1500_drysoil = VWC1500_dry - VWC1500_TinWt;
VWC1500 = ((VWC1500_wetsoil-VWC1500_drysoil)/VWC1500_drysoil)*VWC_bd; % decimal 

VWC1500_Table = table(treatments_lab,reps_lab,VWC1500,'VariableNames',{'Treatment','Rep','VWC1500'});

% Tables VWC @ 1500 individual treatments
filter_idx = (treatments_lab == 1);
VWC1500_clay = VWC1500(filter_idx);

filter_idx = (treatments_lab == 2);
VWC1500_glass = VWC1500(filter_idx);

filter_idx = (treatments_lab == 3);
VWC1500_sand = VWC1500(filter_idx);

filter_idx = (treatments_lab == 4);
VWC1500_shale = VWC1500(filter_idx);

filter_idx = (treatments_lab == 5);
VWC1500_shell = VWC1500(filter_idx);

% Means VWC @ -33 kPa
Means_VWC1500_All = groupsummary(VWC1500_Table,'Treatment','mean','VWC1500');
Means_VWC1500_Reps = groupsummary(VWC1500_Table,["Treatment","Rep"],"mean","VWC1500");

% Std Dev VWC33
std_VWC1500clay = std(VWC1500_clay);
std_VWC1500glass = std(VWC1500_glass);
std_VWC1500sand = std(VWC1500_sand);
std_VWC1500shale = std(VWC1500_shale);
std_VWC1500shell = std(VWC1500_shell);

% Available Water Capacity; AWC = (VWC33 - VWC1500) *LengthSoil
FillMix_length = 45; % cm 
AWC = (VWC33 - VWC1500);
AWC_length = (VWC33 - VWC1500)*FillMix_length;
AWC_Table = table(treatments_lab,reps_lab,AWC,'VariableNames',{'Treatment','Rep','AWC'});

% Tables AWC individual treatments
filter_idx = (treatments_lab == 1);
AWC_clay = AWC(filter_idx);

filter_idx = (treatments_lab == 2);
AWC_glass = AWC(filter_idx);

filter_idx = (treatments_lab == 3);
AWC_sand = AWC(filter_idx);

filter_idx = (treatments_lab == 4);
AWC_shale = AWC(filter_idx);

filter_idx = (treatments_lab == 5);
AWC_shell = AWC(filter_idx);

% Means VWC @ -33 kPa
Means_AWC_All = groupsummary(AWC_Table,'Treatment','mean','AWC');
Means_AWC_Reps = groupsummary(AWC_Table,["Treatment","Rep"],"mean","AWC");

% Std Dev VWC33
std_AWCclay = std(AWC_clay);
std_AWCglass = std(AWC_glass);
std_AWCsand = std(AWC_sand);
std_AWCshale = std(AWC_shale);
std_AWCshell = std(AWC_shell);

% Ksat; Ks = (V * L)/(A * t * h)
ksat_V = table2array(dataKsat(:,14));
ksat_t = table2array(dataKsat(:,13)); % minutes
ksat_tanki = table2array(dataKsat(:,7));
ksat_tankf = table2array(dataKsat(:,11));
ksat_sampi = table2array(dataKsat(:,8));
ksat_sampf = table2array(dataKsat(:,12));
ksat_tankh = (ksat_tanki + ksat_tankf)/2;
ksat_samph = (ksat_sampi + ksat_sampf)/2;
ksat_h = ksat_sampf - ksat_tankf;

%Ksat = (core_length .* Ksat_V)/(Ksat_t * Ksat_h .* core_area);
ksat_num = length_core .* ksat_V;
ksat_txh = (ksat_t.*ksat_h);
ksat_denom = area_core .* ksat_txh;
ksat_cmmin = ksat_num./ksat_denom; % cm min-1
ksat_cmhr = ksat_cmmin*60; % cm hr-1
ksat_inhr = ksat_cmhr/2.54; % in hr-1
ksat_cmday = ksat_cmmin*1440; % cm day-1

ksat_Table = table(ksat_treatments,ksat_rep,ksat_run,ksat_cmhr,'VariableNames',{'Treatment','Rep','Run','Ksat'});

% Tables ksat individual treatments
filter_idx = (ksat_treatmentscat == '1');
ksat_clay = ksat_cmhr(filter_idx);

filter_idx = (ksat_treatmentscat == '2');
ksat_glass = ksat_cmhr(filter_idx);

filter_idx = (ksat_treatmentscat == '3');
ksat_sand = ksat_cmhr(filter_idx);

filter_idx = (ksat_treatmentscat == '4');
ksat_shale = ksat_cmhr(filter_idx);

filter_idx = (ksat_treatmentscat == '5');
ksat_shell = ksat_cmhr(filter_idx);

% Means Ksat
Means_ksat_All = groupsummary(ksat_Table,'Treatment','mean','Ksat');
Means_ksat_Reps = groupsummary(ksat_Table,["Treatment","Rep"],"mean","Ksat");

% Std Dev Ksat
std_ksatclay = std(ksat_clay);
std_ksatglass = std(ksat_glass);
std_ksatsand = std(ksat_sand);
std_ksatshale = std(ksat_shale);
std_ksatshell = std(ksat_shell);

% Bulk density packed cores for Ksat and tempe cell
% bd = oven dry weight / total volume
% use bulk density calculated in tempe cell procedure
BD = VWC_bd;
BD_Table = table(treatments_lab,reps_lab,BD,'VariableNames',{'Treatment','Rep','BulkDensity'});

% Tables BD individual treatments
filter_idx = (treatments_labcat == '1');
BD_clay = BD(filter_idx);

filter_idx = (treatments_labcat == '2');
BD_glass = BD(filter_idx);

filter_idx = (treatments_labcat == '3');
BD_sand = BD(filter_idx);

filter_idx = (treatments_labcat == '4');
BD_shale = BD(filter_idx);

filter_idx = (treatments_labcat == '5');
BD_shell = BD(filter_idx);

% Means BD
Means_BD_All = groupsummary(BD_Table,'Treatment','mean','BulkDensity');
Means_BD_Reps = groupsummary(BD_Table,["Treatment","Rep"],"mean","BulkDensity");

% Std Dev BD
std_BDclay = std(BD_clay);
std_BDglass = std(BD_glass);
std_BDsand = std(BD_sand);
std_BDshale = std(BD_shale);
std_BDshell = std(BD_shell);

% effluent pH 
pH_dataonly = table2array(datapH(:,5));
pH = [];
for i = 1:2:length(pH_dataonly)
    means = mean(pH_dataonly(i:i+1));
    pH = [pH,means];
end

pH_column = pH';
pH_Table = table(treatments, run, pH_column,'VariableNames',{'Treatment','Run','pH'});

% Tables pH individual treatments
filter_idx = (treatments == '1');
pH_clay = pH_column(filter_idx);

filter_idx = (treatments == '2');
pH_glass = pH_column(filter_idx);

filter_idx = (treatments == '3');
pH_sand = pH_column(filter_idx);

filter_idx = (treatments == '4');
pH_shale = pH_column(filter_idx);

filter_idx = (treatments == '5');
pH_shell = pH_column(filter_idx);

% Means pH
pH_Table.Treatment = double(pH_Table.Treatment);
pH_Table.Run = double(pH_Table.Run);
Means_pH_All = groupsummary(pH_Table,'Treatment','mean','pH');
Means_pH_Runs = groupsummary(pH_Table,["Treatment","Run"],"mean","pH");

% Std Dev pH
std_pHclay = std(pH_clay);
std_pHglass = std(pH_glass);
std_pHsand = std(pH_sand);
std_pHshale = std(pH_shale);
std_pHshell = std(pH_shell);



%% Part II
% Time series Pollutant concentration for each treatment
% load in pollutant results
% NO3-N
filename='Nitrate_Analysis_clean.xlsx';
sheetname = 'Matlab1to6';
dataNO3_R1to6= readtable(filename,'Sheet', sheetname);
NO3absorbance = table2array(dataNO3_R1to6(:,6));
NO3slope = table2array(dataNO3_R1to6(:,7));
NO3intercept = table2array(dataNO3_R1to6(:,8));

% NH4-N
filename='Ammonium_Analysis565_clean.xlsx';
sheetname = 'Matlab1to6';
dataNH4_R1to6= readtable(filename,'Sheet', sheetname);
NH4absorbance = table2array(dataNH4_R1to6(:,6));
NH4slope = table2array(dataNH4_R1to6(:,7));
NH4intercept = table2array(dataNH4_R1to6(:,8));

% ICP – P, Cu, Zn, Pb
filename = 'ICP_Results_clean.xlsx';
sheetname = 'Matlab1to6';
opts = detectImportOptions(filename, 'Sheet',sheetname,'TextType','string');
opts = setvartype (opts,{'P','P_SD','Pb','Pb_SD'},'string');
dataICP_R1to6= readtable(filename,opts);

% extract treatments and runs from ICP data 
treatment_ICP = dataICP_R1to6{2:91,"Treatment"};
treatments_unique = unique(treatment_ICP);
treatments_ICP = categorical(treatment_ICP); % change from cell to categorical
run_ICP = table2array(dataICP_R1to6(2:91,"Run"));
runs_unique = unique(run_ICP);
%runs_ICP = categorical(run_ICP);
reps_run = [1 2 3 1 2 3 1 2 3 1 2 3 1 2 3 1 2 3]';
runs_treatment = [1 1 1 2 2 2 3 3 3 4 4 4 5 5 5 6 6 6]';

% Calculate data
% extract NO3-N concentration from dataNO3 table 
NO31to6 = (NO3absorbance - NO3intercept)./NO3slope;
NO3_concentration1to6 = []; 

for i = 1:3:length(NO31to6)
    means = mean(NO31to6(i:i+2)); 
    NO3_concentration1to6 = [NO3_concentration1to6, means];
end 

NO3_vertical1to6 = NO3_concentration1to6'; % convert horizontal vector to vertical vector
NO3_Table1to6 = table(treatment_ICP,run_ICP,NO3_vertical1to6,'VariableNames',{'Treatment','Run','NO3'});

% Means NO3-N – treatments
Means_NO3_All1to6 = groupsummary(NO3_Table1to6,'Treatment','mean','NO3');
Means_NO3_All1to6.Properties.VariableNames{'mean_NO3'} = 'NO3';
Means_NO3_Runs1to6 = groupsummary(NO3_Table1to6,["Treatment","Run"],"mean","NO3");
Means_NO3_Runs1to6.Properties.VariableNames{'mean_NO3'} = 'NO3';

% Tables NO3 individual treatments
filter_idx = (treatments_ICP == '1');
NO3_clay1to6 = NO3_vertical1to6(filter_idx);

filter_idx = (treatments_ICP == '2');
NO3_glass1to6 = NO3_vertical1to6(filter_idx);

filter_idx = (treatments_ICP == '3');
NO3_sand1to6 = NO3_vertical1to6(filter_idx);

filter_idx = (treatments_ICP == '4');
NO3_shale1to6 = NO3_vertical1to6(filter_idx);

filter_idx = (treatments_ICP == '5');
NO3_shell1to6 = NO3_vertical1to6(filter_idx);

% Tables NO3 individual runs within treatments
NO3_clayR1 = NO3_clay1to6(1:3);
NO3_clayR2 = mean(NO3_clay1to6(4:6));
NO3_clayR3 = mean(NO3_clay1to6(7:9));
NO3_clayR4 = mean(NO3_clay1to6(10:12));
NO3_clayR5 = mean(NO3_clay1to6(13:15));
NO3_clayR6 = mean(NO3_clay1to6(16:18));
meanNO3_clay = [0.1322,NO3_clayR2,NO3_clayR3,NO3_clayR4,NO3_clayR5,NO3_clayR6]';

NO3_glassR1 = mean(NO3_glass1to6(1:3));
NO3_glassR2 = mean(NO3_glass1to6(4:6));
NO3_glassR3 = mean(NO3_glass1to6(7:9));
NO3_glassR4 = mean(NO3_glass1to6(10:12));
NO3_glassR5 = mean(NO3_glass1to6(13:15));
NO3_glassR6 = mean(NO3_glass1to6(16:18));
meanNO3_glass = [NO3_glassR1,NO3_glassR2,NO3_glassR3,NO3_glassR4,NO3_glassR5,NO3_glassR6]';

NO3_sandR1 = mean(NO3_sand1to6(1:3));
NO3_sandR2 = mean(NO3_sand1to6(4:6)); 
NO3_sandR3 = mean(NO3_sand1to6(7:9));
NO3_sandR4 = mean(NO3_sand1to6(10:12));
NO3_sandR5 = mean(NO3_sand1to6(13:15));
NO3_sandR6 = mean(NO3_sand1to6(16:18));
meanNO3_sand = [NO3_sandR1,NO3_sandR2,NO3_sandR3,NO3_sandR4,NO3_sandR5,NO3_sandR6]';

NO3_shaleR1 = mean(NO3_shale1to6(1:3));
NO3_shaleR2 = mean(NO3_shale1to6(4:6));
NO3_shaleR3 = mean(NO3_shale1to6(7:9));
NO3_shaleR4 = mean(NO3_shale1to6(10:12));
NO3_shaleR5 = mean(NO3_shale1to6(13:15));
NO3_shaleR6 = mean(NO3_shale1to6(16:18));
meanNO3_shale = [NO3_shaleR1,NO3_shaleR2,NO3_shaleR3,NO3_shaleR4,NO3_shaleR5,NO3_shaleR6]';

NO3_shellR1 = mean(NO3_shell1to6(1:3));
NO3_shellR2 = mean(NO3_shell1to6(4:6));
NO3_shellR3 = mean(NO3_shell1to6(7:9));
NO3_shellR4 = mean(NO3_shell1to6(10:12));
NO3_shellR5 = mean(NO3_shell1to6(13:15));
NO3_shellR6 = mean(NO3_shell1to6(16:18));
meanNO3_shell = [NO3_shellR1,NO3_shellR2,NO3_shellR3,NO3_shellR4,NO3_shellR5,NO3_shellR6]';

% NH4
NH41to6 = (NH4absorbance - NH4intercept)./NH4slope;
NH41to6_concentration = []; 

for i = 1:2:length(NH41to6)
    means = mean(NH41to6(i:i+1)); 
    NH41to6_concentration = [NH41to6_concentration, means];
end 

NH41to6_vertical = NH41to6_concentration'; % convert horizontal vector to vertical vector
NH41to6_Table = table(treatment_ICP,run_ICP,NH41to6_vertical,'VariableNames',{'Treatment','Run','NH4'});

% Means NH4-N – treatments
Means_NH4_All1to6 = groupsummary(NH41to6_Table,'Treatment','mean','NH4');
Means_NH4_All1to6.Properties.VariableNames{'mean_NH4'} = 'NH4';
Means_NH4_Runs1to6 = groupsummary(NH41to6_Table,["Treatment","Run"],"mean","NH4");
Means_NH4_Runs1to6.Properties.VariableNames{'mean_NH4'} = 'NH4';

% Tables NH4 individual treatments 
filter_idx = (treatments_ICP == '1');
NH4_clay1to6 = NH41to6_vertical(filter_idx);

filter_idx = (treatments_ICP == '2');
NH4_glass1to6 = NH41to6_vertical(filter_idx);

filter_idx = (treatments_ICP == '3');
NH4_sand1to6 = NH41to6_vertical(filter_idx);

filter_idx = (treatments_ICP == '4');
NH4_shale1to6 = NH41to6_vertical(filter_idx);

filter_idx = (treatments_ICP == '5');
NH4_shell1to6 = NH41to6_vertical(filter_idx);

% Tables NH4 individual runs within treatments
NH4_clayR1 = NH4_clay1to6(1:3);
NH4_clayR2 = mean(NH4_clay1to6(4:6));
NH4_clayR3 = mean(NH4_clay1to6(7:9));
NH4_clayR4 = mean(NH4_clay1to6(10:12));
NH4_clayR5 = mean(NH4_clay1to6(13:15));
NH4_clayR6 = mean(NH4_clay1to6(16:18));
meanNH4_clay = [0.8980,NH4_clayR2,NH4_clayR3,NH4_clayR4,NH4_clayR5,NH4_clayR6]';

NH4_glassR1 = mean(NH4_glass1to6(1:3));
NH4_glassR2 = mean(NH4_glass1to6(4:6));
NH4_glassR3 = mean(NH4_glass1to6(7:9));
NH4_glassR4 = mean(NH4_glass1to6(10:12));
NH4_glassR5 = mean(NH4_glass1to6(13:15));
NH4_glassR6 = mean(NH4_glass1to6(16:18));
meanNH4_glass = [NH4_glassR1,NH4_glassR2,NH4_glassR3,NH4_glassR4,NH4_glassR5,NH4_glassR6]';

NH4_sandR1 = mean(NH4_sand1to6(1:3));
NH4_sandR2 = mean(NH4_sand1to6(4:6));
NH4_sandR3 = mean(NH4_sand1to6(7:9));
NH4_sandR4 = mean(NH4_sand1to6(10:12));
NH4_sandR5 = mean(NH4_sand1to6(13:15));
NH4_sandR6 = mean(NH4_sand1to6(16:18));
meanNH4_sand = [NH4_sandR1,NH4_sandR2,NH4_sandR3,NH4_sandR4,NH4_sandR5,NH4_sandR6]';

NH4_shaleR1 = mean(NH4_shale1to6(1:3));
NH4_shaleR2 = mean(NH4_shale1to6(4:6));
NH4_shaleR3 = mean(NH4_shale1to6(7:9));
NH4_shaleR4 = mean(NH4_shale1to6(10:12));
NH4_shaleR5 = mean(NH4_shale1to6(13:15));
NH4_shaleR6 = mean(NH4_shale1to6(16:18));
meanNH4_shale = [NH4_shaleR1,NH4_shaleR2,NH4_shaleR3,NH4_shaleR4,NH4_shaleR5,NH4_shaleR6]';

NH4_shellR1 = mean(NH4_shell1to6(1:3));
NH4_shellR2 = mean(NH4_shell1to6(4:6));
NH4_shellR3 = mean(NH4_shell1to6(7:9));
NH4_shellR4 = mean(NH4_shell1to6(10:12));
NH4_shellR5 = mean(NH4_shell1to6(13:15));
NH4_shellR6 = mean(NH4_shell1to6(16:18));
meanNH4_shell = [NH4_shellR1,NH4_shellR2,NH4_shellR3,NH4_shellR4,NH4_shellR5,NH4_shellR6]';

% P
% Extract P data from dataICP_R1to6
P_1to6 = table2array(dataICP_R1to6(2:91,"P")); % unites: ng/ml

% recalculate all values under detection limit (i.e. "<") to 65% DL
for i = 1:length(P_1to6)
    if contains(P_1to6(i),"<")
        val = str2double(erase(P_1to6(i),"<"));
        P_1to6(i) = string(val*0.65);
    end    
end

%convert ng ml-1 to mg l-1
P_ng1to6 = str2double(P_1to6); % units: ng/ml
P_mg1to6 = (P_ng1to6/1000000)*1000; % units: mg/l
P_Table1to6 = table(treatment_ICP,run_ICP,P_mg1to6,'VariableNames',{'Treatment','Run','P'});

% Tables P individual treatments
filter_idx = (treatments_ICP == '1');
P_clay1to6 = P_mg1to6(filter_idx);

filter_idx = (treatments_ICP == '2');
P_glass1to6 = P_mg1to6(filter_idx);

filter_idx = (treatments_ICP == '3');
P_sand1to6 = P_mg1to6(filter_idx);

filter_idx = (treatments_ICP == '4');
P_shale1to6 = P_mg1to6(filter_idx);

filter_idx = (treatments_ICP == '5');
P_shell1to6 = P_mg1to6(filter_idx);

% Means P – treatments 
Means_P_All1to6 = groupsummary(P_Table1to6,'Treatment','mean','P');
Means_P_All1to6.Properties.VariableNames{'mean_P'} = 'P';
Means_P_Runs1to6 = groupsummary(P_Table1to6,["Treatment","Run"],"mean","P");
Means_P_Runs1to6.Properties.VariableNames{'mean_P'} = 'P';
 
% Tables P runs within individual treatments
P_clayR1 = P_clay1to6(1:3);
P_clayR2 = mean(P_clay1to6(4:6));
P_clayR3 = mean(P_clay1to6(7:9));
P_clayR4 = mean(P_clay1to6(10:12));
P_clayR5 = mean(P_clay1to6(13:15));
P_clayR6 = mean(P_clay1to6(16:18));
meanP_clay = [0.4344,P_clayR2,P_clayR3,P_clayR4,P_clayR5,P_clayR6]';

P_glassR1 = mean(P_glass1to6(1:3));
P_glassR2 = mean(P_glass1to6(4:6));
P_glassR3 = mean(P_glass1to6(7:9));
P_glassR4 = mean(P_glass1to6(10:12));
P_glassR5 = mean(P_glass1to6(13:15));
P_glassR6 = mean(P_glass1to6(16:18));
meanP_glass = [P_glassR1,P_glassR2,P_glassR3,P_glassR4,P_glassR5,P_glassR6]';

P_sandR1 = mean(P_sand1to6(1:3));
P_sandR2 = mean(P_sand1to6(4:6));
P_sandR3 = mean(P_sand1to6(7:9));
P_sandR4 = mean(P_sand1to6(10:12));
P_sandR5 = mean(P_sand1to6(13:15));
P_sandR6 = mean(P_sand1to6(16:18));
meanP_sand = [P_sandR1,P_sandR2,P_sandR3,P_sandR4,P_sandR5,P_sandR6]';

P_shaleR1 = mean(P_shale1to6(1:3));
P_shaleR2 = mean(P_shale1to6(4:6));
P_shaleR3 = mean(P_shale1to6(7:9));
P_shaleR4 = mean(P_shale1to6(10:12));
P_shaleR5 = mean(P_shale1to6(13:15));
P_shaleR6 = mean(P_shale1to6(16:18));
meanP_shale = [P_shaleR1,P_shaleR2,P_shaleR3,P_shaleR4,P_shaleR5,P_shaleR6]';

P_shellR1 = mean(P_shell1to6(1:3));
P_shellR2 = mean(P_shell1to6(4:6));
P_shellR3 = mean(P_shell1to6(7:9));
P_shellR4 = mean(P_shell1to6(10:12));
P_shellR5 = mean(P_shell1to6(13:15));
P_shellR6 = mean(P_shell1to6(16:18));
meanP_shell = [P_shellR1,P_shellR2,P_shellR3,P_shellR4,P_shellR5,P_shellR6]';

% Cu
% Extract Cu data from dataICP_R1to6 
Cu_ng1to6 = table2array(dataICP_R1to6(2:91,"Cu"));
%convert ng ml-1 to mg l-1
Cu_mg1to6 = (Cu_ng1to6/1000000)*1000; % units: mg/l
Cu_Table1to6 = table(treatment_ICP,run_ICP,Cu_mg1to6,'VariableNames',{'Treatment','Run','Cu'});

% Tables Cu individual treatments
filter_idx = (treatments_ICP == '1');
Cu_clay1to6 = Cu_mg1to6(filter_idx);

filter_idx = (treatments_ICP == '2');
Cu_glass1to6 = Cu_mg1to6(filter_idx);

filter_idx = (treatments_ICP == '3');
Cu_sand1to6 = Cu_mg1to6(filter_idx);

filter_idx = (treatments_ICP == '4');
Cu_shale1to6 = Cu_mg1to6(filter_idx);

filter_idx = (treatments_ICP == '5');
Cu_shell1to6 = Cu_mg1to6(filter_idx);

% Means Cu – treatments 
Means_Cu_All1to6 = groupsummary(Cu_Table1to6,'Treatment','mean','Cu');
Means_Cu_All1to6.Properties.VariableNames{'mean_Cu'} = 'Cu';
Means_Cu_Runs1to6 = groupsummary(Cu_Table1to6,["Treatment","Run"],"mean","Cu");
Means_Cu_Runs1to6.Properties.VariableNames{'mean_Cu'} = 'Cu';

% Tables Cu individual runs within treatments
Cu_clayR1 = Cu_clay1to6(1:3);
Cu_clayR2 = mean(Cu_clay1to6(4:6));
Cu_clayR3 = mean(Cu_clay1to6(7:9));
Cu_clayR4 = mean(Cu_clay1to6(10:12));
Cu_clayR5 = mean(Cu_clay1to6(13:15));
Cu_clayR6 = mean(Cu_clay1to6(16:18));
meanCu_clay = [0.0372,Cu_clayR2,Cu_clayR3,Cu_clayR4,Cu_clayR5,Cu_clayR6]'; 

Cu_glassR1 = mean(Cu_glass1to6(1:3));
Cu_glassR2 = mean(Cu_glass1to6(4:6));
Cu_glassR3 = mean(Cu_glass1to6(7:9));
Cu_glassR4 = mean(Cu_glass1to6(10:12));
Cu_glassR5 = mean(Cu_glass1to6(13:15));
Cu_glassR6 = mean(Cu_glass1to6(16:18));
meanCu_glass = [Cu_glassR1,Cu_glassR2,Cu_glassR3,Cu_glassR4,Cu_glassR5,Cu_glassR6]';

Cu_sandR1 = mean(Cu_sand1to6(1:3));
Cu_sandR2 = mean(Cu_sand1to6(4:6));
Cu_sandR3 = mean(Cu_sand1to6(7:9));
Cu_sandR4 = mean(Cu_sand1to6(10:12));
Cu_sandR5 = mean(Cu_sand1to6(13:15));
Cu_sandR6 = mean(Cu_sand1to6(16:18));
meanCu_sand = [Cu_sandR1,Cu_sandR2,Cu_sandR3,Cu_sandR4,Cu_sandR5,Cu_sandR6]';

Cu_shaleR1 = mean(Cu_shale1to6(1:3));
Cu_shaleR2 = mean(Cu_shale1to6(4:6));
Cu_shaleR3 = mean(Cu_shale1to6(7:9));
Cu_shaleR4 = mean(Cu_shale1to6(10:12));
Cu_shaleR5 = mean(Cu_shale1to6(13:15));
Cu_shaleR6 = mean(Cu_shale1to6(16:18));
meanCu_shale = [Cu_shaleR1,Cu_shaleR2,Cu_shaleR3,Cu_shaleR4,Cu_shaleR5,Cu_shaleR6]';

Cu_shellR1 = mean(Cu_shell1to6(1:3));
Cu_shellR2 = mean(Cu_shell1to6(4:6));
Cu_shellR3 = mean(Cu_shell1to6(7:9));
Cu_shellR4 = mean(Cu_shell1to6(10:12));
Cu_shellR5 = mean(Cu_shell1to6(13:15));
Cu_shellR6 = mean(Cu_shell1to6(16:18));
meanCu_shell = [Cu_shellR1,Cu_shellR2,Cu_shellR3,Cu_shellR4,Cu_shellR5,Cu_shellR6]';


% Zn
% Extract Zn data from dataICP_Rto6
Zn_ng1to6 = table2array(dataICP_R1to6(2:91,"Zn"));
%convert ng ml-1 to mg l-1
Zn_mg1to6 = (Zn_ng1to6/1000000)*1000; % units: mg/l
Zn_Table1to6 = table(treatment_ICP,run_ICP,Zn_mg1to6,'VariableNames',{'Treatment','Run','Zn'});

% Tables Zn individual treatments
filter_idx = (treatments_ICP == '1');
Zn_clay1to6 = Zn_mg1to6(filter_idx);

filter_idx = (treatments_ICP == '2');
Zn_glass1to6 = Zn_mg1to6(filter_idx);

filter_idx = (treatments_ICP == '3');
Zn_sand1to6 = Zn_mg1to6(filter_idx);

filter_idx = (treatments_ICP == '4');
Zn_shale1to6 = Zn_mg1to6(filter_idx);

filter_idx = (treatments_ICP == '5');
Zn_shell1to6 = Zn_mg1to6(filter_idx);

% Means Zn – treatments 
Means_Zn_All1to6 = groupsummary(Zn_Table1to6,'Treatment','mean','Zn');
Means_Zn_All1to6.Properties.VariableNames{'mean_Zn'} = 'Zn';
Means_Zn_Runs1to6 = groupsummary(Zn_Table1to6,["Treatment","Run"],"mean","Zn");
Means_Zn_Runs1to6.Properties.VariableNames{'mean_Zn'} = 'Zn';

% Tables Zn individual runs within treatments 
Zn_clayR1 = Zn_clay1to6(1:3);
Zn_clayR2 = mean(Zn_clay1to6(4:6));
Zn_clayR3 = mean(Zn_clay1to6(7:9));
Zn_clayR4 = mean(Zn_clay1to6(10:12));
Zn_clayR5 = mean(Zn_clay1to6(13:15));
Zn_clayR6 = mean(Zn_clay1to6(16:18));
meanZn_clay = [11.7133,Zn_clayR2,Zn_clayR3,Zn_clayR4,Zn_clayR5,Zn_clayR6]';

Zn_glassR1 = mean(Zn_glass1to6(1:3));
Zn_glassR2 = mean(Zn_glass1to6(4:6));
Zn_glassR3 = mean(Zn_glass1to6(7:9));
Zn_glassR4 = mean(Zn_glass1to6(10:12));
Zn_glassR5 = mean(Zn_glass1to6(13:15));
Zn_glassR6 = mean(Zn_glass1to6(16:18));
meanZn_glass = [Zn_glassR1,Zn_glassR2,Zn_glassR3,Zn_glassR4,Zn_glassR5,Zn_glassR6]';

Zn_sandR1 = mean(Zn_sand1to6(1:3));
Zn_sandR2 = mean(Zn_sand1to6(4:6));
Zn_sandR3 = mean(Zn_sand1to6(7:9));
Zn_sandR4 = mean(Zn_sand1to6(10:12));
Zn_sandR5 = mean(Zn_sand1to6(13:15));
Zn_sandR6 = mean(Zn_sand1to6(16:18));
meanZn_sand = [Zn_sandR1,Zn_sandR2,Zn_sandR3,Zn_sandR4,Zn_sandR5,Zn_sandR6]';

Zn_shaleR1 = mean(Zn_shale1to6(1:3));
Zn_shaleR2 = mean(Zn_shale1to6(4:6));
Zn_shaleR3 = mean(Zn_shale1to6(7:9));
Zn_shaleR4 = mean(Zn_shale1to6(10:12));
Zn_shaleR5 = mean(Zn_shale1to6(13:15));
Zn_shaleR6 = mean(Zn_shale1to6(16:18));
meanZn_shale = [Zn_shaleR1,Zn_shaleR2,Zn_shaleR3,Zn_shaleR4,Zn_shaleR5,Zn_shaleR6]';

Zn_shellR1 = mean(Zn_shell1to6(1:3));
Zn_shellR2 = mean(Zn_shell1to6(4:6));
Zn_shellR3 = mean(Zn_shell1to6(7:9));
Zn_shellR4 = mean(Zn_shell1to6(10:12));
Zn_shellR5 = mean(Zn_shell1to6(13:15));
Zn_shellR6 = mean(Zn_shell1to6(16:18));
meanZn_shell = [Zn_shellR1,Zn_shellR2,Zn_shellR3,Zn_shellR4,Zn_shellR5,Zn_shellR6]';

% Pb
Pb1to6 = table2array(dataICP_R1to6(2:91,"Pb"));

% recalculate all values under detection limit (i.e. "<") to 65% DL
for i = 1:length(Pb1to6)
    if contains(Pb1to6(i),"<")
        val = str2double(erase(Pb1to6(i),"<"));
        Pb1to6(i) = string(val*0.65);
    end    
end

%convert ng ml-1 to mg l-1
Pb_ng1to6 = str2double(Pb1to6);
Pb_mg1to6 = (Pb_ng1to6/1000000)*1000; % units: mg/l
Pb_Table1to6 = table(treatment_ICP,run_ICP,Pb_mg1to6,'VariableNames',{'Treatment','Run','Pb'});

% Tables Pb individual treatments
filter_idx = (treatments_ICP == '1');
Pb_clay1to6 = Pb_mg1to6(filter_idx);

filter_idx = (treatments_ICP == '2');
Pb_glass1to6 = Pb_mg1to6(filter_idx);

filter_idx = (treatments_ICP == '3');
Pb_sand1to6 = Pb_mg1to6(filter_idx);

filter_idx = (treatments_ICP == '4');
Pb_shale1to6 = Pb_mg1to6(filter_idx);

filter_idx = (treatments_ICP == '5');
Pb_shell1to6 = Pb_mg1to6(filter_idx);

% Means Pb – treatments 
Means_Pb_All1to6 = groupsummary(Pb_Table1to6,'Treatment','mean','Pb');
Means_Pb_All1to6.Properties.VariableNames{'mean_Pb'} = 'Pb';
Means_Pb_Runs1to6 = groupsummary(Pb_Table1to6,["Treatment","Run"],"mean","Pb");
Means_Pb_Runs1to6.Properties.VariableNames{'mean_Pb'} = 'Pb';

% Tables Pb individual treatments within runs
Pb_clayR1 = Pb_clay1to6(1:3);
Pb_clayR2 = mean(Pb_clay1to6(4:6));
Pb_clayR3 = mean(Pb_clay1to6(7:9));
Pb_clayR4 = mean(Pb_clay1to6(10:12));
Pb_clayR5 = mean(Pb_clay1to6(13:15));
Pb_clayR6 = mean(Pb_clay1to6(16:18));
meanPb_clay = [0.0446,Pb_clayR2,Pb_clayR3,Pb_clayR4,Pb_clayR5,Pb_clayR6]';

Pb_glassR1 = mean(Pb_glass1to6(1:3));
Pb_glassR2 = mean(Pb_glass1to6(4:6));
Pb_glassR3 = mean(Pb_glass1to6(7:9));
Pb_glassR4 = mean(Pb_glass1to6(10:12));
Pb_glassR5 = mean(Pb_glass1to6(13:15));
Pb_glassR6 = mean(Pb_glass1to6(16:18));
meanPb_glass = [Pb_glassR1,Pb_glassR2,Pb_glassR3,Pb_glassR4,Pb_glassR5,Pb_glassR6]';

Pb_sandR1 = mean(Pb_sand1to6(1:3));
Pb_sandR2 = mean(Pb_sand1to6(4:6));
Pb_sandR3 = mean(Pb_sand1to6(7:9));
Pb_sandR4 = mean(Pb_sand1to6(10:12));
Pb_sandR5 = mean(Pb_sand1to6(13:15));
Pb_sandR6 = mean(Pb_sand1to6(16:18));
meanPb_sand = [Pb_sandR1,Pb_sandR2,Pb_sandR3,Pb_sandR4,Pb_sandR5,Pb_sandR6]';

Pb_shaleR1 = mean(Pb_shale1to6(1:3));
Pb_shaleR2 = mean(Pb_shale1to6(4:6));
Pb_shaleR3 = mean(Pb_shale1to6(7:9));
Pb_shaleR4 = mean(Pb_shale1to6(10:12));
Pb_shaleR5 = mean(Pb_shale1to6(13:15));
Pb_shaleR6 = mean(Pb_shale1to6(16:18));
meanPb_shale = [Pb_shaleR1,Pb_shaleR2,Pb_shaleR3,Pb_shaleR4,Pb_shaleR5,Pb_shaleR6]';

Pb_shellR1 = mean(Pb_shell1to6(1:3));
Pb_shellR2 = mean(Pb_shell1to6(4:6));
Pb_shellR3 = mean(Pb_shell1to6(7:9));
Pb_shellR4 = mean(Pb_shell1to6(10:12));
Pb_shellR5 = mean(Pb_shell1to6(13:15));
Pb_shellR6 = mean(Pb_shell1to6(16:18));
meanPb_shell = [Pb_shellR1,Pb_shellR2,Pb_shellR3,Pb_shellR4,Pb_shellR5,Pb_shellR6]';


% Create Matrices for each pollutant 
runs_unique = [1,2,3,4,5,6];

NO3Matrix = NaN(length(runs_unique), length(treatments_unique));% Initialize a matrix to store NO3 concentrations
for i = 1:length(treatments_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(runs_unique)
        % Find the corresponding row in the table
        idx = (Means_NO3_Runs1to6.Treatment == treatments_unique(i)) & (Means_NO3_Runs1to6.Run == runs_unique(j));
        % Extract the mean infiltration rate for this treatment and run
        if any(idx)
            NO3Matrix(j, i) = Means_NO3_Runs1to6.NO3(idx);
        end
    end
end

NH4Matrix = NaN(length(runs_unique), length(treatments_unique));% Initialize a matrix to store NH4 concentration 
for i = 1:length(treatments_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(runs_unique)
        % Find the corresponding row in the table
        idx = (Means_NH4_Runs1to6.Treatment == treatments_unique(i)) & (Means_NH4_Runs1to6.Run == runs_unique(j));
        % Extract the mean infiltration rate for this treatment and run
        if any(idx)
            NH4Matrix(j, i) = Means_NH4_Runs1to6.NH4(idx);
        end
    end
end

PMatrix = NaN(length(runs_unique), length(treatments_unique));% Initialize a matrix to store infiltration rates
for i = 1:length(treatments_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(runs_unique)
        % Find the corresponding row in the table
        idx = (Means_P_Runs1to6.Treatment == treatments_unique(i)) & (Means_P_Runs1to6.Run == runs_unique(j));
        % Extract the mean infiltration rate for this treatment and run
        if any(idx)
            PMatrix(j, i) = Means_P_Runs1to6.P(idx);
        end
    end
end

CuMatrix = NaN(length(runs_unique), length(treatments_unique));% Initialize a matrix to store infiltration rates
for i = 1:length(treatments_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(runs_unique)
        % Find the corresponding row in the table
        idx = (Means_Cu_Runs1to6.Treatment == treatments_unique(i)) & (Means_Cu_Runs1to6.Run == runs_unique(j));
        % Extract the mean infiltration rate for this treatment and run
        if any(idx)
            CuMatrix(j, i) = Means_Cu_Runs1to6.Cu(idx);
        end
    end
end

ZnMatrix = NaN(length(runs_unique), length(treatments_unique));% Initialize a matrix to store infiltration rates
for i = 1:length(treatments_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(runs_unique)
        % Find the corresponding row in the table
        idx = (Means_Zn_Runs1to6.Treatment == treatments_unique(i)) & (Means_Zn_Runs1to6.Run == runs_unique(j));
        % Extract the mean infiltration rate for this treatment and run
        if any(idx)
            ZnMatrix(j, i) = Means_Zn_Runs1to6.Zn(idx);
        end
    end
end

PbMatrix = NaN(length(runs_unique), length(treatments_unique));% Initialize a matrix to store infiltration rates
for i = 1:length(treatments_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(runs_unique)
        % Find the corresponding row in the table
        idx = (Means_Pb_Runs1to6.Treatment == treatments_unique(i)) & (Means_Pb_Runs1to6.Run == runs_unique(j));
        % Extract the mean infiltration rate for this treatment and run
        if any(idx)
            PbMatrix(j, i) = Means_Pb_Runs1to6.Pb(idx);
        end
    end
end



% Create stacked plots
% Clay
ClayTable = array2table([NO3Matrix(:,1),NH4Matrix(:,1),PMatrix(:,1),CuMatrix(:,1),ZnMatrix(:,1),PbMatrix(:,1)]);
header = {'NO3','NH4','P','Cu','Zn','Pb'};
ClayTable.Properties.VariableNames = header;
stackedplot(ClayTable);

tiledlayout(6,1);
nexttile
plot(ClayTable,'NO3');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('NO_{3}-N clay (mg L^{-1})');
ytickformat('%.3f');
%yline(2.53,'LineWidth',2);
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,9,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(3,3.50,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(4,1.50,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(5,1.25,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(5.90,1.25,0.75,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(ClayTable,'NH4');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('NH_{4}-N clay (mg L^{-1})');
ytickformat('%.3f');
%yline(4.61,'LineWidth',2);
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.125,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(3,0.125,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.125,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,0.125,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.90,0.125,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(ClayTable,'P');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('P clay (mg L^{-1})');
ytickformat('%.3f');
%yline(0.640,'LineWidth',2);
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.1,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,0.1,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.1,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,0.1,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.90,0.1,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(ClayTable,'Cu');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('Cu clay (mg L^{-1})');
ytickformat('%.3f');
%yline(0.08,'LineWidth',2);
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.015,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,0.015,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.0125,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,0.0125,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.90,0.0125,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(ClayTable,'Zn');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('Zn clay (mg L^{-1})');
ytickformat('%.3f');
%yline(0.48,'LineWidth',2);
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,2,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,1.5,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,1.5,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,1.5,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.90,1.5,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(ClayTable,'Pb');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('Pb clay (mg L^{-1})');
ytickformat('%.3f');
%yline(0.02,'LineWidth',2);
xlabel('Run Number');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'ab';
text(2,0.01,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,0.01,txt,'FontSize',12,'FontWeight','bold');
txt = 'abc';
text(4,0.01,txt,'FontSize',12,'FontWeight','bold');
txt = 'bc';
text(5,0.01,txt,'FontSize',12,'FontWeight','bold');
txt = 'c';
text(5.90,0.01,txt,'FontSize',12,'FontWeight','bold');

% Glass
GlassTable = array2table([NO3Matrix(:,2),NH4Matrix(:,2),PMatrix(:,2),CuMatrix(:,2),ZnMatrix(:,2),PbMatrix(:,2)]);
header = {'NO3','NH4','P','Cu','Zn','Pb'};
GlassTable.Properties.VariableNames = header;
stackedplot(GlassTable);

tiledlayout(6,1);
nexttile
plot(GlassTable,'NO3');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('NO_{3}-N glass (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,4,txt,'FontSize',12,'FontWeight','bold');
txt = 'ab';
text(3,2.5,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(4,1.25,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(5,1,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(5.90,1.25,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(GlassTable,'NH4');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('NH_{4}-N glass (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.5,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,2.0,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.5,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,0.5,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.90,0.5,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(GlassTable,'P');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('P glass (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.13,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,0.12,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.125,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,0.18,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.9,0.23,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(GlassTable,'Cu');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('Cu glass (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.014,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,0.015,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.015,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,0.015,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.90,0.013,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(GlassTable,'Zn');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('Zn glass (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,1.25,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,0.75,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.75,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,0.75,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.90,0.75,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(GlassTable,'Pb');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('Pb glass (mg L^{-1})');
ytickformat('%.3f');
xlabel('Run Number');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.0011,txt,'FontSize',12,'FontWeight','bold');
txt = 'abc';
text(3,0.0008,txt,'FontSize',12,'FontWeight','bold');
txt = 'ab';
text(4,0.001,txt,'FontSize',12,'FontWeight','bold');
txt = 'bc';
text(5,0.00065,txt,'FontSize',12,'FontWeight','bold');
txt = 'c';
text(5.90,0.0007,txt,'FontSize',12,'FontWeight','bold');

% Sand
SandTable = array2table([NO3Matrix(:,3),NH4Matrix(:,3),PMatrix(:,3),CuMatrix(:,3),ZnMatrix(:,3),PbMatrix(:,3)]);
header = {'NO3','NH4','P','Cu','Zn','Pb'};
SandTable.Properties.VariableNames = header;
stackedplot(SandTable);

tiledlayout(6,1);
nexttile
plot(SandTable,'NO3');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('NO_{3}-N sand (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,2,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,1,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(4,0.75,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(5,0.75,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(5.90,1,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(SandTable,'NH4');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('NH_{4}-N sand (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.05,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(3,0.125,txt,'FontSize',12,'FontWeight','bold');
txt = 'ab';
text(4,0.075,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,0.05,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.90,0.05,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(SandTable,'P');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('P sand (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,2.70,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,2.25,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,2.45,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,2.45,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.9,2.40,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(SandTable,'Cu');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('Cu sand (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.009,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,0.008,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.009,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,0.011,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.90,0.012,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(SandTable,'Zn');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('Zn sand (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.45,txt,'FontSize',12,'FontWeight','bold');
txt = 'ab';
text(3,0.45,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(4,0.45,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(5,0.45,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(5.90,0.45,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(SandTable,'Pb');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('Pb sand (mg L^{-1})');
ytickformat('%.3f');
xlabel('Run Number');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.0009,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,0.001,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.0014,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(5,0.0035,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(5.90,0.0045,txt,'FontSize',12,'FontWeight','bold');

% Shale
ShaleTable = array2table([NO3Matrix(:,4),NH4Matrix(:,4),PMatrix(:,4),CuMatrix(:,4),ZnMatrix(:,4),PbMatrix(:,4)]);
header = {'NO3','NH4','P','Cu','Zn','Pb'};
ShaleTable.Properties.VariableNames = header;
stackedplot(ShaleTable);

tiledlayout(6,1);
nexttile
plot(ShaleTable,'NO3');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('NO_{3}-N shale (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,4.2,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(3,2.4,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(4,1,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(5,0.75,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(5.90,0.75,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(ShaleTable,'NH4');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('NH_{4}-N shale (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.05,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,0.125,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.06,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,0.05,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.90,0.05,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(ShaleTable,'P');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('P shale (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.03,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,0.02,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.02,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,0.04,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.9,0.08,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(ShaleTable,'Cu');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('Cu shale (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.009,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,0.009,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.0075,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,0.0075,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.90,0.008,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(ShaleTable,'Zn');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('Zn shale (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,1,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,0.75,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.5,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,0.5,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.90,0.5,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(ShaleTable,'Pb');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('Pb shale (mg L^{-1})');
ytickformat('%.3f');
xlabel('Run Number');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.002,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,0.002,txt,'FontSize',12,'FontWeight','bold');
txt = 'ab';
text(4,0.002,txt,'FontSize',12,'FontWeight','bold');
txt = 'ab';
text(5,0.002,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(5.90,0.002,txt,'FontSize',12,'FontWeight','bold');

% Shell 
ShellTable = array2table([NO3Matrix(:,5),NH4Matrix(:,5),PMatrix(:,5),CuMatrix(:,5),ZnMatrix(:,5),PbMatrix(:,5)]);
header = {'NO3','NH4','P','Cu','Zn','Pb'};
ShellTable.Properties.VariableNames = header;
stackedplot(ShellTable);

tiledlayout(6,1);
nexttile
plot(ShellTable,'NO3');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('NO_{3}-N shell (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'ab';
text(2,5.3,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,3,txt,'FontSize',12,'FontWeight','bold');
txt = 'ab';
text(4,2,txt,'FontSize',12,'FontWeight','bold');
txt = 'ab';
text(5,1,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(5.90,1,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(ShellTable,'NH4');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('NH_{4}-N shell (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'bc';
text(2,0.5,txt,'FontSize',12,'FontWeight','bold');
txt = 'ab';
text(3,0.5,txt,'FontSize',12,'FontWeight','bold');
txt = 'c';
text(4,0.5,txt,'FontSize',12,'FontWeight','bold');
txt = 'bc';
text(5,0.6,txt,'FontSize',12,'FontWeight','bold');
txt = 'bc';
text(5.8,0.5,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(ShellTable,'P');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('P shell (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.06,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,0.06,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.06,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,0.06,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.9,0.13,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(ShellTable,'Cu');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('Cu shell (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.0175,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,0.02,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.0175,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,0.0175,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.90,0.015,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(ShellTable,'Zn');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('Zn shell (mg L^{-1})');
ytickformat('%.3f');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,1.75,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,1,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.75,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5,0.75,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(5.90,0.75,txt,'FontSize',12,'FontWeight','bold');

nexttile
plot(ShellTable,'Pb');
xticks([1 2 3 4 5 6]);
xticklabels({'1','2', '3', '4', '5', '6'});
ylabel('Pb shell (mg L^{-1})');
ytickformat('%.3f');
xlabel('Run Number');
set(gca,'FontSize',9,'FontWeight','bold');
txt = 'a';
text(2,0.001,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(3,0.00115,txt,'FontSize',12,'FontWeight','bold');
txt = 'a';
text(4,0.00115,txt,'FontSize',12,'FontWeight','bold');
txt = 'ab';
text(5,0.00075,txt,'FontSize',12,'FontWeight','bold');
txt = 'b';
text(5.90,0.0007,txt,'FontSize',12,'FontWeight','bold');
