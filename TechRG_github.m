%% Matlab script for the data analysis for technosolraingardens
% Last edited on 1/14/26 by DMW
% Part I calculates variables and treatment means
% Part II tests for normal distribution between treatments
% Part III performs statistical analysis between treatments as presented in manuscript tables
% Part IV tests for normal distribution within runs for each treatment
% Part V statistical analysis within runs for each treatment as presented in Figs 2 through 6
% Part VI calculates treatment means within each run for each pollutant and creates the time series figures (runs 1 to 6) for each treatment for runs 1 through 6 as presented in manuscript figures 


%% Part I
% Import data files
filename = 'Infiltration_Tests_clean.xlsx';
sheetname = 'Matlab2to6';
dataInf = readtable(filename,'Sheet', sheetname);

%read in variables
treatment= table2array(dataInf(2:76,"Treatment"));
treatments = categorical(treatment); % change from cell to categorical
treatments_run = treatment(1:15,:);
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
P_mg = (P_ng/1000000)*1000; % units: mg/l
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
Zn_ng = table2array(dataICP(2:76,"Zn"));
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
% Test for normal distribution between treatments
% use fitlm for 2+ way ANOVA to calculate residuals; test residuals for normal distribution

% time to drain 
% fitlm for T2D all runs combined  
T2D_Table.Treatment = categorical(T2D_Table.Treatment);
T2D_Table.Run = categorical(T2D_Table.Run);
mdl_T2D = fitlm(T2D_Table, 'T2D ~ Treatment + Run');
resid_T2D = mdl_T2D.Residuals.Raw;
qqplot(resid_T2D);

% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swT2D, pvalue_swT2D, W_swT2D] = swtest(resid_T2D, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed


% infiltration rate 
% fitlm for IR all runs combined  
IR_Table.Treatment = categorical(IR_Table.Treatment);
IR_Table.Run = categorical(IR_Table.Run);
mdl_IR = fitlm(IR_Table, 'InfiltrationRate ~ Treatment + Run');
resid_IR = mdl_IR.Residuals.Raw;
qqplot(resid_IR);

% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swIR, pvalue_swIR, W_swIR] = swtest(resid_IR, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% drainage volume @ 15 minutes
% fitlm for DV15min all runs combined  
DV15min_Table.Treatment = categorical(DV15min_Table.Treatment);
DV15min_Table.Run = categorical(DV15min_Table.Run);
mdl_DV15min = fitlm(DV15min_Table, 'Vol15min ~ Treatment + Run');
resid_DV15min = mdl_DV15min.Residuals.Raw;
qqplot(resid_DV15min);

% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swDV15min, pvalue_swDV15min, W_swDV15min] = swtest(resid_DV15min, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% drainage volume @ 1 hour
DV1hr_Table.Treatment = categorical(DV1hr_Table.Treatment);
DV1hr_Table.Run = categorical(DV1hr_Table.Run);
mdl_DV1hr = fitlm(DV1hr_Table, 'Vol1hr ~ Treatment + Run');
resid_DV1hr = mdl_DV1hr.Residuals.Raw;
qqplot(resid_DV1hr);

% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swDV1hr, pvalue_swDV1hr, W_swDV1hr] = swtest(resid_DV1hr, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% drainage volume @ 24 hour
DV24hr_Table.Treatment = categorical(DV24hr_Table.Treatment);
DV24hr_Table.Run = categorical(DV24hr_Table.Run);
mdl_DV24hr = fitlm(DV24hr_Table, 'Vol24hr ~ Treatment + Run');
resid_DV24hr = mdl_DV24hr.Residuals.Raw;
qqplot(resid_DV24hr);

% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swDV24hr, pvalue_swDV24hr, W_swDV24hr] = swtest(resid_DV24hr, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% Water Retained in the Column
WR_Table.Treatment = categorical(WR_Table.Treatment);
WR_Table.Run = categorical(WR_Table.Run);
mdl_WR = fitlm(WR_Table, 'WR ~ Treatment + Run');
resid_WR = mdl_WR.Residuals.Raw;
qqplot(resid_WR);

% shapiro wilk test for normal distribution all runs combined
%swtest.m from matlab exchange 
[H_swWR, pvalue_swWR, W_swWR] = swtest(resid_WR, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% [NO3-N]
NO3_Table.Treatment = categorical(NO3_Table.Treatment);
NO3_Table.Run = categorical(NO3_Table.Run);
mdl_NO3 = fitlm(NO3_Table, 'NO3 ~ Treatment + Run');
resid_NO3 = mdl_NO3.Residuals.Raw;
qqplot(resid_NO3);

% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swNO3, pvalue_swNO3, W_swNO3] = swtest(resid_NO3, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% [NH4-N]
NH4_Table.Treatment = categorical(NH4_Table.Treatment);
NH4_Table.Run = categorical(NH4_Table.Run);
mdl_NH4 = fitlm(NH4_Table, 'NH4 ~ Treatment + Run');
resid_NH4 = mdl_NH4.Residuals.Raw;
qqplot(resid_NH4);

% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swNH4, pvalue_swNH4, W_swNH4] = swtest(resid_NH4, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% [P]
P_Table.Treatment = categorical(P_Table.Treatment);
P_Table.Run = categorical(P_Table.Run);
mdl_P = fitlm(P_Table, 'P ~ Treatment + Run');
resid_P = mdl_P.Residuals.Raw;
qqplot(resid_P);

% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swP, pvalue_swP, W_swP] = swtest(resid_P, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% [Cu]
Cu_Table.Treatment = categorical(Cu_Table.Treatment);
Cu_Table.Run = categorical(Cu_Table.Run);
mdl_Cu = fitlm(Cu_Table, 'Cu ~ Treatment + Run');
resid_Cu = mdl_Cu.Residuals.Raw;
qqplot(resid_Cu);

% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swCu, pvalue_swCu, W_swCu] = swtest(resid_Cu, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% [Zn]
Zn_Table.Treatment = categorical(Zn_Table.Treatment);
Zn_Table.Run = categorical(Zn_Table.Run);
mdl_Zn = fitlm(Zn_Table, 'Zn ~ Treatment + Run');
resid_Zn = mdl_Zn.Residuals.Raw;
qqplot(resid_Zn);

% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swZn, pvalue_swZn, W_swZn] = swtest(resid_Zn, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% [Pb]
Pb_Table.Treatment = categorical(Pb_Table.Treatment);
Pb_Table.Run = categorical(Pb_Table.Run);
mdl_Pb = fitlm(Pb_Table, 'Pb ~ Treatment + Run');
resid_Pb = mdl_Pb.Residuals.Raw;
qqplot(resid_Pb);

% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swPb, pvalue_swPb, W_swPb] = swtest(resid_Pb, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% tempe cell: volumetric water content @ -33 kPa --> field capacity 
VWC33_Table.Treatment = categorical(VWC33_Table.Treatment);
VWC33_Table.Rep = categorical (VWC33_Table.Rep);
mdl_VWC33 = fitlm(VWC33_Table, 'VWC33 ~ Treatment + Rep');
resid_VWC33 = mdl_VWC33.Residuals.Raw;
qqplot(resid_VWC33);

% SW test for normal distribution 
[H_VWC33, pvalue_VWC33, W_VWC33] = swtest(resid_VWC33, 0.05);
% h = 0, residuals are normally distributed

% pressure plate: volumetric water content @ -1500 kPa --> permanent wilting point
VWC1500_Table.Treatment = categorical(VWC1500_Table.Treatment);
VWC1500_Table.Rep = categorical (VWC1500_Table.Rep);
mdl_VWC1500 = fitlm(VWC1500_Table, 'VWC1500 ~ Treatment + Rep');
resid_VWC1500 = mdl_VWC1500.Residuals.Raw;
qqplot(resid_VWC1500);

% SW test for normal distribution 
[H_VWC1500, pvalue_VWC1500, W_VWC1500] = swtest(resid_VWC1500, 0.05);
% h = 0, residuals are normally distributed

% available water content 
AWC_Table.Treatment = categorical(AWC_Table.Treatment);
AWC_Table.Rep = categorical (AWC_Table.Rep);
mdl_AWC = fitlm(AWC_Table, 'AWC ~ Treatment + Rep');
resid_AWC = mdl_AWC.Residuals.Raw;
qqplot(resid_AWC);

% SW test for normal distribution 
[H_AWC, pvalue_AWC, W_AWC] = swtest(resid_AWC, 0.05);
% h = 0, residuals are normally distributed

% Ksat cm hr-1
ksat_Table.Treatment = categorical(ksat_Table.Treatment);
ksat_Table.Rep = categorical (ksat_Table.Rep);
mdl_ksat = fitlm(ksat_Table, 'Ksat ~ Treatment + Rep');
resid_ksat = mdl_ksat.Residuals.Raw;
qqplot(resid_ksat);

% SW test for normal distribution 
[H_ksat, pvalue_ksat, W_ksat] = swtest(resid_ksat, 0.05);
% h = 1, residuals are not normally distributed

% bulk density of packed cores
BD_Table.Treatment = categorical(BD_Table.Treatment);
BD_Table.Rep = categorical (BD_Table.Rep);
mdl_BD = fitlm(BD_Table, 'BulkDensity ~ Treatment + Rep');
resid_BD = mdl_BD.Residuals.Raw;
qqplot(resid_BD);

% SW test for normal distribution 
[H_BD, pvalue_BD, W_BD] = swtest(resid_BD, 0.05);
% h = 0, residuals are normally distributed

% effluent pH
pH_Table.Treatment = categorical(pH_Table.Treatment);
pH_Table.Run = categorical (pH_Table.Run);
mdl_pH = fitlm(pH_Table, 'pH ~ Treatment + Run');
resid_pH = mdl_pH.Residuals.Raw;
qqplot(resid_pH);

% SW test for normal distribution 
[H_pH, pvalue_pH, W_pH] = swtest(resid_pH, 0.05);
% h = 1, residuals are not normally distributed

%% Part III 
% Statistical analysis for difference between treatment means
% kruskal wallis test for non normally distributed data, analagous to one-way ANOVA
% aligned rank transformation (ART) ANOVA for multifactorial interactions on univariate data 

% time for drainage to start -- ART ANOVA
art_T2D = ([treatment, runs_onetofive, T2D]);
raligned_T2D = Aligned_Rank_Transform(art_T2D);
% separate raligned columns
raligned_T2D1 = raligned_T2D(:,1);
raligned_T2D2 = raligned_T2D(:,2);
raligned_T2D3 = raligned_T2D(:,3);

% raligned_T2D1 -- Treatment effect
[pT2D_anova1, tblT2D_anova1, statsT2D_anova1] = anovan(raligned_T2D1,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsT2D_anova1);
% raligned_T2D2 -- Run effect
[pT2D_anova2, tblT2D_anova2, statsT2D_anova2] = anovan(raligned_T2D2,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsT2D_anova2,'Dimension',2);
% raligned_T2D3 -- Treatment x Run effect
[pT2D_anova3, tblT2D_anova3, statsT2D_anova3] = anovan(raligned_T2D3,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsT2D_anova3,'Dimension',[1,2]);

% infiltration rate -- ART ANOVA
art_IR = ([treatment, runs_onetofive, IR]);
raligned_IR = Aligned_Rank_Transform(art_IR);
% separate raligned columns
raligned_IR1 = raligned_IR(:,1);
raligned_IR2 = raligned_IR(:,2);
raligned_IR3 = raligned_IR(:,3);

% raligned_IR1
[pIR_anova1, tblIR_anova1, statsIR_anova1] = anovan(raligned_IR1,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsIR_anova1);
% raligned_IR2
[pIR_anova2, tblIR_anova2, statsIR_anova2] = anovan(raligned_IR2,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsIR_anova2,'Dimension',2);
% raligned_IR3
[pIR_anova3, tblIR_anova3, statsIR_anova3] = anovan(raligned_IR3,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsIR_anova3,'Dimension',[1,2]);

% drainage volume @ 15 min -- ART ANOVA
art_DV15min = ([treatment, runs_onetofive, DV15min]);
raligned_DV15min = Aligned_Rank_Transform(art_DV15min);
% separate raligned columns
raligned_DV15min1 = raligned_DV15min(:,1);
raligned_DV15min2 = raligned_DV15min(:,2);
raligned_DV15min3 = raligned_DV15min(:,3);

% raligned_DV15min1
[pDV15min_anova1, tblDV15min_anova1, statsDV15min_anova1] = anovan(raligned_DV15min1,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsDV15min_anova1);
% raligned_DV15min2
[pDV15min_anova2, tblDV15min_anova2, statsDV15min_anova2] = anovan(raligned_DV15min2,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsDV15min_anova2,'Dimension',2);
% raligned_DV15min3
[pDV15min_anova3, tblDV15min_anova3, statsDV15min_anova3] = anovan(raligned_DV15min3,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsDV15min_anova3,'Dimension',[1,2]);

% drainage volume @ 1 hour -- ART ANOVA
art_DV1hr = ([treatment, runs_onetofive, DV1hr]);
raligned_DV1hr = Aligned_Rank_Transform(art_DV1hr);
% separate raligned columns
raligned_DV1hr1 = raligned_DV1hr(:,1);
raligned_DV1hr2 = raligned_DV1hr(:,2);
raligned_DV1hr3 = raligned_DV1hr(:,3);

% raligned_DV1hr1
[pDV1hr_anova1, tblDV1hr_anova1, statsDV1hr_anova1] = anovan(raligned_DV1hr1,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsDV1hr_anova1);
% raligned_DV1hr2
[pDV1hr_anova2, tblDV1hr_anova2, statsDV1hr_anova2] = anovan(raligned_DV1hr2,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsDV1hr_anova2,'Dimension',2);
% raligned_DV1hr3
[pDV1hr_anova3, tblDV1hr_anova3, statsDV1hr_anova3] = anovan(raligned_DV1hr3,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsDV1hr_anova3,'Dimension',[1,2]);

% drainage volume @ 24 hours -- ART ANOVA
art_DV24hr = ([treatment, runs_onetofive, DV24hr]);
raligned_DV24hr = Aligned_Rank_Transform(art_DV24hr);
% separate raligned columns
raligned_DV24hr1 = raligned_DV24hr(:,1);
raligned_DV24hr2 = raligned_DV24hr(:,2);
raligned_DV24hr3 = raligned_DV24hr(:,3);

% raligned_DV24hr1
[pDV24hr_anova1, tblDV24hr_anova1, statsDV24hr_anova1] = anovan(raligned_DV24hr1,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsDV24hr_anova1);
% raligned_DV24hr2
[pDV24hr_anova2, tblDV24hr_anova2, statsDV24hr_anova2] = anovan(raligned_DV24hr2,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsDV24hr_anova2,'Dimension',2);
% raligned_DV24hr3
[pDV24hr_anova3, tblDV24hr_anova3, statsDV24hr_anova3] = anovan(raligned_DV24hr3,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsDV24hr_anova3,'Dimension',[1,2]);

% water retained in column -- ART ANOVA
art_WR = ([treatment, runs_onetofive, WR]);
raligned_WR = Aligned_Rank_Transform(art_WR);
% separate raligned columns
raligned_WR1 = raligned_WR(:,1);
raligned_WR2 = raligned_WR(:,2);
raligned_WR3 = raligned_WR(:,3);

% raligned_WR1
[pWR_anova1, tblWR_anova1, statsWR_anova1] = anovan(raligned_WR1,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsWR_anova1);
% raligned_WR2
[pWR_anova2, tblWR_anova2, statsWR_anova2] = anovan(raligned_WR2,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsWR_anova2,'Dimension',2);
% raligned_WR3
[pWR_anova3, tblWR_anova3, statsWR_anova3] = anovan(raligned_WR3,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsWR_anova3,'Dimension',[1,2]);

% NO3-N concentration and removal -- ART ANOVA
art_NO3 = ([treatment, runs_onetofive, NO3_vertical]);
raligned_NO3 = Aligned_Rank_Transform(art_NO3);
% separate raligned columns
raligned_NO31 = raligned_NO3(:,1);
raligned_NO32 = raligned_NO3(:,2);
raligned_NO33 = raligned_NO3(:,3);

% raligned_NO31
[pNO3_anova1, tblNO3_anova1, statsNO3_anova1] = anovan(raligned_NO31,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsNO3_anova1);
% raligned_NO32
[pNO3_anova2, tblNO3_anova2, statsNO3_anova2] = anovan(raligned_NO32,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsNO3_anova2,'Dimension',2);
% raligned_NO33
[pNO3_anova3, tblNO3_anova3, statsNO3_anova3] = anovan(raligned_NO33,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
[c] = multcompare(statsNO3_anova3,'Dimension',[1,2]);


% NH4-N concentration and removal -- ART ANOVA 
art_NH4 = ([treatment, runs_onetofive, NH4_vertical]);
raligned_NH4 = Aligned_Rank_Transform(art_NH4);
% separate raligned columns
raligned_NH41 = raligned_NH4(:,1);
raligned_NH42 = raligned_NH4(:,2);
raligned_NH43 = raligned_NH4(:,3);

% raligned_NH41
[pNH4_anova1, tblNH4_anova1, statsNH4_anova1] = anovan(raligned_NH41,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsNH4_anova1);
% raligned_NH42
[pNH4_anova2, tblNH4_anova2, statsNH4_anova2] = anovan(raligned_NH42,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsNH4_anova2,'Dimension',2);
% raligned_NH43
[pNH4_anova3, tblNH4_anova3, statsNH4_anova3] = anovan(raligned_NH43,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsNH4_anova3,'Dimension',[1,2]);

% P concentration and removal -- ART ANOVA
art_P = ([treatment, runs_onetofive, P_mg]);
raligned_P = Aligned_Rank_Transform(art_P);
% separate raligned columns
raligned_P1 = raligned_P(:,1);
raligned_P2 = raligned_P(:,2);
raligned_P3 = raligned_P(:,3);

% raligned_P1
[pP_anova1, tblP_anova1, statsP_anova1] = anovan(raligned_P1,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsP_anova1);
% raligned_P2
[pP_anova2, tblP_anova2, statsP_anova2] = anovan(raligned_P2,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsP_anova2,'Dimension',2);
% raligned_P3
[pP_anova3, tblP_anova3, statsP_anova3] = anovan(raligned_P3,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsP_anova3,'Dimension',[1,2]);

% Cu concentration and removal -- ART ANOVA
art_Cu = ([treatment, runs_onetofive, Cu_mg]);
raligned_Cu = Aligned_Rank_Transform(art_Cu);
% separate raligned columns
raligned_Cu1 = raligned_Cu(:,1);
raligned_Cu2 = raligned_Cu(:,2);
raligned_Cu3 = raligned_Cu(:,3);

% raligned_Cu1
[pCu_anova1, tblCu_anova1, statsCu_anova1] = anovan(raligned_Cu1,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsCu_anova1);
% raligned_Cu2
[pCu_anova2, tblCu_anova2, statsCu_anova2] = anovan(raligned_Cu2,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsCu_anova2,'Dimension',2);
% raligned_Cu3
[pCu_anova3, tblCu_anova3, statsCu_anova3] = anovan(raligned_Cu3,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsCu_anova3,'Dimension',[1,2]);

% Zn concentration and removal -- ART ANOVA
art_Zn = ([treatment, runs_onetofive, Zn_mg]);
raligned_Zn = Aligned_Rank_Transform(art_Zn);
% separate raligned columns
raligned_Zn1 = raligned_Zn(:,1);
raligned_Zn2 = raligned_Zn(:,2);
raligned_Zn3 = raligned_Zn(:,3);

% raligned_Zn1
[pZn_anova1, tblZn_anova1, statsZn_anova1] = anovan(raligned_Zn1,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsZn_anova1);
% raligned_Zn2
[pZn_anova2, tblZn_anova2, statsZn_anova2] = anovan(raligned_Zn2,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
[c] = multcompare(statsZn_anova2,'Dimension',2);
% raligned_Zn3
[pZn_anova3, tblZn_anova3, statsZn_anova3] = anovan(raligned_Zn3,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
[c]=multcompare(statsZn_anova3,'Dimension',[1,2]);

% Pb concentration and removal -- ART ANOVA
art_Pb = ([treatment, runs_onetofive, Pb_mg]);
raligned_Pb = Aligned_Rank_Transform(art_Pb);
% separate raligned columns
raligned_Pb1 = raligned_Pb(:,1);
raligned_Pb2 = raligned_Pb(:,2);
raligned_Pb3 = raligned_Pb(:,3);

% raligned_Pb1
[pPb_anova1, tblPb_anova1, statsPb_anova1] = anovan(raligned_Pb1,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsPb_anova1);
% raligned_Pb2
[pPb_anova2, tblPb_anova2, statsPb_anova2] = anovan(raligned_Pb2,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsPb_anova2,'Dimension',2);
% raligned_Pb3
[pPb_anova3, tblPb_anova3, statsPb_anova3] = anovan(raligned_Pb3,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statsPb_anova3,'Dimension',[1,2]);

% tempe cell: VWC @ -33 kPa  field capacity -- normal distribution
VWC33Matrix = NaN(length(reps_unique), length(treatments_unique));% Initialize a matrix to store NO3 concentrations

for i = 1:length(treatments_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(reps_unique)
        % Find the corresponding row in the table
        idx = (VWC33_Table.Treatment == treatments_unique(i)) & (VWC33_Table.Rep == reps_unique(j));
        % Extract the mean VWC33 for this treatment and rep
        if any(idx)
            VWC33Matrix(j, i) = VWC33_Table.VWC33(idx);
        end
    end
end

[p_VWC33, stats_VWC33, tbl_VWC33] = anova1(VWC33Matrix);
multcompare(tbl_VWC33);

% pressure plate: VWC @ -1500 kPa  permanent wilting point -- normal distribution
VWC1500Matrix = NaN(length(reps_unique), length(treatments_unique));% Initialize a matrix to store NO3 concentrations

for i = 1:length(treatments_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(reps_unique)
        % Find the corresponding row in the table
        idx = (treatments_lab == treatments_unique(i)) & (reps_lab == reps_unique(j));
        % Extract the mean VWC33 for this treatment and rep
        if any(idx)
            VWC1500Matrix(j, i) = VWC1500_Table.VWC1500(idx);
        end
    end
end

[p_VWC1500, stats_VWC1500, tbl_VWC1500] = anova1(VWC1500Matrix);
multcompare(tbl_VWC1500);


% available water capacity -- normal distribution 
AWCMatrix = NaN(length(reps_unique), length(treatments_unique));% Initialize a matrix to store NO3 concentrations

for i = 1:length(treatments_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(reps_unique)
        % Find the corresponding row in the table
        idx = (treatments_lab == treatments_unique(i)) & (reps_lab == reps_unique(j));
        % Extract the mean VWC33 for this treatment and rep
        if any(idx)
            AWCMatrix(j, i) = AWC_Table.AWC(idx);
        end
    end
end

[p_AWC, stats_AWC, tbl_AWC] = anova1(AWCMatrix);
multcompare(tbl_AWC);

% Ksat cm hr-1 -- kwtest (three repeats - three replicates - one run --> use the mean of three repeats for each soil sample)
[p_ksat,tbl_ksat,stats_ksat] = kruskalwallis(ksat_Table.Ksat,ksat_Table.Treatment); 
multcompare(stats_ksat);

% Bulk density of packed cores -- normal distribution

BDMatrix = NaN(length(reps_unique), length(treatments_unique));% Initialize a matrix to store NO3 concentrations

for i = 1:length(treatments_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(reps_unique)
        % Find the corresponding row in the table
        idx = (treatments_lab == treatments_unique(i)) & (reps_lab == reps_unique(j));
        % Extract the mean VWC33 for this treatment and rep
        if any(idx)
            BDMatrix(j, i) = BD_Table.BulkDensity(idx);
        end
    end
end

[p_BD, stats_BD, tbl_BD] = anova1(BDMatrix);
multcompare(tbl_BD);

% effluent pH -- ART ANOVA
art_pH = ([treatment, runs_onetofive, pH_column]);
raligned_pH = Aligned_Rank_Transform(art_pH);
% separate raligned columns
raligned_pH1 = raligned_pH(:,1);
raligned_pH2 = raligned_pH(:,2);
raligned_pH3 = raligned_pH(:,3);

% raligned_pH1
[ppH_anova1, tblpH_anova1, statspH_anova1] = anovan(raligned_pH1,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statspH_anova1);
% raligned_pH2
[ppH_anova2, tblpH_anova2, statspH_anova2] = anovan(raligned_pH2,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statspH_anova2,'Dimension',2);
% raligned_pH3
[ppH_anova3, tblpH_anova3, statspH_anova3] = anovan(raligned_pH3,{treatment,run},'Model','interaction','Varnames',{'Treatment','Run'});
multcompare(statspH_anova3,'Dimension',[1,2]);


%% Part IV 
% Test for normal distribution between runs (2-6) within treatments
runs2to6 = [2 2 2 3 3 3 4 4 4 5 5 5 6 6 6]'; 
reps2to6 = [1 2 3 1 2 3 1 2 3 1 2 3 1 2 3]';

% NO3-N
NO3clay = table(runs2to6,reps2to6,NO3_clay,'VariableNames',{'Run','Rep','NO3'});
NO3glass =table(runs2to6,reps2to6,NO3_glass,'VariableNames',{'Run','Rep','NO3'});
NO3sand = table(runs2to6,reps2to6,NO3_sand,'VariableNames',{'Run','Rep','NO3'});
NO3shale =table(runs2to6,reps2to6,NO3_shale,'VariableNames',{'Run','Rep','NO3'});
NO3shell = table(runs2to6,reps2to6,NO3_shell,'VariableNames',{'Run','Rep','NO3'});

% fitlm for NO3 - Clay all runs combined
mdl_clayNO3 = fitlm(NO3clay,'NO3 ~ Run');
resid_clayNO3 = mdl_clayNO3.Residuals.Raw;
%qqplot(resid_clayNO3);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swclayNO3, pvalue_swclayNO3, W_swclayNO3] = swtest(resid_clayNO3, 0.05);
% h = 0, pvalue >0.05 = residuals are normally distributed

% fitlm for NO3 - Glass all runs combined
mdl_glassNO3 = fitlm(NO3glass, 'NO3 ~ Run');
resid_glassNO3 = mdl_glassNO3.Residuals.Raw;
%qqplot(resid_glassNO3);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swglassNO3, pvalue_swglassNO3, W_swglassNO3] = swtest(resid_glassNO3, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% fitlm for NO3 - Sand all runs combined
mdl_sandNO3 = fitlm(NO3sand, 'NO3 ~ Run');
resid_sandNO3 = mdl_sandNO3.Residuals.Raw;
%qqplot(resid_sandNO3);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swsandNO3, pvalue_swsandNO3, W_swsandNO3] = swtest(resid_sandNO3, 0.05);
% h = 0, pvalue >0.05 = residuals are normally distributed

% fitlm for NO3 - Shale all runs combined
mdl_shaleNO3 = fitlm(NO3shale, 'NO3 ~ Run');
resid_shaleNO3 = mdl_shaleNO3.Residuals.Raw;
%qqplot(resid_shaleNO3);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swshaleNO3, pvalue_swshaleNO3, W_swshaleNO3] = swtest(resid_shaleNO3, 0.05);
% h = 0, pvalue >0.05 = residuals are normally distributed

% fitlm for NO3 - Shell all runs combined
mdl_shellNO3 = fitlm(NO3shell, 'NO3');
resid_shellNO3 = mdl_shellNO3.Residuals.Raw;
%qqplot(resid_shellNO3);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swshellNO3, pvalue_swshellNO3, W_swshellNO3] = swtest(resid_shellNO3, 0.05);
% h = 0, pvalue >0.05 = residuals are normally distributed

%NH4
NH4clay = table(runs2to6,reps2to6,NH4_clay,'VariableNames',{'Run','Rep','NH4'});
NH4glass =table(runs2to6,reps2to6,NH4_glass,'VariableNames',{'Run','Rep','NH4'});
NH4sand = table(runs2to6,reps2to6,NH4_sand,'VariableNames',{'Run','Rep','NH4'});
NH4shale =table(runs2to6,reps2to6,NH4_shale,'VariableNames',{'Run','Rep','NH4'});
NH4shell = table(runs2to6,reps2to6,NH4_shell,'VariableNames',{'Run','Rep','NH4'});

% fitlm for NH4 - Clay all runs combined
mdl_clayNH4 = fitlm(NH4clay, 'NH4 ~ Run');
resid_clayNH4 = mdl_clayNH4.Residuals.Raw;
%qqplot(resid_clayNH4);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swclayNH4, pvalue_swclayNH4, W_swclayNH4] = swtest(resid_clayNH4, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% fitlm for NH4 - Glass all runs combined
mdl_glassNH4 = fitlm(NH4glass, 'NH4 ~ Run');
resid_glassNH4 = mdl_glassNH4.Residuals.Raw;
%qqplot(resid_glassNH4);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swglassNH4, pvalue_swglassNH4, W_swglassNH4] = swtest(resid_glassNH4, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% fitlm for NH4 - Sand all runs combined
mdl_sandNH4 = fitlm(NH4sand, 'NH4 ~ Run');
resid_sandNH4 = mdl_sandNH4.Residuals.Raw;
%qqplot(resid_sandNH4);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swsandNH4, pvalue_swsandNH4, W_swsandNH4] = swtest(resid_sandNH4, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% fitlm for NH4 - Shale all runs combined
mdl_shaleNH4 = fitlm(NH4shale, 'NH4 ~ Run');
resid_shaleNH4 = mdl_shaleNH4.Residuals.Raw;
%qqplot(resid_shaleNH4);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swshaleNH4, pvalue_swshaleNH4, W_swshaleNH4] = swtest(resid_shaleNH4, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% fitlm for NH4 - Shell all runs combined
mdl_shellNH4 = fitlm(NH4shell, 'NH4 ~ Run');
resid_shellNH4 = mdl_shellNH4.Residuals.Raw;
%qqplot(resid_shellNH4);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swshellNH4, pvalue_swshellNH4, W_swshellNH4] = swtest(resid_shellNH4, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% P
Pclay = table(runs2to6,reps2to6,P_clay,'VariableNames',{'Run','Rep','P'});
Pglass =table(runs2to6,reps2to6,P_glass,'VariableNames',{'Run','Rep','P'});
Psand = table(runs2to6,reps2to6,P_sand,'VariableNames',{'Run','Rep','P'});
Pshale =table(runs2to6,reps2to6,P_shale,'VariableNames',{'Run','Rep','P'});
Pshell = table(runs2to6,reps2to6,P_shell,'VariableNames',{'Run','Rep','P'});

% fitlm for P - Clay all runs combined
mdl_clayP = fitlm(Pclay, 'P ~ Run');
resid_clayP = mdl_clayP.Residuals.Raw;
%qqplot(resid_clayP);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swclayP, pvalue_swclayP, W_swclayP] = swtest(resid_clayP, 0.05);
% h = 0, pvalue >0.05 = residuals are normally distributed

% fitlm for P - Glass all runs combined
mdl_glassP = fitlm(Pglass, 'P ~ Run');
resid_glassP = mdl_glassP.Residuals.Raw;
%qqplot(resid_glassP);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swglassP, pvalue_swglassP, W_swglassP] = swtest(resid_glassP, 0.05);
% h = 0, pvalue >0.05 = residuals are normally distributed

% fitlm for P - Sand all runs combined
mdl_sandP = fitlm(Psand, 'P ~ Run');
resid_sandP = mdl_sandP.Residuals.Raw;
%qqplot(resid_sandP);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swsandP, pvalue_swsandP, W_swsandP] = swtest(resid_sandP, 0.05);
% h = 0, pvalue >0.05 = residuals are normally distributed

% fitlm for P - Shale all runs combined
mdl_shaleP = fitlm(Pshale, 'P ~ Run');
resid_shaleP = mdl_shaleP.Residuals.Raw;
%qqplot(resid_shaleP);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swshaleP, pvalue_swshaleP, W_swshaleP] = swtest(resid_shaleP, 0.05);
% h = 1, pvalue >0.05 = residuals are normally distributed

% fitlm for P - Shell all runs combined
mdl_shellP = fitlm(Pshell, 'P ~ Run');
resid_shellP = mdl_shellP.Residuals.Raw;
%qqplot(resid_shellP);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swshellP, pvalue_swshellP, W_swshellP] = swtest(resid_shellP, 0.05);
% h = 0, pvalue >0.05 = residuals are normally distributed


% Cu
Cuclay = table(runs2to6,reps2to6,Cu_clay,'VariableNames',{'Run','Rep','Cu'});
Cuglass =table(runs2to6,reps2to6,Cu_glass,'VariableNames',{'Run','Rep','Cu'});
Cusand = table(runs2to6,reps2to6,Cu_sand,'VariableNames',{'Run','Rep','Cu'});
Cushale =table(runs2to6,reps2to6,Cu_shale,'VariableNames',{'Run','Rep','Cu'});
Cushell = table(runs2to6,reps2to6,Cu_shell,'VariableNames',{'Run','Rep','Cu'});

% fitlm for Cu - Clay all runs combined
mdl_clayCu = fitlm(Cuclay, 'Cu ~ Run');
resid_clayCu = mdl_clayCu.Residuals.Raw;
%qqplot(resid_clayCu);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swclayCu, pvalue_swclayCu, W_swclayCu] = swtest(resid_clayCu, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% fitlm for Cu - Glass all runs combined
mdl_glassCu = fitlm(Cuglass, 'Cu ~ Run');
resid_glassCu = mdl_glassCu.Residuals.Raw;
%qqplot(resid_glassCu);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swglassCu, pvalue_swglassCu, W_swglassCu] = swtest(resid_glassCu, 0.05);
% h = 0, pvalue >0.05 = residuals are normally distributed

% fitlm for Cu - Sand all runs combined
mdl_sandCu = fitlm(Cusand, 'Cu ~ Run');
resid_sandCu = mdl_sandCu.Residuals.Raw;
%qqplot(resid_sandCu);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swsandCu, pvalue_swsandCu, W_swsandCu] = swtest(resid_sandCu, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% fitlm for Cu - Shale all runs combined
mdl_shaleCu = fitlm(Cushale, 'Cu ~ Run');
resid_shaleCu = mdl_shaleCu.Residuals.Raw;
%qqplot(resid_shaleCu);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swshaleCu, pvalue_swshaleCu, W_swshaleCu] = swtest(resid_shaleCu, 0.05);
% h = 0, pvalue >0.05 = residuals are normally distributed

% fitlm for Cu - Shell all runs combined
mdl_shellCu = fitlm(Cushell, 'Cu ~ Run');
resid_shellCu = mdl_shellCu.Residuals.Raw;
%qqplot(resid_shellCu);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swshellCu, pvalue_swshellCu, W_swshellCu] = swtest(resid_shellCu, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% Zn
Znclay = table(runs2to6,reps2to6,Zn_clay,'VariableNames',{'Run','Rep','Zn'});
Znglass =table(runs2to6,reps2to6,Zn_glass,'VariableNames',{'Run','Rep','Zn'});
Znsand = table(runs2to6,reps2to6,Zn_sand,'VariableNames',{'Run','Rep','Zn'});
Znshale =table(runs2to6,reps2to6,Zn_shale,'VariableNames',{'Run','Rep','Zn'});
Znshell = table(runs2to6,reps2to6,Zn_shell,'VariableNames',{'Run','Rep','Zn'});

% fitlm for Zn - Clay all runs combined
mdl_clayZn = fitlm(Znclay, 'Zn ~ Run');
resid_clayZn = mdl_clayZn.Residuals.Raw;
%qqplot(resid_clayZn);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swclayZn, pvalue_swclayZn, W_swclayZn] = swtest(resid_clayZn, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% fitlm for Zn - Glass all runs combined
mdl_glassZn = fitlm(Znglass, 'Zn ~ Run');
resid_glassZn = mdl_glassZn.Residuals.Raw;
%qqplot(resid_glassZn);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swglassZn, pvalue_swglassZn, W_swglassZn] = swtest(resid_glassZn, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% fitlm for Zn - Sand all runs combined
mdl_sandZn = fitlm(Znsand, 'Zn ~ Run');
resid_sandZn = mdl_sandZn.Residuals.Raw;
%qqplot(resid_sandZn);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swsandZn, pvalue_swsandZn, W_swsandZn] = swtest(resid_sandZn, 0.05);
% h = 1, pvalue <0.05 = residuals are not normally distributed

% fitlm for Zn - Shale all runs combined
mdl_shaleZn = fitlm(Znshale, 'Zn ~ Run');
resid_shaleZn = mdl_shaleZn.Residuals.Raw;
%qqplot(resid_shaleZn);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swshaleZn, pvalue_swshaleZn, W_swshaleZn] = swtest(resid_shaleZn, 0.05);
% h = 0, pvalue >0.05 = residuals are normally distributed

% fitlm for Zn - Shell all runs combined
mdl_shellZn = fitlm(Znshell, 'Zn ~ Run');
resid_shellZn = mdl_shellZn.Residuals.Raw;
%qqplot(resid_shellZn);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swshellZn, pvalue_swshellZn, W_swshellZn] = swtest(resid_shellZn, 0.05);
% h = 0, pvalue <0.05 = residuals are normally distributed

% Pb
Pbclay = table(runs2to6,reps2to6,Pb_clay,'VariableNames',{'Run','Rep','Pb'});
Pbglass =table(runs2to6,reps2to6,Pb_glass,'VariableNames',{'Run','Rep','Pb'});
Pbsand = table(runs2to6,reps2to6,Pb_sand,'VariableNames',{'Run','Rep','Pb'});
Pbshale =table(runs2to6,reps2to6,Pb_shale,'VariableNames',{'Run','Rep','Pb'});
Pbshell = table(runs2to6,reps2to6,Pb_shell,'VariableNames',{'Run','Rep','Pb'});

% fitlm for Pb - Clay all runs combined
mdl_clayPb = fitlm(Pbclay, 'Pb ~ Run');
resid_clayPb = mdl_clayPb.Residuals.Raw;
%qqplot(resid_clayPb);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swclayPb, pvalue_swclayPb, W_swclayPb] = swtest(resid_clayPb, 0.05);
% h = 0, pvalue >0.05 = residuals are normally distributed

% fitlm for Pb - Glass all runs combined
mdl_glassPb = fitlm(Pbglass, 'Pb ~ Run');
resid_glassPb = mdl_glassPb.Residuals.Raw;
%qqplot(resid_glassPb);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swglassPb, pvalue_swglassPb, W_swglassPb] = swtest(resid_glassPb, 0.05);
% h = 0, pvalue >0.05 = residuals are normally distributed

% fitlm for Pb - Sand all runs combined
mdl_sandPb = fitlm(Pbsand, 'Pb ~ Run');
resid_sandPb = mdl_sandPb.Residuals.Raw;
%qqplot(resid_sandPb);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swsandPb, pvalue_swsandPb, W_swsandPb] = swtest(resid_sandPb, 0.05);
% h = 0, pvalue >0.05 = residuals are normally distributed

% fitlm for Pb - Shale all runs combined
mdl_shalePb = fitlm(Pbshale, 'Pb ~ Run');
resid_shalePb = mdl_shalePb.Residuals.Raw;
%qqplot(resid_shalePb);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swshalePb, pvalue_swshalePb, W_swshalePb] = swtest(resid_shalePb, 0.05);
% h = 0, pvalue >0.05 = residuals are normally distributed

% fitlm for Pb - Shell all runs combined
mdl_shellPb = fitlm(Pbshell, 'Pb ~ Run');
resid_shellPb = mdl_shellPb.Residuals.Raw;
%qqplot(resid_shellPb);
% shapiro wilk test for normal distribution all runs combined
% swtest.m from matlab exchange 
[H_swshellPb, pvalue_swshellPb, W_swshellPb] = swtest(resid_shellPb, 0.05);
% h = 0, pvalue >0.05 = residuals are normally distributed

%% Part V
% Statistical signifigance for differences between runs within treatments
% NO3 Clay - normal

% NO3 Glass - not normal
[p_kwglassNO3,tbl_kwglassNO3,stats_kwglassNO3] = kruskalwallis(NO3glass_Table.NO3,NO3glass_Table.Run); 
multcompare(stats_kwglassNO3);

% NO3 Sand - normal
SandNO3Matrix = NaN(length(reps_unique), length(runs_unique));% Initialize a matrix to store NO3 concentrations
for i = 1:length(runs_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(reps_unique)
        % Find the corresponding row in the table
        idx = (NO3sand_Table.Run == runs_unique(i)) & (NO3sand_Table.Rep == reps_unique(j));
        % Extract the mean VWC33 for this treatment and rep
        if any(idx)
            SandNO3Matrix(j, i) = NO3sand_Table.NO3(idx);
        end
    end
end
[p_sandNO3, stats_sandNO3, tbl_sandNO3] = anova1(SandNO3Matrix);
multcompare(tbl_sandNO3);

% NO3 Shale - normal
ShaleNO3Matrix = NaN(length(reps_unique), length(runs_unique));% Initialize a matrix to store NO3 concentrations
for i = 1:length(runs_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(reps_unique)
        % Find the corresponding row in the table
        idx = (NO3shale_Table.Run == runs_unique(i)) & (NO3shale_Table.Rep == reps_unique(j));
        % Extract the mean VWC33 for this treatment and rep
        if any(idx)
            ShaleNO3Matrix(j, i) = NO3shale_Table.NO3(idx);
        end
    end
end
[p_shaleNO3, stats_shaleNO3, tbl_shaleNO3] = anova1(ShaleNO3Matrix);
multcompare(tbl_shaleNO3);

% NO3 Shell - normal 
ShellNO3Matrix = NaN(length(reps_unique), length(runs_unique));% Initialize a matrix to store NO3 concentrations
for i = 1:length(runs_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(reps_unique)
        % Find the corresponding row in the table
        idx = (NO3shell_Table.Run == runs_unique(i)) & (NO3shell_Table.Rep == reps_unique(j));
        % Extract the mean VWC33 for this treatment and rep
        if any(idx)
            ShellNO3Matrix(j, i) = NO3shell_Table.NO3(idx);
        end
    end
end
[p_shellNO3, stats_shellNO3, tbl_shellNO3] = anova1(ShellNO3Matrix);
multcompare(tbl_shellNO3);

% NH4 Clay - not normal
[p_kwclayNH4,tbl_kwclayNH4,stats_kwclayNH4] = kruskalwallis(ClayNH4_Table.NH4,ClayNH4_Table.Run); 
multcompare(stats_kwclayNH4);

% NH4 Glass - not normally distributed
[p_kwglassNH4,tbl_kwglassNH4,stats_kwglassNH4] = kruskalwallis(NH4glass.NH4,NH4glass.Run); 
multcompare(stats_kwglassNH4);

% NH4 Sand - not normal
[p_kwsandNH4,tbl_kwsandNH4,stats_kwsandNH4] = kruskalwallis(NH4sand_Table.NH4,NH4sand_Table.Run); 
multcompare(stats_kwsandNH4);

% NH4 Shale - not normal
[p_kwshaleNH4,tbl_kwshaleNH4,stats_kwshaleNH4] = kruskalwallis(NH4shale_Table.NH4,NH4shale_Table.Run); 
multcompare(stats_kwshaleNH4);

% NH4 Shell - not normal
[p_kwshellNH4,tbl_kwshellNH4,stats_kwshellNH4] = kruskalwallis(NH4shell_Table.NH4,NH4shell_Table.Run); 
multcompare(stats_kwshellNH4); 

% P Clay - normal

% P Glass - normal
GlassPMatrix = NaN(length(reps_unique), length(runs_unique));% Initialize a matrix to store NO3 concentrations
for i = 1:length(runs_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(reps_unique)
        % Find the corresponding row in the table
        idx = (Pglass.Run == runs_unique(i)) & (Pglass.Rep == reps_unique(j));
        % Extract the mean VWC33 for this treatment and rep
        if any(idx)
            GlassPMatrix(j, i) = Pglass.P(idx);
        end
    end
end
[p_glassP, stats_glassP, tbl_glassP] = anova1(GlassPMatrix);
multcompare(tbl_glassP);

% P Sand - normal
SandPMatrix = NaN(length(reps_unique), length(runs_unique));% Initialize a matrix to store NO3 concentrations
for i = 1:length(runs_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(reps_unique)
        % Find the corresponding row in the table
        idx = (Psand_Table.Run == runs_unique(i)) & (Psand_Table.Rep == reps_unique(j));
        % Extract the mean VWC33 for this treatment and rep
        if any(idx)
            SandPMatrix(j, i) = Psand_Table.P(idx);
        end
    end
end
[p_sandP, stats_sandP, tbl_sandP] = anova1(SandPMatrix);
multcompare(tbl_sandP);

% P Shale - not normal

% P Shell - normal
ShellPMatrix = NaN(length(reps_unique), length(runs_unique));% Initialize a matrix to store NO3 concentrations
for i = 1:length(runs_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(reps_unique)
        % Find the corresponding row in the table
        idx = (Pshell_Table.Run == runs_unique(i)) & (Pshell_Table.Rep == reps_unique(j));
        % Extract the mean VWC33 for this treatment and rep
        if any(idx)
            ShellPMatrix(j, i) = Pshell_Table.P(idx);
        end
    end
end
[p_shellP, stats_shellP, tbl_shellP] = anova1(ShellPMatrix);
multcompare(tbl_shellP);

% Cu Clay - not normal

% Cu Glass - normal
GlassCuMatrix = NaN(length(reps_unique), length(runs_unique));% Initialize a matrix to store NO3 concentrations
for i = 1:length(runs_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(reps_unique)
        % Find the corresponding row in the table
        idx = (Cuglass_Table.Run == runs_unique(i)) & (Cuglass_Table.Rep == reps_unique(j));
        % Extract the mean VWC33 for this treatment and rep
        if any(idx)
            GlassCuMatrix(j, i) = Cuglass_Table.Cu(idx);
        end
    end
end
[p_glassCu, stats_glassCu, tbl_glassCu] = anova1(GlassCuMatrix);
multcompare(tbl_glassCu);

% Cu Sand - not normal

% Cu Shale - normal

% Cu Shell - not normal

% Zn Clay - not normal
[p_kwclayZn,tbl_kwclayZn,stats_kwclayZn] = kruskalwallis(ClayZn_Table.Zn,ClayZn_Table.Run); 
multcompare(stats_kwclayZn);

% Zn Glass - not normal
[p_kwglassZn,tbl_kwglassZn,stats_kwglassZn] = kruskalwallis(Znglass_Table.Zn,Znglass_Table.Run); 
multcompare(stats_kwglassZn);

% Zn Sand - not normal

% Zn Shale - normal
ShaleZnMatrix = NaN(length(reps_unique), length(runs_unique));% Initialize a matrix to store NO3 concentrations
for i = 1:length(runs_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(reps_unique)
        % Find the corresponding row in the table
        idx = (Znshale_Table.Run == runs_unique(i)) & (Znshale_Table.Rep == reps_unique(j));
        % Extract the mean VWC33 for this treatment and rep
        if any(idx)
            ShaleZnMatrix(j, i) = Znshale_Table.Zn(idx);
        end
    end
end
[p_shaleZn, stats_shaleZn, tbl_shaleZn] = anova1(ShaleZnMatrix);
multcompare(tbl_shaleZn);

% Zn Shell - normal

% Pb Clay - normal
ClayPbMatrix = NaN(length(reps), length(runsICPUnique));% Initialize a matrix to store NO3 concentrations
for i = 1:length(runsICPUnique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(reps)
        % Find the corresponding row in the table
        idx = (ClayPb_Table.Run == runsICPUnique(i)) & (ClayPb_Table.Rep == reps(j));
        % Extract the mean VWC33 for this treatment and rep
        if any(idx)
            ClayPbMatrix(j, i) = ClayPb_Table.Pb(idx);
        end
    end
end
[p_clayPb, stats_clayPb, tbl_clayPb] = anova1(ClayPbMatrix);
multcompare(tbl_clayPb);

% Pb Glass - normal
GlassPbMatrix = NaN(length(reps_unique), length(runs_unique));% Initialize a matrix to store NO3 concentrations
for i = 1:length(runs_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(reps_unique)
        % Find the corresponding row in the table
        idx = (Pbglass_Table.Run == runs_unique(i)) & (Pbglass_Table.Rep == reps_unique(j));
        % Extract the mean VWC33 for this treatment and rep
        if any(idx)
            GlassPbMatrix(j, i) = Pbglass_Table.Pb(idx);
        end
    end
end
[p_glassPb, stats_glassPb, tbl_glassPb] = anova1(GlassPbMatrix);
multcompare(tbl_glassPb);

% Pb Sand - normal

% Pb Shale - normal
ShalePbMatrix = NaN(length(reps_unique), length(runs_unique));% Initialize a matrix to store NO3 concentrations
for i = 1:length(runs_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(reps_unique)
        % Find the corresponding row in the table
        idx = (Pbshale_Table.Run == runs_unique(i)) & (Pbshale_Table.Rep == reps_unique(j));
        % Extract the mean VWC33 for this treatment and rep
        if any(idx)
            ShalePbMatrix(j, i) = Pbshale_Table.Pb(idx);
        end
    end
end
[p_shalePb, stats_shalePb, tbl_shalePb] = anova1(ShalePbMatrix);
multcompare(tbl_shalePb);

% Pb Shell - normal
ShellPbMatrix = NaN(length(reps_unique), length(runs_unique));% Initialize a matrix to store NO3 concentrations
for i = 1:length(runs_unique) % Loop through treatments and runs to populate the matrix
    for j = 1:length(reps_unique)
        % Find the corresponding row in the table
        idx = (Pbshell_Table.Run == runs_unique(i)) & (Pbshell_Table.Rep == reps_unique(j));
        % Extract the mean VWC33 for this treatment and rep
        if any(idx)
            ShellPbMatrix(j, i) = Pbshell_Table.Pb(idx);
        end
    end
end
[p_shellPb, stats_shellPb, tbl_shellPb] = anova1(ShellPbMatrix);
multcompare(tbl_shellPb);


%% Part VI
% Time series Pollutant concentration for each treatment Runs 1 through 6
% Figures 2-3 in manuscript
% load in pollutant results

Log_MeanNO3_RunsAll = log(MeanNO3_RunsAll);
Log_MeanNH4_RunsAll = log(MeanNH4_RunsAll);
Log_MeanP_RunsAll = log(MeanP_RunsAll);
Log_MeanCu_RunsAll = log(MeanCu_RunsAll);
Log_MeanZn_RunsAll = log(MeanZn_RunsAll);
Log_MeanPb_RunsAll = log(MeanPb_RunsAll);

Log_Total_NO3_N =log(1.55); % mg L-1

Log_Total_NH4_N = log(1.55); % mg L-1

Log_Total_P = log(0.448); % mg L-1

Log_Total_Cu = log(0.0310); % mg L-1

Log_Total_Zn = log(0.230); % mg L-1

Log_Total_Pb = log(0.0150); % mg L-1

% NO3
b = bar(Log_MeanNO3_RunsAll);
nRuns = size(Log_MeanNO3_RunsAll,1); % Number of runs
grays = repmat(linspace(0.2,1.0,nRuns+1)',1,3);  % Create grayscale shades from light to dark
% Apply different gray to each run
for k = 1:nRuns+1 
    b(k).FaceColor = grays(k,:); 
end
xticklabels({'Clay','Glass','Sand','Shale','Shell'});
ylabel('log NO_{3}-N (mg L^{-1})');
yline(Log_Total_NO3_N, 'LineWidth',2,'Color', 'k');
ylim([-7 3]);
%legend({'Run 1','Run 2','Run 3','Run 4','Run 5','Run 6'}, 'Location', 'southeast');
set(gca,'FontSize',12,'FontWeight','bold');
ax = gca;
ax.Box = 'off';
% Clay
txt = 'c';
text(0.75,2.6,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'b';
text(0.88,1.15,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(1.02,-1.15,txt,'FontSize',12,'FontWeight','bold'); %R4
txt = 'a';
text(1.15,-4.1,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(1.29,-2.25,txt,'FontSize',12,'FontWeight','bold');%R6
% Glass
txt = 'b';
text(1.75,1.55,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'ab';
text(1.87,0.75,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'ab';
text(1.92,-1.33,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(2.15,-6.0,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'ab';
text(2.28,-1.42,txt,'FontSize',12,'FontWeight','bold');%R6
% Sand
txt = 'b';
text(2.75,0.75,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'a';
text(2.9,-1.4,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(3.02,-4.7,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(3.16,-3.75,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(3.3,-1.95,txt,'FontSize',12,'FontWeight','bold');%R6
% Shale
txt = 'c';
text(3.75,1.5,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'b';
text(3.9,0.95,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(4.01,-1.2,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(4.16,-2.95,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(4.3,-2.1,txt,'FontSize',12,'FontWeight','bold');%R6
% Shell
txt = 'c';
text(4.75,1.8,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'b';
text(4.875,1.1,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(5.01,-0.75,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(5.16,-3.15,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(5.3,-2.0,txt,'FontSize',12,'FontWeight','bold');%R6


% NH4
b = bar(Log_MeanNH4_RunsAll);
nRuns = size(Log_MeanNH4_RunsAll,1); % Number of runs
grays = repmat(linspace(0.2,1.0,nRuns+1)',1,3);  % Create grayscale shades from light to dark
% Apply different gray to each run
for k = 1:nRuns+1 
    b(k).FaceColor = grays(k,:); 
end
xticklabels({'Clay','Glass','Sand','Shale','Shell'});
ylabel('log NH_{4}-N (mg L^{-1})');
yline(Log_Total_NH4_N, 'LineWidth',2,'Color', 'k');
ylim ([-6.5 2]);
%legend({'Run 1','Run 2','Run 3','Run 4','Run 5','Run 6'}, 'Location', 'northeast');
set(gca,'FontSize',12,'FontWeight','bold');
ax = gca;
ax.Box = 'off';
%Clay
txt = 'a';
text(0.76,-5.785,txt,'FontSize',12,'FontWeight','bold'); %R2
txt = 'b';
text(0.89,-4.14,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'b';
text(1.02,-3.5,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'ab';
text(1.05,-4.45,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'ab';
text(1.25,-4.6,txt,'FontSize',12,'FontWeight','bold');%R6
% Glass
txt = 'a';
text(1.76,-4.85,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'a';
text(1.89,1.15,txt,'FontSize',12,'FontWeight','bold'); %R3
txt = 'a';
text(2.02,-4.15,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(2.15,-3.1,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(2.29,-4.55,txt,'FontSize',12,'FontWeight','bold');%R6
% Sand
txt = 'a';
text(2.76,-5.81,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'b';
text(2.89,-2.6,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'ab';
text(2.92,-3.725,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'ab';
text(3.06,-4.42,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'ab';
text(3.26,-4.65,txt,'FontSize',12,'FontWeight','bold');%R6
% Shale
txt = 'a';
text(3.75,-5.8,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'b';
text(3.88,-1.9,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'ab';
text(3.92,-4.2,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'ab';
text(4.065,-4.47,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'ab';
text(4.28,-4.65,txt,'FontSize',12,'FontWeight','bold');%R6
% Shell
txt = 'a';
text(4.76,-5.8,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'b';
text(4.90,-2.7,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'b';
text(5.02,-2.0,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'ab';
text(5.05,-4.42,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'ab';
text(5.275,-4.55,txt,'FontSize',12,'FontWeight','bold');%R6

% P
b = bar(Log_MeanP_RunsAll);
nRuns = size(Log_MeanP_RunsAll,1); % Number of runs
grays = repmat(linspace(0.2,1.0,nRuns+1)',1,3);  % Create grayscale shades from light to dark
% Apply different gray to each run
for k = 1:nRuns+1 
    b(k).FaceColor = grays(k,:); 
end
xticklabels({'Clay','Glass','Sand','Shale','Shell'});
ylabel('log P (mg L^{-1})');
yline(Log_Total_P, 'LineWidth',2,'Color', 'k');
set(gca,'FontSize',12,'FontWeight','bold');
ax = gca;
ax.Box = 'off';
%Clay
txt = 'a';
text(0.75,-5.95,txt,'FontSize',12,'FontWeight','bold'); %R2
txt = 'a';
text(0.89,-6.05,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(1.02,-5.16,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(1.15,-5.16,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(1.29,-5.16,txt,'FontSize',12,'FontWeight','bold');%R6
% Glass
txt = 'a';
text(1.75,-2.47,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'a';
text(1.89,-2.69,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(2.02,-2.42,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(2.145,-2.0,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(2.30,-1.59,txt,'FontSize',12,'FontWeight','bold');%R6
% Sand
txt = 'a';
text(2.75,1.36,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'b';
text(2.89,0.96,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'ab';
text(2.96,1.3,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'ab';
text(3.14,1.61,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'ab';
text(3.28,1.1,txt,'FontSize',12,'FontWeight','bold');%R6
% Shale
txt = 'ab';
text(3.65,-4.25,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'a';
text(3.87,-8.5,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'ab';
text(3.9875,-5.185,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'ab';
text(4.145,-4.1,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'b';
text(4.285,-2.9,txt,'FontSize',12,'FontWeight','bold');%R6
% Shell
txt = 'a';
text(4.75,-3.53,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'a';
text(4.89,-3.5,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(5.02,-3.38,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(5.16,-3.48,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'b';
text(5.3,-2.22,txt,'FontSize',12,'FontWeight','bold');%R6


% Cu
b = bar(Log_MeanCu_RunsAll);
nRuns = size(Log_MeanCu_RunsAll,1); % Number of runs
grays = repmat(linspace(0.2,1.0,nRuns+1)',1,3);  % Create grayscale shades from light to dark
% Apply different gray to each run
for k = 1:nRuns+1 
    b(k).FaceColor = grays(k,:); 
end
xticklabels({'Clay','Glass','Sand','Shale','Shell'});
ylabel('log Cu (mg L^{-1})');
yline(Log_Total_Cu, 'LineWidth',2,'Color', 'k');
set(gca,'FontSize',12,'FontWeight','bold');
ax = gca;
ax.Box = 'off';
%Clay
txt = 'a';
text(0.75,-4.85,txt,'FontSize',12,'FontWeight','bold'); %R2
txt = 'a';
text(0.89,-4.73,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(1.01,-4.95,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(1.15,-5.05,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(1.29,-5.22,txt,'FontSize',12,'FontWeight','bold');%R6
% Glass
txt = 'a';
text(1.75,-4.5,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'a';
text(1.87,-4.41,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(2.01,-4.36,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(2.15,-4.36,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(2.29,-4.55,txt,'FontSize',12,'FontWeight','bold');%R6
% Sand
txt = 'ab';
text(2.67,-5.0,txt,'FontSize',11,'FontWeight','bold');%R2
txt = 'a';
text(2.88,-5.155,txt,'FontSize',11,'FontWeight','bold');%R3
txt = 'ab';
text(3.01,-5.0,txt,'FontSize',11,'FontWeight','bold');%R4
txt = 'ab';
text(3.15,-4.7,txt,'FontSize',11,'FontWeight','bold');%R5
txt = 'b';
text(3.285,-4.5,txt,'FontSize',11,'FontWeight','bold');%R6
% Shale
txt = 'a';
text(3.75,-5.05,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'a';
text(3.89,-4.95,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(4.02,-5.25,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(4.15,-5.22,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(4.29,-5.15,txt,'FontSize',12,'FontWeight','bold');%R6
% Shell
txt = 'a';
text(4.75,-4.4,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'a';
text(4.89,-4.15,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(5.02,-4.4,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(5.155,-4.48,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(5.29,-4.65,txt,'FontSize',12,'FontWeight','bold');%R6



%Zn
b = bar(Log_MeanZn_RunsAll);
nRuns = size(Log_MeanZn_RunsAll,1); % Number of runs
grays = repmat(linspace(0.2,1.0,nRuns+1)',1,3);  % Create grayscale shades from light to dark
% Apply different gray to each run
for k = 1:nRuns+1 
    b(k).FaceColor = grays(k,:); 
end
xticklabels({'Clay','Glass','Sand','Shale','Shell'});
ylabel('log Zn (mg L^{-1})');
yline(Log_Total_Zn, 'LineWidth',2,'Color', 'k');
set(gca,'FontSize',12,'FontWeight','bold');
ax = gca;
ax.Box = 'off';
%Clay
txt = 'b';
text(0.75,0.2,txt,'FontSize',12,'FontWeight','bold'); %R2
txt = 'ab';
text(0.77,-1.3,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'ab';
text(0.92,-1.725,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(1.15,-1.925,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'ab';
text(1.28,-1.875,txt,'FontSize',12,'FontWeight','bold');%R6
% Glass
txt = 'b';
text(1.75,-0.5,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'ab';
text(1.75,-2.0,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'ab';
text(1.925,-2.2,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(2.15,-2.4,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'b';
text(2.3,-2.475,txt,'FontSize',12,'FontWeight','bold');%R6
% Sand
txt = 'ab';
text(2.65,-2.57,txt,'FontSize',11,'FontWeight','bold');%R2
txt = 'a';
text(2.89,-3.75,txt,'FontSize',11,'FontWeight','bold');%R3
txt = 'ab';
text(3.0,-3.75,txt,'FontSize',11,'FontWeight','bold');%R4
txt = 'ab';
text(3.15,-2.8,txt,'FontSize',11,'FontWeight','bold');%R5
txt = 'b';
text(3.285,-2.25,txt,'FontSize',11,'FontWeight','bold');%R6
% Shale
txt = 'b';
text(3.75,-0.7,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'a';
text(3.89,-1.6,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(4.02,-1.75,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(4.15,-2.25,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(4.285,-2.375,txt,'FontSize',12,'FontWeight','bold');%R6
% Shell
txt = 'b';
text(4.76,0.25,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'a';
text(4.89,-1.725,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(5.02,-2.25,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(5.15,-2.5,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(5.29,-2.6,txt,'FontSize',12,'FontWeight','bold');%R6

%Pb
b = bar(Log_MeanPb_RunsAll);
nRuns = size(Log_MeanPb_RunsAll,1); % Number of runs
grays = repmat(linspace(0.2,1.0,nRuns+1)',1,3);  % Create grayscale shades from light to dark
% Apply different gray to each run
for k = 1:nRuns+1 
    b(k).FaceColor = grays(k,:); 
end
xticklabels({'Clay','Glass','Sand','Shale','Shell'});
ylabel('log Pb (mg L^{-1})');
yline(Log_Total_Pb, 'LineWidth',2,'Color', 'k');
set(gca,'FontSize',12,'FontWeight','bold');
ax = gca;
ax.Box = 'off';
%Clay
txt = 'a';
text(0.75,-7.65,txt,'FontSize',12,'FontWeight','bold'); %R2
txt = 'a';
text(0.88,-7.5,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(1.02,-8.125,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(1.15,-8.75,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(1.29,-8.625,txt,'FontSize',12,'FontWeight','bold');%R6
% Glass
txt = 'a';
text(1.75,-7.15,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'a';
text(1.88,-7.625,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(2.02,-7.2,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(2.15,-7.9,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(2.29,-7.65,txt,'FontSize',12,'FontWeight','bold');%R6
% Sand
txt = 'a';
text(2.75,-10.2,txt,'FontSize',11,'FontWeight','bold');%R2
txt = 'a';
text(2.89,-8.375,txt,'FontSize',11,'FontWeight','bold');%R3
txt = 'a';
text(3.02,-7.7,txt,'FontSize',11,'FontWeight','bold');%R4
txt = 'ab';
text(3.15,-6.25,txt,'FontSize',11,'FontWeight','bold');%R5
txt = 'b';
text(3.29,-5.6,txt,'FontSize',11,'FontWeight','bold');%R6
% Shale
txt = 'b';
text(3.74,-7.95,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'ab';
text(3.79,-8.49,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(4.02,-8.69,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(4.15,-8.98,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(4.28,-8.9,txt,'FontSize',12,'FontWeight','bold');%R6
% Shell
txt = 'a';
text(4.75,-7.25,txt,'FontSize',12,'FontWeight','bold');%R2
txt = 'a';
text(4.89,-7.17,txt,'FontSize',12,'FontWeight','bold');%R3
txt = 'a';
text(5.02,-7.17,txt,'FontSize',12,'FontWeight','bold');%R4
txt = 'a';
text(5.15,-7.62,txt,'FontSize',12,'FontWeight','bold');%R5
txt = 'a';
text(5.28,-7.975,txt,'FontSize',12,'FontWeight','bold');%R6