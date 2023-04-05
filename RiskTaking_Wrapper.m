%% Run the main analyses and generate the figures from 'Risk taking for potential losses but not gains increases with time of day'
% Published in Scientific Reports (2023)
% Bedder, R.L*, Vaghi, M, Dolan, R.J & Rutledge R.B 

% * Corresponding author rlbedder@princeton.edu

% This code and the prospect theory models can be downloaded from
% https://github.com/RachelBedder/RiskTaking_ScientificReports

% The nain dataset can be downloaded from the repository
%https://datadryad.org/stash/dataset/doi:10.5061/dryad.prr4xgxkk

%% General Notes 

%03/26/2023     A previous version of matlab was used for the original
%paper analyses, and as such some of the numbers are rounded slightly
%different for the time of day (e.g., 4th decimal places). 

%% ...download the data from dryad

load('/Users/rb1725/Downloads/Rutledge_GBE_risk_data_with_TOD.mat')
clear depData

%% ...extract participants from UK (0-7) and USA (400)

locationID = extractfield(subjData,'location');

subjData = subjData((locationID>=0 & locationID<=7) | locationID==400);
locationID = locationID((locationID>=0 & locationID<=7) | locationID==400);

locationID(locationID<400) = 1;
locationID(locationID==400) = 2;

%% ... initialise some matrices for subject by play number

maxPlay             = 164; %...in the UK and USA sample, this is the maximum times a participant plays
allData.timeOfDay   = nan(length(subjData),maxPlay);
allData.version     = allData.timeOfDay;
allData.yearDay     = allData.timeOfDay;
allData.gambles{1}  = allData.timeOfDay;
allData.gambles{2}  = allData.timeOfDay;
allData.gambles{3}  = allData.timeOfDay;
allData.gambles{4}  = allData.timeOfDay;

%% ...extract relevant data from each subject and their plays

for subj = 1:length(subjData)
    
    subjPlays   =   subjData(subj);
    noPlays     =   length(subjPlays.data);
    
    for np = 1:noPlays;
        
        data   =    subjPlays.data{np};
                
        allData.timeOfDay(subj,np)  =   subjPlays.timeOfDay(np);
        allData.version(subj,np)    =   subjPlays.designVersion(np);
        allData.yearDay(subj,np)    =   subjPlays.dayNumber(np);
        
        allData.gambles{1}(subj,np) =   sum(data(data(:,3)>0,7))/sum(data(:,3)>0);
        allData.gambles{2}(subj,np) =   sum(data(data(:,3)<0,7))/sum(data(:,3)<0);
        allData.gambles{3}(subj,np) =   sum(data(data(:,3)==0,7))/sum(data(:,3)==0);
        allData.gambles{4}(subj,np) =   sum(data(:,7))/size(data,1);
    
    end
end
   
allData.location    =   locationID';
allData.gender      =   extractfield(subjData,'isFemale')';
allData.age         =   extractfield(subjData,'age')';

% convert the versions
allData.version(allData.version<=2)     =   1;
allData.version(allData.version>2)      =   2;

%% ...transform the time of day so 6am is the lowest value (0.25)

transTOD    =   allData.timeOfDay;

hour        =   1/24;
minTime     =   hour*6;

transTOD(locationID==2,:)   =   transTOD(locationID==2,:)-(hour*6.5);
transTOD(transTOD<=minTime) =   transTOD(transTOD<=minTime)+1;

allData.timeOfDay   =   transTOD;

%% ...load the prospect theory fits

load('allFits.mat')
allData.betas       =   allFits.betas;

%% ...run the script for each of the figures in the paper

%nSample  =  [10000]; %...for the permutated p-values (n 10,000 in the original paper)
nSample = [2]

[allStats.fig2] = RiskTaking_Figure2(allData,nSample);
[allStats.fig3] = RiskTaking_Figure3(allData,nSample);
[allStats.fig4] = RiskTaking_Figure4(allData,nSample);
