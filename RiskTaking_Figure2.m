function [statistics] = RiskTaking_Figure2(allData,nSample)

rand('seed',sum(100*clock));

%...extract the demographic vectors you need from the main struct
location        =       allData.location;
timeOfDay       =       allData.timeOfDay;
gender          =       allData.gender;
age             =       allData.age;
version         =       allData.version;
yearDay         =       allData.yearDay;

%%...make an index of the cells in runData
i.gain                =       1;
i.loss                =       2;
i.mix                 =       3;
i.riskGain            =       4;
i.riskLoss            =       5;
i.lambda              =       6;
i.mu                  =       7;

%...extract the data to run statistical analyses on with time of day
runData             =       [allData.gambles(1:3) allData.betas];
runData{i.lambda}   =       log(runData{i.lambda});

%...set the order you want to run the data in (swap these around if you
%want to run comparisons of effect sizes)
runOrder = [i.loss i.gain i.mix i.riskLoss i.riskGain i.lambda i.mu];

%...index the groups
i.PARTS     =           length(allData.timeOfDay);
i.UK        =   1;      kLabel{i.UK}    = 'UK';            kIdx{i.UK}     =  find(true(i.PARTS,1) & location==1);
i.US        =   2;      kLabel{i.US}    = 'USA';           kIdx{i.US}     =  find(location==2);
i.Male      =   3;      kLabel{i.Male} 	= 'Male';          kIdx{i.Male}   =  find(gender==0);
i.Female    =   4;      kLabel{i.Female}= 'Female';        kIdx{i.Female} =  find(gender==1);
i.Young     =   5;      kLabel{i.Young} = 'Young';         kIdx{i.Young}  =  find(age<=3) ;
i.Old       =   6;      kLabel{i.Old}   = 'Old';           kIdx{i.Old}    =  find(age>=4);
i.All       =   7;      kLabel{i.All}   = 'All';           kIdx{i.All}    =  find(location<3);
i.Ratio     =   8;      kLabel{i.Ratio}    = 'Ratio';      kIdx{i.Ratio}     =  find(version==1);
i.Uncorr    =   9;      kLabel{i.Uncorr}   = 'Uncorr';     kIdx{i.Uncorr}    =  find(version==2);

%...set the colours to be plotted in
colours{i.loss}       =       [[1 .4 .4];[1 0 0];[.6 0 0]];
colours{i.gain}       =       [[.4 1 .4];[0 .4 0];[0 .2 0]];
colours{i.mix}        =       [[.4 .7 1];[0 0 1];[0 .3 .6]];
colours{i.riskLoss}   =       [[153 76 0];[255 128 0];[255 153 51]]/255;
colours{i.riskGain}   =       [[0 102 102];[0 204 204];[0 255 255]]/255;
colours{i.lambda}     =       [[76 0 153];[127 0 255];[178 102 255]]/255;

%...make the 3.5 hour daily time difference bins for the figure (no
%analysis is based on this)
timeBins          =     [0 .1449 .2908 .4367 .5828];

%% For each participant find the first two eligible plays

%...add the time parameters
earlyTime         =     0.25+(0.0417*2);
lateTime          =     0.75+(0.0417*4);

withinTime  =  nan(i.PARTS,2);
withinGam   =  nan(i.PARTS,2);
withinPlay  =  nan(i.PARTS,2);

for count = runOrder; %...run in order for datasets given above (i.e loss trials, gain trials lambda fits)
    
    for part = 1:i.PARTS %...run for each participant
        
        pTime           =       timeOfDay(part,:);
        pData           =       runData{count}(part,:);
        pDay            =       yearDay(part,:);
        pVerse          =       version(part,:);

        
        %% ...find the first eligible play
        earlyPlay =  find(pTime >earlyTime & pTime <lateTime,1,'first');
        
        if ~isempty(earlyPlay) %...only run for those with a first eligible play
            
            timeWithin(1)  =       pTime(earlyPlay);
            gamWithin(1)   =       pData(earlyPlay);
            
            %% ...find the next eligible play
            
            %...transform all the data relative to first play
            diffTime      =       pTime-timeWithin;
            diffData      =       pData-gamWithin;
            diffDay       =       pDay-pDay(earlyPlay);
            
            %...find all other eligible plays
            timeIdx       =       find(pTime>earlyTime & pTime<lateTime); %...must be later than earliest time (8am) and earlier than the latest time (10pm)
            dayIdx        =       find(diffDay>=1); %...must be at least one day after the first play
            verseIdx      =       find(pVerse==pVerse(earlyPlay)); %...must be the same version as the first play
            
            allEligible   =       intersect(timeIdx,dayIdx); %...find the plays that are both a eligible time and a different day
            allEligible   =       intersect(allEligible,verseIdx); %...of those above plays find those that haven't changed version
            
            if ~isempty(allEligible)
                
                if timeWithin<pTime(min(allEligible)); %...if the first play is an earlier time, return in order [first eligible, second eligible]
                    
                    withinTime(part,:) =     [timeWithin pTime(min(allEligible))];
                    withinGam(part,:)  =     [gamWithin pData(min(allEligible))];
                    withinPlay(part,:) =     [earlyPlay min(allEligible)];
                    
                else %...if the second play is an earlier time, return in order [second eligible, first eligible]
                    
                    withinTime(part,:) =   [pTime(min(allEligible)) timeWithin];
                    withinGam(part,:)  =   [pData(min(allEligible)) gamWithin];
                    withinPlay(part,:) =  [min(allEligible) earlyPlay];
                    
                end
                
            end
            
        else %...don't run for those without eligible plays
            
        end
        
    end
    
    newGam(:,count)  = withinGam(:,2)-withinGam(:,1); %...for each data set, get the gambling difference
    saveGam{count} = withinGam; %...this preserves everything in the two column format to plot later
    
end

newTime(:,1)         = withinTime(:,2)-withinTime(:,1); %...get the time difference

%...re-index the demographic subgroups based on eligibility for within
%participant
for k  = 1:length(kIdx);
    
    kIdx{k}     = intersect(kIdx{k},find(~isnan(newTime)));
    
end



%% Run the analysis for each participant and each demongraphic group n times, including bootstrapping and permutation tests

for n = 1:nSample; %...iterate for n samples for permutation tests

for k = 1:length(kIdx); %...run for each subgroup of participants
    
    idx            =   kIdx{k};
    
    for count = runOrder;
        
            %...randomly permutate the gamble data for each participant (preserving
            %each participants data together)
            reData                  =       newGam(idx(randperm(length(idx))),:);
            rSample{k}(n,count)     =       corr(newTime(idx),reData(:,i.loss),'type','pearson');
                        
            clear reData
            
            %...sample without replacement for bootstrapped effect sizes
            bootParts                  =      datasample(idx,length(idx));
            reBoot{k}(n,count)         =      corr(newTime(bootParts),newGam(bootParts,count),'type','pearson'); clear bootParts
            
    end
end
              if mod(n,1000)==0
               disp(['***GBE_Figure2 ',int2str(n),'for index' int2str(k),'***'])
    end
     
end 

%%

for k = 1:length(kIdx); %...run for each subgroup of participants
    
    idx            =   kIdx{k};
    
    for count = runOrder;
        
        [Rstat(k,count) pValue(k,count)]    =       corr(newTime(idx),newGam(idx,count),'type','pearson');%...get real effect size and p value
        pPerm(k,count)                    =       sum(Rstat(k,count)>rSample{k}(:,count))/nSample;%...generate a permutated p value
        pPerm2(k,count)                 =       (sum(abs(Rstat(k,count))<abs(rSample{k}(:,count)))/nSample)  ;    %...get the permutated p value
 
        sortBoot                          =       sort(reBoot{k}(:,count));
        lowerBound(k,count)               =       Rstat(k,count) - sortBoot(ceil(nSample*0.025));%...get the bootstrapped lower confidence bound
        upperBound(k,count)               =       sortBoot(ceil(nSample*0.975)) - Rstat(k,count);%...get the bootstrapped upper confidence bound
        
        clear sortBoot
        
        %...generate means and standard errors for each time of day bin to plot
        for tb = 1:length(timeBins)-1;
            
            todIdx                         =      intersect(idx,find(newTime>=timeBins(tb) & newTime<=timeBins(tb+1)));
            meanData{k}(tb,count)          =      mean(newGam(todIdx,count));
            seData{k}(tb,count)            =      std(newGam(todIdx,count))/sqrt(length(idx));
            
        end
        
        
    end; clear count
end; clear k idx

%...make tables of all relevant statistics for each gender group
stats               =   {'Rstat','P value','Permu P value','Permu P 2','BS Lower Bound','BS Upper Bound', 'N Participants'}';
for k = 1:length(kIdx)
    
    riskGain        =       [Rstat(k,i.riskGain) pValue(k,i.riskGain) pPerm(k,i.riskGain) pPerm2(k,i.riskGain) lowerBound(k,i.riskGain) upperBound(k,i.riskGain) length(kIdx{k})]';
    riskLoss        =       [Rstat(k,i.riskLoss) pValue(k,i.riskLoss) pPerm(k,i.riskLoss) pPerm2(k,i.riskLoss) lowerBound(k,i.riskLoss) upperBound(k,i.riskLoss) length(kIdx{k})]';
    lambda          =       [Rstat(k,i.lambda) pValue(k,i.lambda) pPerm(k,i.lambda)  pPerm2(k,i.lambda) lowerBound(k,i.lambda) upperBound(k,i.lambda) length(kIdx{k})]';
    mu              =       [Rstat(k,i.mu) pValue(k,i.mu) pPerm(k,i.mu) pPerm2(k,i.mu) lowerBound(k,i.mu) upperBound(k,i.mu) length(kIdx{k})]';
    gainPlays       =       [Rstat(k,i.gain) pValue(k,i.gain) pPerm(k,i.gain) pPerm2(k,i.gain) lowerBound(k,i.gain) upperBound(k,i.gain) length(kIdx{k})]';
    lossPlays       =       [Rstat(k,i.loss) pValue(k,i.loss) pPerm(k,i.loss) pPerm2(k,i.loss) lowerBound(k,i.loss) upperBound(k,i.loss) length(kIdx{k})]';
    mixPlays        =       [Rstat(k,i.mix) pValue(k,i.mix) pPerm(k,i.mix)  pPerm2(k,i.mix) lowerBound(k,i.mix) upperBound(k,i.mix) length(kIdx{k})]';

    
    statistics.(strcat('WithinSubjects_',kLabel{k})) = table(stats,riskGain,riskLoss,lambda,mu,lossPlays,gainPlays,mixPlays);clear riskGain riskLoss lambda mu
    
end; clear k rSample reBoot


%% Generate the figures

%% ...generate figure 4a
figure; hold on
bar([1],Rstat(i.All,i.gain),'FaceColor',colours{i.gain}(2,:),'EdgeAlpha',0)
bar([2],Rstat(i.All,i.loss),'FaceColor',colours{i.loss}(2,:),'EdgeAlpha',0)
bar([3],Rstat(i.All,i.mix),'FaceColor',colours{i.mix}(2,:),'EdgeAlpha',0)
errorbar([1:3],Rstat(i.All,[i.gain i.loss i.mix]),lowerBound(i.All,[i.gain i.loss i.mix]),upperBound(i.All,[i.gain i.loss i.mix]),'k','LineStyle', 'none')
ylabel('Within-subject effect size');ylim([-0.06 0.12])
set(gca,'xtick',[1 2 3],'xticklabel',{'Gain','Loss','Mix'},'ytick',-.06:0.02:0.12);xlim([0 4])
axis square
legend({'Gain','Loss','Mix'},'location','eastoutside')
title([kLabel{i.All},' N = ',num2str(length(kIdx{i.All}))])
pbaspect([1 2 1])

%% ...generate figure 4b
figure; hold on
bar([1],Rstat(i.All,i.riskGain),'FaceColor',colours{i.riskGain}(2,:),'EdgeAlpha',0)
bar([2],Rstat(i.All,i.riskLoss),'FaceColor',colours{i.riskLoss}(2,:),'EdgeAlpha',0)
bar([3],Rstat(i.All,i.lambda),'FaceColor',colours{i.lambda}(2,:),'EdgeAlpha',0)
errorbar([1:3],Rstat(i.All,[i.riskGain i.riskLoss i.lambda]),lowerBound(i.All,[i.riskGain i.riskLoss i.lambda]),upperBound(i.All,[i.riskGain i.riskLoss i.lambda]),'k','LineStyle', 'none')
ylabel('Within-subject effect size');ylim([-0.12 0.06])
set(gca,'xtick',[1 2 3],'xticklabel',{'\alpha+','\alpha-','log(\lambda)'},'ytick',-.12:0.02:0.06);xlim([0 4])
axis square
legend({'\alpha+','\alpha-','log(\lambda)'},'location','eastoutside')
title([kLabel{i.All},' N = ',num2str(length(kIdx{i.All}))])
pbaspect([1 2 1])

%% Illustrative plot
earlyTime        =     0.25+(0.0417*2);    %...from 8am
midTime          =     0.75+(0.0417*-3);   %...3pm (halfway between 8am and 10pm)
lateTime         =     0.75+(0.0417*4);    %...10pm

newTime          =   withinTime(intersect(find(~isnan(withinTime)),kIdx{i.All})); 
tmpTime = withinTime(intersect(find(~isnan(withinTime)),kIdx{i.All}),:);
tmpTime = tmpTime(:,2)-tmpTime(:,1);


for count = [1:6] %...plot for the betas
    
    earlyIdx    =   newTime>earlyTime & newTime<midTime;
    lateIdx     =   newTime>midTime & newTime<lateTime;
    
    timeDiff = [5];
    tmpEarly = tmpTime<(0.0417*timeDiff);
    tmpLate = tmpTime>=(0.0417*timeDiff);
    
    newBeta       =   saveGam{count}(intersect(find(~isnan(saveGam{count})),kIdx{i.All}),:);
    tmpBeta = newBeta(:,2)-newBeta(:,1); %...gets the within difference
    
    meanBeta(count,[1 2])     = [mean(tmpBeta(tmpEarly)) mean(tmpBeta(tmpLate))];
    errorBeta(count,[1 2])    = [std(tmpBeta(tmpEarly)) std(tmpBeta(tmpLate))]./sqrt(sum([tmpEarly tmpLate]));
    
end



figure;
subplot(1,3,1)
plot([0 3],[0 0],'LineWidth',1,'Color','k'); hold on
errorbar([1 2],meanBeta(i.riskGain,:)-meanBeta(i.riskGain,1),errorBeta(i.riskGain,:),'Color',colours{i.riskGain}(2,:),'LineWidth',3); hold on
xlim([0 3]);set(gca,'xtick',[1 2],'xticklabel',{'0-5','>5'});
xlabel('Time Difference (h)')
ylim([-0.06 0.06]);
ylabel('Gain sensitivity \alpha+ (relative to 0-5 (h))','FontSize',15)
set(gca,'ytick',[-0.06:0.03:0.06])
title([kLabel{i.All},' N = ',num2str(length(kIdx{i.All}))])
xtickangle(45)
% 
subplot(1,3,2)
plot([0 3],[0 0],'LineWidth',1,'Color','k'); hold on
errorbar([1 2],meanBeta(i.riskLoss,:)-meanBeta(i.riskLoss,1),errorBeta(i.riskLoss,:),'Color',colours{i.riskLoss}(2,:),'LineWidth',3); hold on
xlim([0 3]);set(gca,'xtick',[1 2],'xticklabel',{'0-5','>5'});
xlabel('Time Difference (h)')
ylim([-0.06 0.06]);
ylabel('Loss sensitivity \alpha- (relative to 0-5 (h))','FontSize',15);
set(gca,'ytick',[-0.06:0.03:0.06])
title([kLabel{i.All},' N = ' ,num2str(length(kIdx{i.All}))])
xtickangle(45)
% 
subplot(1,3,3)
plot([0 3],[0 0],'LineWidth',1,'Color','k'); hold on
errorbar([1 2],meanBeta(i.lambda,:)-meanBeta(i.lambda,1),errorBeta(i.lambda,:),'Color',colours{i.lambda}(2,:),'LineWidth',3); hold on
xlim([0 3]);set(gca,'xtick',[1 2],'xticklabel',{'0-5','>5'});
xlabel('Time Difference (h)')
ylim([-0.12 0.12]);
ylabel('Loss aversion log(\lambda) (relative to 0-5 (h))','FontSize',15)
set(gca,'ytick',[-0.12:0.06:0.12])
title([kLabel{i.All},' N = ',num2str(length(kIdx{i.All}))])
xtickangle(45)


end
