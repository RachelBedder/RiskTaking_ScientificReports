function [statistics] = RiskTaking_Figure4(allData,nSample)

rand('seed',sum(100*clock));

%...extract the demographic vectors you need from the main struct
location        =       allData.location;
timeOfDay       =       allData.timeOfDay(:,1);
gender          =       allData.gender;
age             =       allData.age;
version         =       allData.version(:,1);

%%...make an index of the cells in betaData
i.riskGain            =       1;
i.riskLoss            =       2;
i.lambda              =       3;
i.mu                  =       4;

%...extract the data to run statistical analyses on with time of day

betas           =       allData.betas;
betaData        =       [betas{i.riskGain}(:,1) betas{i.riskLoss}(:,1) log(betas{i.lambda}(:,1)) betas{i.mu}(:,1)];

%...set the order you want to run the data in (swap these around if you
%want to run comparisons of effect sizes)
runOrder = [i.riskLoss i.riskGain i.lambda i.mu];


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
colours{i.riskLoss}     =       [[153 76 0];[255 128 0];[255 153 51]]/255;
colours{i.riskGain}       =       [[0 102 102];[0 204 204];[0 255 255]]/255;
colours{i.lambda}     =       [[76 0 153];[127 0 255];[178 102 255]]/255;



%...make the 4 hour daily time bins
timeBins          =     [.25 .50 .75 1 1.25];

%% Run the analysis for each participant and each demongraphic group n times, including bootstrapping and permutation tests

for n = 1:nSample; %...iterate for n samples for permutation tests
    
    for k = 1:length(kIdx); %...run for each subgroup of participants
        
        idx            =   kIdx{k};
 
        for count = runOrder
            
            %...randomly permutate the gamble data for each participant (keeping
            %their loss, gain and mixed percentage together)
            reBeta                     =       betaData(idx(randperm(length(idx))),:);
            rSample{k}(n,count)        =       corr(timeOfDay(idx),reBeta(:,count),'type','pearson'); clear reBeta
            
            %...sample without replacement for bootstrapped effect sizes
            bootParts                  =      datasample(idx,length(idx));
            reBoot{k}(n,count)         =      corr(timeOfDay(bootParts),betaData(bootParts,count),'type','pearson');
            
        end
        
      
    end; clear k idx
    
      if mod(n,1000)==0
            disp(['***GBE_Figure4 ',int2str(n),'***'])
        end
        
end
    
    %%...generate the new statistics and add them to a table
    for k = 1:length(kIdx);
        
                for count = runOrder

        
        idx            =   kIdx{k};
        
        [Rstat(k,count) pValue(k,count)] =       corr(timeOfDay(idx),betaData(idx,count),'type','pearson');%...get real effect size and p value
        pPerm(k,count)                  =       sum(Rstat(k,count)>rSample{k}(:,count))/nSample;%...generate a permutated p value
        pPerm2(k,count)                 =       (sum(abs(Rstat(k,count))<abs(rSample{k}(:,count)))/nSample)  ;    %...get the permutated p value
        
        sortBoot                         =       sort(reBoot{k}(:,count));
        lowerBound(k,count)              =       Rstat(k,count) - sortBoot(ceil(nSample*0.025));%...get the bootstrapped lower confidence bound
        upperBound(k,count)              =       sortBoot(ceil(nSample*0.975)) - Rstat(k,count);%...get the bootstrapped upper confidence bound
        
        clear sortBoot
        
        %...generate means and standard errors for each time of day bin to plot
        for tod = 1:length(timeBins)-1
            
            todIdx                          =      intersect(idx,find(timeOfDay>=timeBins(tod) & timeOfDay<=timeBins(tod+1)));
            meanBeta{k}(tod,count)          =      mean(betaData(todIdx,count));
            seBeta{k}(tod,count)            =      std(betaData(todIdx,count))/sqrt(length(todIdx));
            
        end
        
        
    end; clear kidx
end; clear count

%% ...make tables of all relevant statistics for each gender group
stats               =   {'Rstat','P value','Permu P value','Permu P 2 side','BS Lower Bound','BS Upper Bound'}';
for k = 1:length(kIdx)
    
    riskGain        =       [Rstat(k,i.riskGain) pValue(k,i.riskGain) pPerm(k,i.riskGain) pPerm2(k,i.riskGain) lowerBound(k,i.riskGain) upperBound(k,i.riskGain)]';
    riskLoss        =       [Rstat(k,i.riskLoss) pValue(k,i.riskLoss) pPerm(k,i.riskLoss) pPerm2(k,i.riskLoss) lowerBound(k,i.riskLoss) upperBound(k,i.riskLoss)]';
    lambda          =       [Rstat(k,i.lambda) pValue(k,i.lambda) pPerm(k,i.lambda) pPerm2(k,i.lambda) lowerBound(k,i.lambda) upperBound(k,i.lambda)]';
    mu              =       [Rstat(k,i.mu) pValue(k,i.mu) pPerm(k,i.mu) pPerm2(k,i.mu) lowerBound(k,i.mu) upperBound(k,i.mu)]';
    
    statistics.(strcat('BetaTimeofDay_',kLabel{k})) = table(stats,riskGain,riskLoss,lambda,mu);clear riskGain riskLoss lambda mu
    
end; clear k


%% ...generate figure 4a
% 
dataA = i.Female;
dataB = i.Male;

figure;
plot([0 5],[0 0],'k'); hold on
e1 = errorbar([1:4]-0.1,meanBeta{dataA}(:,i.riskLoss)-meanBeta{dataA}(1,i.riskLoss),seBeta{dataA}(:,i.riskLoss),'Color',colours{i.riskLoss}(2,:),'LineWidth',3);
e2 = errorbar([1:4]+0.1,meanBeta{dataB}(:,i.riskLoss)-meanBeta{dataB}(1,i.riskLoss),seBeta{dataB}(:,i.riskLoss),'Color',colours{i.riskLoss}(2,:),'LineWidth',3,'LineStyle','--');
set(gca,'xtick',[1:4],'xticklabel',{'6am','Midday','6pm','Midnight'},'ytick',-0.08:0.02:0.02)
ylim([-0.08 0.02]);ylabel('Loss sensitivity \alpha-')
xlim([0 5]);xlabel('Time of day');
axis square
legend([e1 e2],{'Female','Male'})
xtickangle(45)
set(gca,'TickDir','out');

saveas(gcf,'GBE_Figure4a_ScientificReports.svg')


% 
%% ...generate figure 4b

figure
plot([0 5],[0 0],'k'); hold on
e1 = errorbar([1:4]-0.1,meanBeta{dataA}(:,i.lambda)-meanBeta{dataA}(1,i.lambda),seBeta{dataA}(:,i.lambda),'Color',colours{i.lambda}(2,:),'LineWidth',3);
e2 = errorbar([1:4]+0.1,meanBeta{dataB}(:,i.lambda)-meanBeta{dataA}(1,i.lambda),seBeta{dataB}(:,i.lambda),'Color',colours{i.lambda}(2,:),'LineWidth',3,'LineWidth',3,'LineStyle','--');
set(gca,'xtick',[1:4],'xticklabel',{'6am','Midday','6pm','Midnight'},'ytick',-0.04:0.04:0.12)
ylim([-0.04 0.12]);ylabel('Loss aversion log(\lambda)')
xlim([0 5]);xlabel('Time of day');
axis square
legend([e1 e2],{'Female','Male'})
xtickangle(45)
set(gca,'TickDir','out');

saveas(gcf,'GBE_Figure4b_ScientificReports.svg')


%% ...generate figure 4c

figure; hold on
bar([1],Rstat(dataA,i.riskGain),'FaceColor',colours{i.riskGain}(2,:),'EdgeAlpha',0)
bar([2],Rstat(dataA,i.riskLoss),'FaceColor',colours{i.riskLoss}(2,:),'EdgeAlpha',0)
bar([3],Rstat(dataA,i.lambda),'FaceColor',colours{i.lambda}(2,:),'EdgeAlpha',0)
bar([5],Rstat(dataB,i.riskGain),'FaceColor',colours{i.riskGain}(2,:),'EdgeAlpha',0)
bar([6],Rstat(dataB,i.riskLoss),'FaceColor',colours{i.riskLoss}(2,:),'EdgeAlpha',0)
bar([7],Rstat(dataB,i.lambda),'FaceColor',colours{i.lambda}(2,:),'EdgeAlpha',0)
errorbar([1:3],Rstat(dataA,[1:3]),lowerBound(dataA,[1:3]),upperBound(dataA,[1:3]),'k','LineStyle', 'none')
errorbar([5:7],Rstat(dataB,[1:3]),lowerBound(dataB,[1:3]),upperBound(dataB,[1:3]),'k','LineStyle', 'none')
ylabel('Effect sizes for time of day');ylim([-0.06 0.06])
set(gca,'xtick',[2 6],'xticklabel',{[kLabel{dataA},' N = ',num2str(length(kIdx{dataA}))],[kLabel{dataB},' N = ',num2str(length(kIdx{dataB}))]},...
    'ytick',-0.06:0.02:0.06);xtickangle(45);xlim([0 8])
axis square
xtickangle(45)


end
