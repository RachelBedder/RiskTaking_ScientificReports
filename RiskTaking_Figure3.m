function [statistics] = RiskTaking_Figure3(allData,nSample)

rand('seed',sum(100*clock));

%...extract the vectors you need from the main struct
location        =       allData.location;
timeOfDay       =       allData.timeOfDay(:,1);
gender          =       allData.gender;
age             =       allData.age;
version         =       allData.version(:,1);

%%...make an index of the cells in runData
i.gain                =       1;     
i.loss                =       2;     
i.mix                 =       3;     

%...extract the data to run statistical analyses on with time of day
runData         =       [allData.gambles{1}(:,1) allData.gambles{2}(:,1) allData.gambles{3}(:,1)];

%...set the order you want to run the data in (swap these around if you
%want to run comparisons of effect sizes
runOrder = [i.loss i.gain i.mix];

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

%...make the 4 hour daily time bins for the figure (no analysis is based on
%this)
timeBins          =     [.25 .50 .75 1 1.25];


%% Run the analysis for each participant and each demongraphic group n times, including bootstrapping and permutation tests

for n = 1:nSample; %...iterate for n samples for permutation tests
    
    for k = 1:length(kIdx); %...run for each subgroup of participants
        
        idx            =   kIdx{k};

        for count = runOrder; %...run in order for datasets given above (i.e loss trials, gain trials lambda fits)
                           
            %...randomly permutate the gamble data for each participant (keeping
            %their loss, gain and mixed percentage together)
            reData                     =       runData(idx(randperm(length(idx))),:);       
            rSample{k}(n,count)        =       corr(timeOfDay(idx),reData(:,count),'type','pearson'); clear reData
                        
            %...sample without replacement for bootstrapped effect sizes
            bootParts                  =      datasample(idx,length(idx));
            reBoot{k}(n,count)         =      corr(timeOfDay(bootParts),runData(bootParts,count),'type','pearson');
        end
        
    end; clear k idx
    
    if mod(n,1000)==0
        disp(['***GBE_Figure3 ',int2str(n),'for index','***'])
    end
    
    
end
%%

for k = 1:length(kIdx);
    
    for count = runOrder;
        
        idx            =   kIdx{k};

        [Rstat(k,count) pValue(k,count)]  =       corr(timeOfDay(idx),runData(idx,count),'type','pearson'); clear xData yData           %...get real effect size and p value
        pPerm(k,count)                  =       (sum(Rstat(k,count)<rSample{k}(:,count))/nSample)  ;    %...get the permutated p value
        pPerm2(k,count)                 =       (sum(abs(Rstat(k,count))<abs(rSample{k}(:,count)))/nSample)  ;    %...get the permutated p value
        
        sortBoot                        =       sort(reBoot{k}(:,count));
        lowerBound(k,count)             =       Rstat(k,count) - sortBoot(ceil(nSample*0.025));     %...get the bootstrapped lower confidence bound
        upperBound(k,count)             =       sortBoot(ceil(nSample*0.975)) - Rstat(k,count);     %...get the bootstrapped upper confidence bound
        
        
        clear sortBoot
        
        %...generate means and standard errors for each time of day bin to plot
        for tod = 1:length(timeBins)-1
            
            todIdx                          =      intersect(idx,find(timeOfDay>=timeBins(tod) & timeOfDay<=timeBins(tod+1)));
            meanGam{k}(tod,count)          =      mean(runData(todIdx,count));
            seGam{k}(tod,count)            =      std(runData(todIdx,count))/sqrt(length(todIdx));
            
         end
        

    end
end; clear k idx


%...make tables of all relevant statistics for each gender group
stats               =   {'Rstat','P value','Permu P value','Permu P 2 side','BS Lower Bound','BS Upper Bound','N Participants'}';
for k = 1:length(kIdx);
    
    gainPlays       =       [Rstat(k,i.gain) pValue(k,i.gain) pPerm(k,i.gain) pPerm2(k,i.gain) lowerBound(k,i.gain) upperBound(k,i.gain) length(kIdx{k})]';
    lossPlays       =       [Rstat(k,i.loss) pValue(k,i.loss) pPerm(k,i.loss) pPerm2(k,i.loss) lowerBound(k,i.loss) upperBound(k,i.loss) length(kIdx{k})]';
    mixPlays        =       [Rstat(k,i.mix) pValue(k,i.mix) pPerm(k,i.mix) pPerm2(k,i.mix) lowerBound(k,i.mix) upperBound(k,i.mix) length(kIdx{k})]';
    
    statistics.(strcat('GambingTimeofDay_',kLabel{k})) = table(stats,gainPlays,lossPlays,mixPlays);clear LossPlays GainPlays MixPlays
    
end; clear k

% 
% j = 0;
% for k = [7 5 4];
%     
%     j = j+1;
%     
%     groupidx{j,1}     =   kLabel{k};
%     N(j,1)            =   length(kIdx{k})
%     gainPlays_R(j,1)  =   Rstat(k,i.gain);      gainPlays_pVal(j,1)  =   pPerm2(k,i.gain);
%     lossPlays_R(j,1)  =   Rstat(k,i.loss);      lossPlays_pVal(j,1)  =   pPerm2(k,i.loss);
%     mixPlays_R(j,1)  =   Rstat(k,i.mix);        mixPlays_pVal(j,1)  =   pPerm2(k,i.mix);
%     gainPlays(j,1)  =   Rstat(k,i.gain_evMax); gain_evMax_pVal(j,1)  =   pPerm2(k,i.gain_evMax);
%     loss_evMax(j,1) =   Rstat(k,i.loss_evMax); loss_evMax_pVal(j,1)  =   pPerm2(k,i.loss_evMax);
%     mix_evMax(j,1) =   Rstat(k,i.mix_evMax); mix_evMax_pVal(j,1)  =   pPerm2(k,i.mix_evMax);
%     gain_points(j,1)  =   Rstat(k,i.gain_points); gain_points_pVal(j,1)  =   pPerm2(k,i.gain_points);
%     loss_points(j,1) =   Rstat(k,i.loss_points); loss_points_pVal(j,1)  =   pPerm2(k,i.loss_points);
%     mix_points(j,1) =   Rstat(k,i.mix_points); mix_points_pVal(j,1)  =   pPerm2(k,i.mix_points);
% 
% end
% % table(groupidx,N,gainPlays_R,gainPlays_pVal,lossPlays_R,lossPlays_pVal,mixPlays_R,mixPlays_pVal,gain_evMax,gain_evMax_pVal,loss_evMax,loss_evMax_pVal,mix_evMax,mix_evMax_pVal,points,points_pVal)
% table(groupidx,N,gainPlays,gain_evMax_pVal,loss_evMax,loss_evMax_pVal,mix_evMax,mix_evMax_pVal,...
%     gain_points,gain_points_pVal,loss_points,loss_points_pVal,mix_points,mix_points_pVal)



%% ...generate figure 3a

figure;
plot([0 5],[0 0],'k'); hold on
e1 = errorbar([1:4]-0.2,meanGam{i.All}(:,i.gain)-meanGam{i.All}(1,i.gain),seGam{i.All}(:,i.gain),'Color',colours{i.gain}(2,:),'LineWidth',4);
e3 = errorbar([1:4]+0.2,meanGam{i.All}(:,i.mix)-meanGam{i.All}(1,i.mix),seGam{i.All}(:,i.mix),'Color',colours{i.mix}(2,:),'LineWidth',4);
e2 = errorbar([1:4],meanGam{i.All}(:,i.loss)-meanGam{i.All}(1,i.loss),seGam{i.All}(:,i.loss),'Color',colours{i.loss}(2,:),'LineWidth',4);
set(gca,'xtick',[1:4],'xticklabel',{'6am','Midday','6pm','Midnight'})
set(gca,'ytick',[-0.02:0.02:0.06],'yticklabel',{'-2%','0%','2%','4%','6%'})
ylim([-0.02 0.06]);ylabel('Percent Gambles chosen')
xlim([0 5]);xlabel('Time of day');
axis square
legend([e1 e2 e3],{'Gain','Loss','Mix'},'location','eastoutside')
title([kLabel{i.All},' N = ',num2str(length(kIdx{i.All}))]);
xtickangle(45)
set(gca,'TickDir','out');


%% ...generate figure 3b
dataA = i.Female;
dataB = i.Male;


figure; hold on
bar([1],Rstat(dataA,i.gain),'FaceColor',colours{i.gain}(2,:),'EdgeAlpha',0)
bar([2],Rstat(dataA,i.loss),'FaceColor',colours{i.loss}(2,:),'EdgeAlpha',0)
bar([3],Rstat(dataA,i.mix),'FaceColor',colours{i.mix}(2,:),'EdgeAlpha',0)
bar([5],Rstat(dataB,i.gain),'FaceColor',colours{i.gain}(2,:),'EdgeAlpha',0)
bar([6],Rstat(dataB,i.loss),'FaceColor',colours{i.loss}(2,:),'EdgeAlpha',0)
bar([7],Rstat(dataB,i.mix),'FaceColor',colours{i.mix}(2,:),'EdgeAlpha',0)
errorbar([1:3],Rstat(dataA,:),lowerBound(dataA,:),upperBound(dataA,:),'k','LineStyle', 'none')
errorbar([5:7],Rstat(dataB,:),lowerBound(dataB,:),upperBound(dataB,:),'k','LineStyle', 'none')
ylabel('Between-subject effect sizes');ylim([-0.03 0.06])
set(gca,'xtick',[2 6],'xticklabel',{[kLabel{dataA},' N = ',num2str(length(kIdx{dataA}))],[kLabel{dataB},' N = ',num2str(length(kIdx{dataB}))]...
    },'ytick',-0.03:0.03:0.06);xtickangle(45);xlim([0 8]);
axis square
set(gca,'TickDir','out');



end
