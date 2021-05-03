% load pEMv2 results` 
%%Load results
period=0.5;
pix=0.117;%0.117; %pixel to um convertion 0.113 Casellas

 splitX = results.X;
 trackInfo = results.trackInfo;
 splitIndex = results.trackInfo.splitIndex;
 splitLength = trackInfo.splitLength;
 numFeatures = trackInfo.numFeatures;
 vacf_exp = trackInfo.vacf_exp;
 BIC = results.BIC;
 state = results.state;
 
 
 posteriorProb=results.posteriorProb;
 
     % estimated state sequence`

     [MAX,est_stateIndex] = max(posteriorProb,[],2);

    %MSD
        numLags = size(splitX{1},1)-1;

        [est_stateMSD,est_stateMSDerror,est_stateN] = AverageMSD(splitX,est_stateIndex,numLags);
        numStates=results.optimalSize;
        
         % plot MSD
    figure; hold on; box on;
    colorSet = hsv(numStates);
    timeLags = 1:numLags;
    timeLags=timeLags*period;
    legendname = 'h = legend(';
    for i = 1:numStates
        stderror = est_stateMSDerror(i,:)./sqrt(est_stateN(i));
        %errorbar(timeLags,est_stateMSD(i,:)*0.117*0.117,stderror*0.117*0.117,'color',colorSet(i,:),'linewidth',1.5);
        errorbar(timeLags,est_stateMSD(i,:)*pix*pix,stderror*pix*pix,'color',colorSet(i,:),'linewidth',1.5);

        legendname = [legendname '''State ' num2str(i) ''','];
    end
    set(gca,'fontsize',20,'linewidth',2);
    eval([legendname(1:end-1) ');']);
    set(h,'box','off','location','northwest','fontsize',20);
    xlabel('Time lags (steps)','fontsize',20);
    ylabel('MSD (\mum^2)','fontsize',20);

    
        est_popFraction = zeros(1,numStates);
    for i = 1:numStates
        est_popFraction(i) = sum(est_stateIndex == i)/length(est_stateIndex);
    end
    
    disp(['Est population fraction: ' num2str(est_popFraction)]);
    analytics.est_popFraction = est_popFraction;

