% dataleft=EEG.data;
% dataright=EEG.data;

config.chanloc=EEG.chanlocs;
%% PSD analysis
% BCI2000 power spectrue density analysis
modelOrder=18+round(fs/100);
binwidth=1;
numPerbin=10;
specstart=0;
specstop=40;
specbin=round((specstop-specstart)/binwidth)+1;

para=[modelOrder specstart specstop binwidth numPerbin 2 fs];
condleft=[];condright=[];

% condoneind=posmrk(label==1);
% condtwoind=posmrk(label==-1);
for ind=1:size(dataleft,3)
    [tspec1,freqbin]=mem(transpose(squeeze(double(dataleft(:,625:1000,ind)))),para);
    condleft(:,:,end+1)=mean(tspec1,3);
end
condleft=condleft(:,:,2:end);
for ind=1:size(dataright,3)        
    [tspec1,freqbin]=mem(transpose(squeeze(double(dataright(:,625:1000,ind)))),para);
    condright(:,:,end+1)=mean(tspec1,3); 
end
condright=condright(:,:,2:end);
ressq = calc_rsqu(double(condleft), double(condright), 1);
r1=mean(condleft,3);
r2=mean(condright,3);

figure
data2plot=ressq';
data2plot=cat(2, data2plot, zeros(size(data2plot, 1), 1));
data2plot=cat(1, data2plot, zeros(1, size(data2plot, 2)));
xData=freqbin-binwidth/2;
xData(end+1) = xData(end) + diff(xData(end-1:end));
surf(xData,1:size(data2plot,1),data2plot);
view(2);
axis tight
colormap jet;
colorbar;
%-------------------------------------------------------
freqbin2plot=1:size(ressq,1);% choose frebins to plot topographic
num2plot=length(freqbin2plot);
if num2plot==1
    figure
    topoplot(ressq(freqbin2plot,:),config.chanloc);
    title(['r^2 topographic of ',num2str(freqbin(freqbin2plot)),'Hz']);
    drawnow;
    colormap jet;
    colorbar;
end
if num2plot>1
    SIZEBOX=150;
    rownum=6;
    colnum=6;
    maxnum=rownum*colnum;
    nbfig=ceil(num2plot/(rownum*colnum));
    
    for figind=1:nbfig
        curplotnum=(num2plot>figind*maxnum)*maxnum+(1-(num2plot>figind*maxnum))*mod(num2plot,maxnum);
        rowcols=[round(sqrt(curplotnum)),ceil(sqrt(curplotnum))];
        if figind>1
            rowcols=[rownum colnum];
        end
        curfig=figure;
        pos=get(curfig,'Position');
        posx = max(0, pos(1)+(pos(3)-SIZEBOX*rowcols(2))/2);
        posy = pos(2)+pos(4)-SIZEBOX*rowcols(1);
        set(curfig,'Position', [posx posy  SIZEBOX*rowcols(2)  SIZEBOX*rowcols(1)]);
        for i=1:curplotnum
            curaxes=subplot(rowcols(1),rowcols(2),i,'Parent',curfig);
            tmpobj=topoplot(ressq(i+maxnum*(figind-1),:),config.chanloc,'maplimits','maxmin');
            title(['r^2 topographic of ',num2str(freqbin(i+maxnum*(figind-1))),'Hz']);
            drawnow;
        end
    end
    
    colormap jet;
    
    ColorbarHandle = cbar(0,0,get(gca, 'clim'));
    pos = get(ColorbarHandle,'position');  % move left & shrink to match head size
    set(ColorbarHandle,'position',[pos(1)-.05 pos(2)+0.13 pos(3)*0.8 pos(4)-0.26]);
end