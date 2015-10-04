function [ output_args ] = PlotDraw( Data, Answer, Label, RecordNum )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

MIN = min(Data);
MAX = max(Data);
Domain = 1:length(Data);

figure()
subplot(3,1,1);
plot(Domain, Data,'color',[0,0,0.8]);

TempStr = 'Raw Volatility of ';
title([TempStr RecordNum],'FontSize',20)
xlabel('time (s) downsampled', 'FontSize', 16)
ylabel('Log(RR_i) - Log(RR_i-1)')
ylim([-1.5,1.5])





subplot(3,1,2);
hold on

for i = 1:length(Data),
    if Label(i) == 1,
        plot([i,i],[MIN-1 MAX+1],'color',[1.0,0.8,0.8],'LineWidth',0.01)  
    end
end
plot(Domain, Data,'color',[0,0,0.8]);
TempStr = 'True regimes of ';
title([TempStr RecordNum],'FontSize',20)
xlabel('time (s) downsampled', 'FontSize', 16)
ylabel('Log(RR_i) - Log(RR_i-1)')
ylim([-1.5,1.5])
hold off 




subplot(3,1,3);
hold on

for i = 1:length(Data),
    if Answer(i) == 1,
        plot([i,i],[MIN-1 MAX+1],'color',[1.0,0.8,0.8],'LineWidth',0.01)  
    end
end
plot(Domain, Data,'color',[0,0,0.8]);
hold off 
TempStr = 'Estimated regimes of ';
title([TempStr RecordNum],'FontSize',20)
xlabel('time (s), downsampled', 'FontSize', 16)
ylabel('Log(RR_i) - Log(RR_i-1)')
ylim([-1.5,1.5])







