function [DownLogRecord,Label ,Answer, PerformanceInt, PerformanceFloat] = Main( RecordNum )

%%%%%%%%%%%%%%% File READ %%%%%%%%%%%%
Extension = '_RR.csv';
FileName = [RecordNum Extension];
out = textread(FileName, '%s', 'whitespace',',');
HansRR = out(4:3:end);
Ans = out(6:3:end);
Label = [];
RR = [];
for i = 1:length(HansRR),
    val = cell2mat(HansRR(i));
    val = str2num(val);
    RR = [RR;val];
    LabelVal = cell2mat(Ans(i));
    if length(LabelVal) == 4',
        if LabelVal == 'AFIB',
            Label = [Label,1];
        else
            Label = [Label,0];
            123123123
        end
    else
        Label = [Label,0];
    end
end

%%%%%%%%%%%%%%% Downsampling %%%%%%%%%%%%

Freq = 20; 
DownLogRecord = DownLog(RR,Freq);
Label = Label(2:1:end);
Label = Label(1:Freq:end);

% figure(1);
% StrTemp = 'Raw volatility of ';
% title([StrTemp FileName]);
% plot(DownLogRecord);


%%%%%%%%%%%%%%% Markov Regime Switching GARCH %%%%%%%%%%%%
ORDERS = [1,1,1];
reg = 2;
flag = 2; 

[thetahat results struct] = swgarchest(DownLogRecord,flag,ORDERS,reg);



%%%%%%%%%%%%%%% Obtain Smoothing Param %%%%%%%%%%%%
Answer = ResultExtract(results.smooth);

%%%%%%%%%%%%%%% PLOT %%%%%%%%%%%%
PlotDraw(DownLogRecord,Answer, Label, RecordNum);



%%%%%%%%%%%%%%% Result %%%%%%%%%%%%
TP = 0;
TN = 0; 
FP = 0; 
FN = 0; 

for i = 1:length(Answer),
    MyAns = Answer(i);
    TrueAns = Label(i);
    if TrueAns == 0 && MyAns == 0,
        TN = TN+1;
    elseif TrueAns == 0 && MyAns == 1,
        FP = FP + 1; 
    elseif TrueAns == 1 && MyAns == 1, 
        TP = TP + 1;
    elseif TrueAns == 1 && MyAns == 0,
        FN = FN + 1;
    end
end
RecordNum
PerformanceInt = [TP, TN, FP, FN]
Acc = (TN + TP) / (TN + TP + FN + FP);
Se = TP / (TP + FN);
Spe = TN / (TN + FP);
Pp = TP / (TP + FP);

PerformanceFloat = [Acc, Se, Spe, Pp]
        


end

