function [ output_args ] = PlotDraw( Data, Answer )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

figure()
hold on

for i = 1:length(Data),
    if Answer(i) == 0,
        plot(i,Data(i),'bo')
    else,
        plot(i,Data(i),'ro')
    end
end




