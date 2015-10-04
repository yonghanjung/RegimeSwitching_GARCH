function [ DownSample ] = DownLog( Data, SampleFreq )
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% DownSample 
LogRR = log(Data(2:1:end)) - log(Data(1:1:end-1));
DownSample = downsample(LogRR,SampleFreq);



end

