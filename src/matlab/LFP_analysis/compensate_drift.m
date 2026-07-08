function [Compensated]= compensate_drift(data,SF);

FreqRange = [0.1];
[b,a]=cheby2 (2,20,FreqRange./SF/2,'low');
Filtered=filtfilt(b,a,data);
Compensated=data-Filtered;
