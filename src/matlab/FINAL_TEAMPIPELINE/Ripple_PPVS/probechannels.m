function [channellist]=probechannels(probetype,channel);
%Select all the channels that are oon 1 shank of a seen spike
if probetype == 5
    if channel <= 12
        channellist=[1,2,3,4,5,6,7,8,9,10,11,12];
        
    elseif channel  <= 24 & channel >  12;
           channellist= [13,14,15,16,17,18,19,20,21,22,23,24];
     elseif channel <= 40 & channel >  24;
        channellist=[25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40];
    elseif channel <= 52 & channel >  40;
        channellist=[41,42,43,44,45,46,47,48,49,50,51,52];
    elseif  channel >  52;
        channellist=[53,54,55,56,57,58,59,60,61,62,63,64];
    end
elseif probetype == 64
    if channel <= 8
        channellist = [1:8];
    elseif channel  <= 16 & channel >  8;
        channellist = [9:16];
            elseif channel  <= 24 & channel > 16;
        channellist = [17:24];
            elseif channel  <= 32 & channel >  24;
        channellist = [25:32];
        elseif channel  <= 40 & channel >  32;
        channellist = [33:40];
            elseif channel  <= 48 & channel >  40;
        channellist = [41:48];
            elseif channel  <= 56 & channel >  48;
        channellist = [49:56];
            elseif channel  <= 64 & channel >  56;
        channellist = [57:64];
    end
elseif probetype == 77 
    if  ismember(channel,[16,21,18,23,20,24,27,30,19,29,17,25,26,32,28,22])
        channellist = [16,21,18,23,20,24,27,30,19,29,17,25,26,32,28,22];
    elseif ismember(channel,[11,1,13,3,15,5,31,7,14,9,12,10,8,6,4,2])
        channellist=[11,1,13,3,15,5,31,7,14,9,12,10,8,6,4,2];
    elseif ismember(channel,[54,64,52,62,50,60,34,58,51,56,53,55,57,59,61,63])
        channellist=[54,64,52,62,50,60,34,58,51,56,53,55,57,59,61,63];
    elseif ismember(channel,[49,44,47,42,45,41,38,35,46,36,48,40,39,33,37,43])
        channellist = [49,44,47,42,45,41,38,35,46,36,48,40,39,33,37,43];
    end
   
    end