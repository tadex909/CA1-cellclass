function [data, TS]=Read_OEP_Binary_altered (filename,chan, start, stop, subsample)
%[Data, TS]=Read_OEP_Binary (filename, chans, start, stop)
% chan is a vector specifying the channels you want to read starting from 1
% it's the channel in the order it appears in the binary file
% not the real a/d channel
% start and stop in seconds
% subsample= the number of elements to skip between samples
% Pierre-Pascal lenck-Santini INMED Marseille


% go to the continuous.dat folder 
% !! si on va dans cont.dat NOW on peut pas ouvrir OEBIN la ligne d'après !! 

opath=pwd;
% cd ..
% cd ..

fname = 'structure.oebin';
fid = fopen(fname); 
raw = fread(fid,inf); 
str = char(raw'); 
fclose(fid); 
val = jsondecode(str);
header=val.continuous(1);
%cd (opath)
% [val]=get_oebinFile ();
SF=header.sample_rate;
nChans=header.num_channels;
D.Header = header;
 
if start*SF==0 %where to start recording
    beg=1;
else
   beg= round(start*SF);
end

    


if max(chan) > nChans
    disp(['there are only ', num2str(nChans,'%d'), 'channels'])
    
end


D.Timestamps = readNPY(strcat(opath,  filesep, 'continuous\Rhythm_FPGA-', num2str(val.continuous.source_processor_id), '.', num2str(val.continuous.source_processor_sub_idx),filesep,'timestamps.npy'));% no idea what units these are
format long g
TS=double(D.Timestamps)./SF;
contFile=fullfile(filename);
file=dir(contFile);
samples=file.bytes/2/header.num_channels;
D.Data=memmapfile(contFile,'Format',{'int16' [header.num_channels samples] 'mapped'});
if stop==-1
    data=double(D.Data.Data.mapped(chan,beg:subsample:end)).*header.channels(1).bit_volts; 
    TS=TS(beg:subsample:end);
else
    data=double(D.Data.Data.mapped(chan,beg:subsample:round(stop*SF))).*header.channels(1).bit_volts; 
    TS=TS(beg:subsample:round(stop*SF));
end