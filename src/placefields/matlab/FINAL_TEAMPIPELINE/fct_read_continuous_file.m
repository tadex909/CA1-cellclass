% -----------------------------
% Written by MARTI Geoffrey
% 12/03/15
% 09/15
% 03/16
% -----------------------------

% This function reads ".continuous" files in the path "file_path" and
% returns the data ("data") over time ("time"), the "header"
% (recording informations) and the number of 1024-integers sub-recordings ("record_nb")

% Cette fonction permet de lire les données d'un fichier ".continuous" situé dans "file_path" et
% de retourner les données ("data") au cours du temps ("time"), le "header"
% (informations sur l'enregistrement) et le nombre de sous-enregistrements de 1024 entiers ("record_nb")


% Construction du fichier ".continuous" :
% 1) Header : Informations sur l'enregistrements
% 1024 caractères, chacun codé en 8 bits = 1 octet, pour une taille totale de 1024*1 = 1024 octets
% 2) timestamp : temps associé aux 1024 points d'un sous-enregistrement
% 1 entier, codé en 64 bits = 8 octets, pour une taille totale de 8 octets
% 3) sample_nb : nombre de samples
% 1 entier, codé en 16 bits = 2 octets, pour une taille totale de 2 octets
% 4) record_ind : numéro de l'enregistrement
% 1 entier, codé en 16 bits = 2 octets, pour une taille totale de 2 octets
% 5) samples : l'ensemble des points acquis lors d'un sous-enregistrement (1024 points)
% 1024 entiers, chacun codé en 16 bits = 2 octets, pour une taille totale
% de 2*1024 = 2048 octets
% 6) record_marker : ?
% 10 entiers, chacun codé en 8 bits = 1 octet, pour une taille totale de 10 octets



function [data, time, header, record_nb] = fct_read_continuous_file(file_path)

fid = fopen(file_path);
[~, file_name] = fileparts(file_path);

% Tailles en octet
file_size = fct_getfilesize(fid);
header_size = 1024*1;
timestamps_size = 1*8;
sample_nb_size = 1*2;
record_ind_size = 1*2;
samples_size = 1024*2;
record_marker_size = 10*1;
record_size = timestamps_size + sample_nb_size + record_ind_size + samples_size + record_marker_size;

string_nb = 1024;


hdr = fread(fid, string_nb, 'char*1');
eval(char(hdr')); 
clear hdr

% Useful to create .dat file, otherwise we lose information
% if strcmp(file_name(1:6), '100_CH')
% data = int16(1); % Be careful, LFP will be coded as 16-bits integers
% end

i = 1;
while ftell(fid) + record_size <= file_size
    timestamp(i) = fread(fid, 1, 'int64', 0, 'l');
    samples_nb(i) = fread(fid, 1, 'uint16',0,'l');
    if i > 1
        samples_nb_sum(i) =  samples_nb(i) + samples_nb_sum(i-1);
    else
        samples_nb_sum(i) = samples_nb(i);
    end

    record_ind(i) = fread(fid, 1, 'uint16');
    data((samples_nb_sum(i) - samples_nb(i) + 1):samples_nb_sum(i)) = fread(fid, samples_nb(i), 'int16', 0, 'b');
    time((samples_nb_sum(i) - samples_nb(i) + 1):samples_nb_sum(i)) = timestamp(i):timestamp(i)+samples_nb(i)-1;

    fread(fid, 10, 'char*1'); 
    i = i + 1;
end

fclose(fid);

record_nb = length(timestamp);

% Conversion à partir du Header
time = time / header.sampleRate; % en secondes
data = data * header.bitVolts; % en microVolts


