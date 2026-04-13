% -----------------------------
% Written by MARTI Geoffrey
% 10/03/15
% -----------------------------



function filesize = fct_getfilesize(fid)

fseek(fid,0,'eof');
filesize = ftell(fid);
fseek(fid,0,'bof');

end