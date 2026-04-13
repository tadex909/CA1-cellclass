%select_your data
fpath = 'N:\Vinca\MATLAB\Data\All_data';
main_folder='N:\Vinca\MATLAB\Data\'
h=readtable(strcat(main_folder,'select_sess.xlsx'))
keep_dat=h.session_name(logical(h.doselect))

%color infos
loadcolmathSFN
load colormaprom.mat

%load data
load(strcat(fpath, '\supool.mat'));

%groups
a17=find(supool.agegrp==17);
a18=find(supool.agegrp==18);
g1=sort([a17;a18]);
a19=find(supool.agegrp==19);
a20=find(supool.agegrp==20);
g2=sort([a19;a20]);
a21=find(supool.agegrp==21);
a22=find(supool.agegrp==22);
g3=sort([a21;a22]);
a23=find(supool.agegrp==23);
a24=find(supool.agegrp==24);
g4=sort([a23;a24]);
g5=find(supool.agegrp==25);

velocity=nan(5,5);
for i = 17:25
%VELOCITY
vel = supool.vel_cond(eval(strcat('a',string(i))));
velocity(i-16,1)=vel
%REWARD FREQUENCY
finderr=[(eval(strcat('a',string(i))))*2;eval(strcat('a',string(i)))*2+1]
rw1=supool.freq_rewardsess(finderr)
end