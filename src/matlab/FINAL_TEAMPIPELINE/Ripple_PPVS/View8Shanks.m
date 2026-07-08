
function varargout = View8Shanks(varargin)
% VIEW8SHANKS MATLAB code for View8Shanks.fig
%      VIEW8SHANKS, by itself, creates a new VIEW8SHANKS or raises the existing
%      singleton*.
%
%      H = VIEW8SHANKS returns the handle to a new VIEW8SHANKS or the handle to
%      the existing singleton*.
%
%      VIEW8SHANKS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VIEW8SHANKS.M with the given input arguments.
%
%      VIEW8SHANKS('Property','Value',...) creates a new VIEW8SHANKS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before View8Shanks_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to View8Shanks_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help View8Shanks

% Last Modified by GUIDE v2.5 02-Feb-2021 18:30:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @View8Shanks_OpeningFcn, ...
    'gui_OutputFcn',  @View8Shanks_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before View8Shanks is made visible.
function View8Shanks_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to View8Shanks (see VARARGIN)

% Choose default command line output for View8Shanks
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global GuiData
chanread=zeros(64,1);
GuiData=struct('Cont_Path','D:\Rachel\RAM10\RAM10_2022-12-13_11-16-53\Record Node 103\experiment1\recording1\continuous\Rhythm_FPGA-100.0',...
    'TStart', 0,...
    'TWidth',1,...
    'Order', [25 29 30 1 2 3 27 6 31 8 4 24 5 22 7 20 23 18 21 16 19 11 17 9 15 13 10 12 14 28 26 32 33 39 37 51 53 55 52 50 56 48 54 46 49 44 47 42 45 58 43 60 41 61 57 34 59 38 62 63 64 35 36 40],...
    'OrderNew', [25 29 30 1 2 3 27 6 31 8 4 24 5 22 7 20 23 18 21 16 19 11 17 9 15 13 10 12 14 28 26 32 33 39 37 51 53 55 52 50 56 48 54 46 49 44 47 42 45 58 43 60 41 61 57 34 59 38 62 63 64 35 36 40],...
    'Order2',[1:64],...
    'SpacingShks',500,...
    'SF',30000,...
    'D',[],...
    'Edge',[],...
    'Pz',[],...
    'SkRd',chanread,...
    'T',[],...
    'EvtTS',[],...
    'RecDur',[],...
    'HighPass',0.1,...
    'LowPass',9000,...
    'Filter',1,...
    'EvtNum', [],...
    'FirstTime',[],...
    'RealChan',strings(64,1),...
    'MemFile',[],...
    'FileType',[],...
    'prefix','100_CH');



% UIWAIT makes View8Shanks wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = View8Shanks_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in radiobutton1.
function radiobutton1_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton1
global GuiData
val=get(hObject,'Value');
if(val==1)
    GuiData.SkRd(1:8)=1;
else
    GuiData.SkRd(1:8)=0;
end

read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth)
PlotAll (handles)

% --- Executes on button press in radiobutton2.
function radiobutton2_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton2
global GuiData
val=get(hObject,'Value');
if(val==1)
    GuiData.SkRd(9:16)=1;
else
    GuiData.SkRd(9:16)=0;
end

read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth)
PlotAll (handles)

% --- Executes on button press in radiobutton3.
function radiobutton3_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton3
global GuiData
val=get(hObject,'Value');
if(val==1)
    GuiData.SkRd(17:24)=1;
else
    GuiData.SkRd(17:24)=0;
end
read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth)
PlotAll (handles)


% --- Executes on button press in radiobutton4.
function radiobutton4_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton4
global GuiData
val=get(hObject,'Value');
if(val==1)
    GuiData.SkRd(25:32)=1;
else
    GuiData.SkRd(25:32)=0;
end
read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth)
PlotAll (handles)

% --- Executes on button press in radiobutton5.
function radiobutton5_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton5
global GuiData
val=get(hObject,'Value');
if(val==1)
    GuiData.SkRd(33:40)=1;
else
    GuiData.SkRd(33:40)=0;
end
read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth)
PlotAll (handles)

% --- Executes on button press in radiobutton6.
function radiobutton6_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton6
global GuiData
val=get(hObject,'Value');
if(val==1)
    GuiData.SkRd(41:48)=1;
else
    GuiData.SkRd(41:48)=0;
end
read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth)
PlotAll (handles)
% --- Executes on button press in radiobutton7.
function radiobutton7_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton7
global GuiData
val=get(hObject,'Value');
if(val==1)
    GuiData.SkRd(49:56)=1;
else
    GuiData.SkRd(49:56)=0;
end
read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth)
PlotAll (handles)
% --- Executes on button press in radiobutton8.
function radiobutton8_Callback(hObject, eventdata, handles)
% hObject    handle to radiobutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobutton8
global GuiData
val=get(hObject,'Value');
if(val==1)
    GuiData.SkRd(57:64)=1;
else
    GuiData.SkRd(57:64)=0;
end
read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth)
PlotAll (handles)


function PlotAll (handles)

global GuiData
hold (handles.axes1,'off')
hold (handles.axes2,'off')
hold (handles.axes3,'off')
hold (handles.axes4,'off')
hold (handles.axes5,'off')
hold (handles.axes6,'off')
hold (handles.axes7,'off')
hold (handles.axes8,'off')





for ch=1:64
    
    if ch<=8
        if GuiData.SkRd(ch)==1
            plot(handles.axes1,GuiData.T,GuiData.D(:,ch)-(GuiData.SpacingShks*(ch-1)),'Color','k')
        else
            cla(handles.axes1)
        end
        hold (handles.axes1,'on')
        grid (handles.axes1,'on')
        
    elseif ch<=16
        if  GuiData.SkRd(ch)==1
            plot(handles.axes2,GuiData.T,GuiData.D(:,ch)-(GuiData.SpacingShks*(ch-1-8)),'Color','k')
        else
            cla(handles.axes2)
        end
        hold (handles.axes2,'on')
        grid (handles.axes2,'on')
    elseif ch<=24
        if  GuiData.SkRd(ch)==1
            plot(handles.axes3,GuiData.T,GuiData.D(:,ch)-(GuiData.SpacingShks*(ch-1-16)),'Color','k')
        else
            cla(handles.axes3)
        end
        hold (handles.axes3,'on')
        grid (handles.axes3,'on')
        grid on
    elseif ch<=32
        if  GuiData.SkRd(ch)==1
            plot(handles.axes4,GuiData.T,GuiData.D(:,ch)-(GuiData.SpacingShks*(ch-1-24)),'Color','k')
        else
            cla(handles.axes4)
        end
        hold (handles.axes4,'on')
        grid (handles.axes4,'on')
        grid on
        set(handles.axes4,'Xlim',[GuiData.T(1) GuiData.T(end)])
    elseif ch<=40
        if GuiData.SkRd(ch)==1
            plot(handles.axes5,GuiData.T,GuiData.D(:,ch)-(GuiData.SpacingShks*(ch-1-32)),'Color','k')
        else
            cla(handles.axes5)
        end
        hold (handles.axes5,'on')
        grid (handles.axes5,'on')
        grid on
        set(handles.axes5,'Xlim',[GuiData.T(1) GuiData.T(end)])
        
    elseif ch<=48
        if GuiData.SkRd(ch)==1
            plot(handles.axes6,GuiData.T,GuiData.D(:,ch)-(GuiData.SpacingShks*(ch-1-40)),'Color','k')
        else
            cla(handles.axes6)
        end
        hold (handles.axes6,'on')
        grid (handles.axes6,'on')
        grid on
        set(handles.axes6,'Xlim',[GuiData.T(1) GuiData.T(end)])
    elseif ch<=56
        if GuiData.SkRd(ch)==1
            plot(handles.axes7,GuiData.T,GuiData.D(:,ch)-(GuiData.SpacingShks*(ch-1-48)),'Color','k')
        else
            cla(handles.axes7)
        end
        hold (handles.axes7,'on')
        grid (handles.axes7,'on')
        grid on
        set(handles.axes7,'Xlim',[GuiData.T(1) GuiData.T(end)])
    elseif ch<=64
        if GuiData.SkRd(ch)==1
            plot(handles.axes8,GuiData.T,GuiData.D(:,ch)-(GuiData.SpacingShks*(ch-1-56)),'Color','k')
        else
            cla(handles.axes8)
        end
        hold (handles.axes8,'on')
        grid (handles.axes8,'on')
        grid on
        set(handles.axes8,'Xlim',[GuiData.T(1) GuiData.T(end)])
    end
end
yscal=[1:8];
yscal=(yscal*GuiData.SpacingShks*-1)+GuiData.SpacingShks;

set(handles.axes1,'YTick',yscal(end:-1:1))
set(handles.axes1,'yTickLabel', GuiData.RealChan(8:-1:1))
% yyaxis right
yticks(yscal(end:-1:1))
% yyaxis left

set(handles.axes2,'YTick',yscal(end:-1:1))
set(handles.axes2,'yTickLabel',  (16:-1:7))
set(handles.axes3,'YTick',yscal(end:-1:1))
set(handles.axes3,'yTickLabel',  (24:-1:17))
set(handles.axes4,'YTick',yscal(end:-1:1))
set(handles.axes4,'yTickLabel',(32:-1:25))
set(handles.axes5,'YTick',yscal(end:-1:1))
set(handles.axes5,'yTickLabel', (40:-1:33))
set(handles.axes6,'YTick',yscal(end:-1:1))
set(handles.axes6,'yTickLabel',(48:-1:41))
set(handles.axes7,'YTick',yscal(end:-1:1))
set(handles.axes7,'yTickLabel',(56:-1:49))
set(handles.axes8,'YTick',yscal(end:-1:1))
set(handles.axes8,'yTickLabel',(64:-1:57))

% if you want to see the real channel names:
% set(handles.axes1,'YTick',yscal(end:-1:1))
% set(handles.axes1,'yTickLabel', GuiData.RealChan(8:-1:1))
% % yyaxis right
% yticks(yscal(end:-1:1))
% % yyaxis left
% 
% set(handles.axes2,'YTick',yscal(end:-1:1))
% set(handles.axes2,'yTickLabel',  GuiData.RealChan(16:-1:7))
% set(handles.axes3,'YTick',yscal(end:-1:1))
% set(handles.axes3,'yTickLabel',  GuiData.RealChan(24:-1:17))
% set(handles.axes4,'YTick',yscal(end:-1:1))
% set(handles.axes4,'yTickLabel', GuiData.RealChan(32:-1:25))
% set(handles.axes5,'YTick',yscal(end:-1:1))
% set(handles.axes5,'yTickLabel', GuiData.RealChan(40:-1:33))
% set(handles.axes6,'YTick',yscal(end:-1:1))
% set(handles.axes6,'yTickLabel', GuiData.RealChan(48:-1:41))
% set(handles.axes7,'YTick',yscal(end:-1:1))
% set(handles.axes7,'yTickLabel', GuiData.RealChan(56:-1:49))
% set(handles.axes8,'YTick',yscal(end:-1:1))
% set(handles.axes8,'yTickLabel', GuiData.RealChan(64:-1:57))

% linkaxes([handles.axes1,handles.axes2,handles.axes3,handles.axes4],'xy');
linkaxes([handles.axes1,handles.axes2,handles.axes3,handles.axes4,handles.axes5,handles.axes6,handles.axes7,handles.axes8],'x');
% set(handles.Title,'string',GuiData.Cont_Path)
set(handles.TStart,'string',GuiData.TStart)
% if ~isempty(GuiData.TTLTS)
% set(handles.EventString,'string',[num2str(length(GuiData.TTLTS),'%d'),' Events'])
% end
drawnow

% --- Executes on button press in LoadFiles.
function LoadFiles_Callback(hObject, eventdata, handles)
% hObject    handle to LoadFiles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global GuiData
% set(handles.Title,'string','WAIT')
drawnow
[name, GuiData.Cont_Path]=uigetfile({'*.dat';'*.continuous'},'Select one File','Multiselect','off');

if contains(name,'.dat')
    GuiData.FileType= '.dat';
    GuiData.prefix='';
    GuiData.Order=[1:64];
    
elseif contains(name,'100_CH')
    GuiData.FileType= '100_CH';
    GuiData.prefix='100_CH';
    GuiData.Order=GuiData.OrderNew;
elseif  contains(name,'100_')
    GuiData.FileType= '100_';
    GuiData.prefix='100_';
    GuiData.Order=GuiData.OrderNew;
end

switch(GuiData.FileType)
    
    case '.dat'
        disp('it is a datfile')
        cd (GuiData.Cont_Path)
        opath=pwd;
        
        fname = 'structure.oebin';
        if isfile(fname)
        fid = fopen(fname);
        raw = fread(fid,inf);
        str = char(raw');
        fclose(fid);
        val = jsondecode(str);
        header=val.continuous(1);
        cd (opath)
        % [val]=get_oebinFile ();
        SF=header.sample_rate;
        nChans=header.num_channels;
        D.Header = header;
        contFile=fullfile(pwd,name);
        s=dir(contFile);
        samples=s.bytes/2/header.num_channels;
        D.Data=memmapfile(contFile,'Format',{'int16' [header.num_channels samples] 'mapped'});
%         s=dir('continuous.dat');
        SizeF = s.bytes/(nChans*2);
        TSdata = readNPY('timestamps.npy');% units are in sample
        GuiData.FirstTime=TSdata(1)/SF;%in seconds
        GuiData.RecDur=SizeF/SF;
        GuiData.SF=SF;
        GuiData.MemFile=D;
        GuiData.BigTS=double(TSdata)./SF;
        else
            errordlg('File Structure.oebin is missing in current directory')
        end
    case '100_CH'
        disp('it is a continuous file with _100_CH')
        NAM=[GuiData.prefix,num2str(GuiData.Order(1),'%d'),'.continuous'];
        [~, timestamps, info] = load_open_ephys_data_faster(NAM);
        GuiData.RealChan (1)=info.header.channel;
        GuiData.RecDur=timestamps(end)-timestamps(1); % in seconds
        GuiData.FirstTime=timestamps(1);
        
    case '100_'
        disp('it is a continuous file with 100_')
        NAM=[GuiData.prefix,num2str(GuiData.Order(1),'%d'),'.continuous'];
        [~, timestamps, info] = load_open_ephys_data_faster(NAM);
        GuiData.RealChan (1)=info.header.channel;
        GuiData.RecDur=timestamps(end)-timestamps(1); % in seconds
        GuiData.FirstTime=timestamps(1);
end


cd (GuiData.Cont_Path)

handles.slider1.Max=GuiData.RecDur;
handles.slider1.Min=0;
set(handles.TextSlider,'String',['Total: ',num2str(handles.slider1.Max,'%4.1f')])
handles.slider1.SliderStep=[(GuiData.TWidth/2)/GuiData.RecDur (GuiData.TWidth)/GuiData.RecDur];

set(handles.Title,'string',GuiData.Cont_Path)
set(handles.Title,'FontSize',8)
read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth,handles)
PlotAll (handles)


function TStart_Callback(hObject, eventdata, handles)
% hObject    handle to TStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of TStart as text
%        str2double(get(hObject,'String')) returns contents of TStart as a double
global GuiData

GuiData.TStart=str2double(get(hObject,'String'));
handles.slider1.Value=GuiData.TStart;
read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth,handles)
PlotAll (handles)

% --- Executes during object creation, after setting all properties.
function TStart_CreateFcn(hObject, eventdata, handles)
% hObject    handle to TStart (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function width_Callback(hObject, eventdata, handles)
% hObject    handle to width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of width as text
%        str2double(get(hObject,'String')) returns contents of width as a double
global GuiData
GuiData.TWidth=str2double(get(hObject,'String'));
handles.slider1.SliderStep=[(GuiData.TWidth/2)/GuiData.RecDur (GuiData.TWidth)/GuiData.RecDur];
read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth,handles)
PlotAll (handles)

% --- Executes during object creation, after setting all properties.
function width_CreateFcn(hObject, eventdata, handles)
% hObject    handle to width (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function spacing_Callback(hObject, eventdata, handles)
% hObject    handle to spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of spacing as text
%        str2double(get(hObject,'String')) returns contents of spacing as a double
global GuiData

GuiData.SpacingShks=str2double(get(hObject,'String'));
read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth)
PlotAll (handles)

% --- Executes during object creation, after setting all properties.
function spacing_CreateFcn(hObject, eventdata, handles)
% hObject    handle to spacing (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function read_all (Start,Stop,handles)
global GuiData
% GuiData.Order=[GuiData.Order,[33:64]];
cd (GuiData.Cont_Path)


switch(GuiData.FileType)
    
    
    case '.dat'
        
        if GuiData.TWidth > 10
            subsample=100;
        elseif GuiData.TWidth  > 1
            subsample=10;
        elseif GuiData.TWidth > 0.5
            subsample=1;
            
        else
            subsample=1;
        end
        
        if GuiData.TStart*GuiData.SF==0 %where to start recording
            beg=1;
        else
            beg= round(GuiData.TStart*GuiData.SF);
        end
        stop=GuiData.TWidth+GuiData.TStart;
        
        data=double(GuiData.MemFile.Data.Data.mapped(1:64,beg:subsample:round(stop*GuiData.SF))).*GuiData.MemFile.Header.channels(1).bit_volts;
        TS=GuiData.BigTS(beg:subsample:round(stop*GuiData.SF)); %chercher le TS dans D!!!
        
        TS=TS-TS(1);
        TS=TS+GuiData.TStart;
        GuiData.T=TS;
        GuiData.D=data';
        if GuiData.Filter %%% make sure the SF is 
            for ch=1:64
                
                format long g
                if GuiData.LowPass>(GuiData.SF/subsample)/2
                    GuiData.LowPass=(GuiData.SF/subsample)/2;
                    disp('!!!! Careful I changed lowpass because it was too high for the subsampling')
                end
                Wn=[GuiData.HighPass GuiData.LowPass]/(GuiData.SF/subsample)/2;
                [b,a]=butter(2,Wn,'bandpass');
                dtmp2=filtfilt(b,a,data(ch,:));
                GuiData.D(:,ch)=dtmp2;
            end
            
            
        end
       GuiData.RealChan=[GuiData.MemFile.Header.channels.source_processor_index];
       GuiData.RealChan=GuiData.RealChan+1;
       
       % if you want the mapped version ask for
        GuiData.RecChan= GuiData.MemFile.Header.channels.recorded_processor_index;
       
    otherwise
        NAM=[GuiData.prefix,num2str(GuiData.Order(1),'%d'),'.continuous'];
        subsample=1;
        if GuiData.SF>=30000
            if Stop-Start>=1000
                subsample=10; %3000 Hz
            end
            if Stop-Start>=10000
                subsample=30; %1000Hz
            end
            
            if Stop-Start>=30000
                subsample=100; %300Hz
            end
        end
        [t, data, ~, ~] = read_cont_PPSub(NAM,GuiData.TStart,GuiData.TStart+GuiData.TWidth,subsample);
        % [timestamps,samples,info,newfreq] = read_cont_PPSub (filename,startT,stopT,subsample)
        
        
        
        t=t-t(1);
        GuiData.T=t+GuiData.TStart;
        GuiData.D=zeros(length(t),32);
        GuiData.Edge=zeros(length(t),32);
        GuiData.Pz=zeros(length(t),2);
        GuiData.D(:,1)=data;
        
        for ch=1:64
            if (GuiData.SkRd(ch)==1)
                
                NAM=[GuiData.prefix,num2str(GuiData.Order(ch),'%d'),'.continuous'];
                [ttmp, dtmp,info,~] = read_cont_PPSub (NAM,Start,Stop,subsample);
                GuiData.RealChan (ch)=info.header.channel;
                %     [ttmp, dtmp, ~] = read_cont_PP(NAM,Start,Stop);
                ttmp=ttmp-ttmp(1);
                ttmp=ttmp+Start;
                if length(ttmp)>size(GuiData.T,1)
                    dtmp(size(GuiData.T,1)+1:end)=[];
                end
                if GuiData.Filter
                    format long g
                    Wn=[GuiData.HighPass GuiData.LowPass]/(GuiData.SF)/2;
                    [b,a]=butter(2,Wn,'bandpass');
                    dtmp2=filtfilt(b,a,dtmp);
                    
                end
                GuiData.D(:,ch)=dtmp2;
            else
                GuiData.D(:,ch)=zeros(size(data))*NaN;
            end
        end
end





function HighPass_Callback(hObject, eventdata, handles)
% hObject    handle to HighPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of HighPass as text
%        str2double(get(hObject,'String')) returns contents of HighPass as a double
global GuiData
GuiData.HighPass=str2double(get(hObject,'String'));
read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth)
PlotAll (handles)


% --- Executes during object creation, after setting all properties.
function HighPass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to HighPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function LowPass_Callback(hObject, eventdata, handles)
% hObject    handle to LowPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of LowPass as text
%        str2double(get(hObject,'String')) returns contents of LowPass as a double
global GuiData

GuiData.LowPass=str2double(get(hObject,'String'));
read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth)
PlotAll (handles)

% --- Executes during object creation, after setting all properties.
function LowPass_CreateFcn(hObject, eventdata, handles)
% hObject    handle to LowPass (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in NextEv.
function NextEv_Callback(hObject, eventdata, handles)
% hObject    handle to NextEv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global GuiData
if GuiData.EvtNum<length(GuiData.EvtTS)
    GuiData.EvtNum=GuiData.EvtNum+1;
    GuiData.TStart=GuiData.EvtTS(GuiData.EvtNum,1)-(GuiData.TWidth/2);
    
    read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth)
    PlotAll (handles)
end

% --- Executes on button press in PrevEv.
function PrevEv_Callback(hObject, eventdata, handles)
% hObject    handle to PrevEv (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global GuiData
if GuiData.EvtNum>1
    GuiData.EvtNum=GuiData.EvtNum-1;
    GuiData.TStart=GuiData.EvtTS(GuiData.EvtNum,1)-(GuiData.TWidth/2);
    
    read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth)
    PlotAll (handles)
end



function gotoEvt_Callback(hObject, eventdata, handles)
% hObject    handle to gotoEvt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gotoEvt as text
%        str2double(get(hObject,'String')) returns contents of gotoEvt as a double
global GuiData

GuiData.EvtNum=str2double(get(hObject,'String'));
GuiData.TStart=GuiData.EvtTS(GuiData.EvtNum,1)-(GuiData.TWidth/2);

read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth)
PlotAll (handles)

% --- Executes during object creation, after setting all properties.
function gotoEvt_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gotoEvt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in LoadEvts.
function LoadEvts_Callback(hObject, eventdata, handles)
% hObject    handle to LoadEvts (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global GuiData
% cd (abfpath)
[Sharpname, abfpath]=uigetfile('*.txt','Select Cluster Sharp Files','Multiselect','off');
cd (abfpath)
[data] = importdata(Sharpname);

GuiData.EvtTS=data(:,1);
GuiData.TStart=data(1,1)-(GuiData.TWidth/2);
GuiData.EvtNum=1;
read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth)
PlotAll (handles)
% Times=Times*1E-6;




% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
global GuiData
GuiData.TStart=get(hObject,'Value');
% disp(num2str(GuiData.TStart,'%5.5f'))

set(handles.TStart,'String',num2str(GuiData.TStart,'%4.1f'))
read_all (GuiData.TStart,GuiData.TStart+GuiData.TWidth)
PlotAll (handles)


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
global GuiData

set(hObject,'Value',0);
set(hObject,'SliderStep',[0.01 0.01]);
