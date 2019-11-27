function varargout = select_noThetaThresh_GUI(varargin)
% SELECT_NOTHETATHRESH_GUI MATLAB code for select_noThetaThresh_GUI.fig
%      SELECT_NOTHETATHRESH_GUI, by itself, creates a new SELECT_NOTHETATHRESH_GUI or raises the existing
%      singleton*.
%
%      H = SELECT_NOTHETATHRESH_GUI returns the handle to a new SELECT_NOTHETATHRESH_GUI or the handle to
%      the existing singleton*.
%
%      SELECT_NOTHETATHRESH_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECT_NOTHETATHRESH_GUI.M with the given input arguments.
%
%      SELECT_NOTHETATHRESH_GUI('Property','Value',...) creates a new SELECT_NOTHETATHRESH_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before select_noThetaThresh_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to select_noThetaThresh_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help select_noThetaThresh_GUI

% Last Modified by GUIDE v2.5 19-Nov-2019 17:59:12

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @select_noThetaThresh_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @select_noThetaThresh_GUI_OutputFcn, ...
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


% --- Executes just before select_noThetaThresh_GUI is made visible.
function select_noThetaThresh_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to select_noThetaThresh_GUI (see VARARGIN)

% Choose default command line output for select_noThetaThresh_GUI
handles.output = hObject;

if length(varargin)==3
    LFPbySweep=varargin{1};
    thetaDiff=varargin{2};
    Fs=varargin{3};
end

[n,x]=histcounts(nanmean(thetaDiff,2),100);
h=plot(nanmean([x(1:end-1); x(2:end)],1),n);

handles.LFPbySweep=LFPbySweep{1};;
handles.thetaDiff=thetaDiff;
handles.thetaDiffPerTrial=nanmean(thetaDiff,2);
handles.h=h;
handles.hist_x=nanmean([x(1:end-1); x(2:end)],1);
handles.hist_y=n;
handles.f=[];
handles.l=[];
handles.Fs=Fs;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes select_noThetaThresh_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = select_noThetaThresh_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

clickLoc=ginput(1);

hist_x=handles.hist_x;
hist_y=handles.hist_y;
thetaDiffPerTrial=handles.thetaDiffPerTrial;
LFPbySweep=handles.LFPbySweep;
f=handles.f;
l=handles.l;
Fs=handles.Fs;

delete(l);
close(f);

disp(clickLoc(1));

hold on;
l=line([clickLoc(1) clickLoc(1)],[0 nanmax(hist_y)],'Color','r');

isLower=find(thetaDiffPerTrial<clickLoc(1) & thetaDiffPerTrial>clickLoc(1)-0.02);
isHigher=find(thetaDiffPerTrial>clickLoc(1) & thetaDiffPerTrial<clickLoc(1)+0.02);

lo=randsample(isLower,1);
hi=randsample(isHigher,1);

i=1;
if ~isempty(lo)
    f(i)=figure();
    plot(0:1/Fs:(length(LFPbySweep(lo,:))-1)*(1/Fs),smooth(LFPbySweep(lo,:),200));
    title('Low theta');
    i=i+1;
end
if ~isempty(hi)
    f(i)=figure();
    plot(0:1/Fs:(length(LFPbySweep(hi,:))-1)*(1/Fs),smooth(LFPbySweep(hi,:),200));
    title('High theta');
    i=i+1;
end

handles.l=l;
handles.f=f;

% Update handles structure
guidata(hObject, handles);

