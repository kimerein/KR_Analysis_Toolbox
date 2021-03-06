function varargout = KR_analysisGUI(varargin)
% KR_ANALYSISGUI M-file for KR_analysisGUI.fig
%      KR_ANALYSISGUI, by itself, creates a new KR_ANALYSISGUI or raises the existing
%      singleton*.
%
%      H = KR_ANALYSISGUI returns the handle to a new KR_ANALYSISGUI or the handle to
%      the existing singleton*.
%
%      KR_ANALYSISGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in KR_ANALYSISGUI.M with the given input arguments.
%
%      KR_ANALYSISGUI('Property','Value',...) creates a new KR_ANALYSISGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before KR_analysisGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to KR_analysisGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help KR_analysisGUI

% Last Modified by GUIDE v2.5 03-Nov-2011 17:07:18

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @KR_analysisGUI_OpeningFcn, ...
                   'gui_OutputFcn',  @KR_analysisGUI_OutputFcn, ...
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


function KR_analysisGUI_OpeningFcn(hObject, eventdata, handles, varargin)

% Set expt struct
expt = varargin{1};
handles.saveExptFcn = varargin{2};

% Add analyis fields, if missing
analysisType = {'overview','orientation','contrast','srf','other'};
handles.exptModified = 0;
for i = 1:length(analysisType)
    if ~isfield(expt.analysis,analysisType{i})
        expt = analysis_def(expt,{analysisType{i}});
        handles.exptModified = 1;
    end
end

if handles.exptModified
    assignin('base','expt',expt);
end

% Update GUI
updateAnalysisGUI(handles,expt);

guidata(hObject, handles);



function varargout = KR_analysisGUI_OutputFcn(hObject, eventdata, handles) 

varargout{1} = handles.exptModified;
varargout{2} = handles.analysisGui;


function updateAnalysisGUI(handles,expt)
tagLabels = {'ov','ori','rg','srf','other'};
textEdits = {'FilesEdit','TypeEdit','TagsEdit','ValuesEdit'};
for i = 5:length(tagLabels)
analysisType = {'overview','orientation','contrast','srf','other'};
    
    for j = 1:length(textEdits)
        tempField = [tagLabels{i} textEdits{j}];
        switch textEdits{j}
            case 'FilesEdit'
                set(handles.(tempField),'String',num2str(expt.analysis.(analysisType{i}).fileInd))
            case {'TypeEdit'}
                str = parseArray(expt.analysis.(analysisType{i}).cond.type);
                set(handles.(tempField),'String',str);
            case 'TagsEdit'
                str = parseArray(expt.analysis.(analysisType{i}).cond.tags);
                set(handles.(tempField),'String',str);
            case 'ValuesEdit'
                str = parseArray(expt.analysis.(analysisType{i}).cond.values);
                set(handles.(tempField),'String',str);
        end
    end
end

setappdata(handles.analysisGui,'expt',expt);

function str = parseArray(array)

str = [];
if iscell(array)
    for i = 1:length(array)
        temp = array{i};
        if isnumeric(temp)
            temp = num2str(temp);
        end
        if i < length(array)
            str = [str temp ';' ' '];
        elseif i == length(array)
            str = [str temp];
        end
    end
else
    str = array;
end

function h = getEditVal(h)
expt = getappdata(h.analysisGui,'expt');
a = expt.analysis;
analysisType = {'overview','orientation','contrast','srf','other'};
tagLabels = {'ov','ori','rg','srf','other'};
textEdits = {'FilesEdit','TypeEdit','TagsEdit','ValuesEdit'};

for i = 5:length(analysisType)
    for j = 1:length(textEdits);
        tempField = [tagLabels{i} textEdits{j}];
        temp = get(h.(tempField),'String');
        temp = parseString(temp);
        switch textEdits{j}
            case 'FilesEdit'
                a.(analysisType{i}).fileInd = str2num(temp);
            case 'TypeEdit'
                a.(analysisType{i}).cond.type = temp;
            case 'TagsEdit'
                if ~iscell(temp)
                    temp2{1} = temp;
                    temp = temp2;
                end
                a.(analysisType{i}).cond.tags = temp;
            case 'ValuesEdit'
                if iscell(temp)
                    for k = 1:length(temp)
                        temp{k} = str2num(temp{k});
                    end
                else
                    temp2{1} = str2num(temp);
                    temp = temp2;
                end
                a.(analysisType{i}).cond.values = temp;
        end     
    end  
end


expt.analysis = a;

setappdata(h.analysisGui,'expt',expt)

function array = parseString(str)

k = strfind(str,';');

if ~isempty(k)
    array{1} = str(1:k-1);
    for i = 1:length(k)
        if i < length(k)
            temp = str((k(i)+1):(k(i+1)-1));
            array{i+1} = strtrim(temp);
        elseif i == length(k)
            temp = str((k(i)+1):end);
            array{i+1} = strtrim(temp);
        end
    end
    
else
    array = str;
end


% Files button
function srfFilesButton_Callback(hObject, eventdata, handles)
fileButtonCallback(hObject,handles)
function otherFilesButton_Callback(hObject, eventdata, handles)
fileButtonCallback(hObject,handles)
function oriFilesButton_Callback(hObject, eventdata, handles)
fileButtonCallback(hObject,handles)
function ovFilesButton_Callback(hObject, eventdata, handles)
fileButtonCallback(hObject,handles)
function rgFilesButton_Callback(hObject, eventdata, handles)
fileButtonCallback(hObject,handles)

function fileButtonCallback(hObject,handles)
expt = getappdata(handles.analysisGui,'expt');
[fileInd ok] = listdlg('ListString',{expt.files.names{:}},'ListSize',[225 300],...
    'Name','Choose files for analysis');
if ok
    tag = get(hObject,'Tag');
    temp = findstr(tag,'FilesButton');
    tempField = [tag(1:temp-1) 'FilesEdit']; 
    set(handles.(tempField),'String',mat2str(fileInd))
end



function makeFigButton_Callback(hObject, eventdata, handles)

% Get expt struct from base
expt = evalin('base','expt');

% Update expt in analysisGui appdata
setappdata(handles.analysisGui,'expt',expt);

% Get analysis type
analysisType = getAnalysisType(handles);
str = setAnalysisTypeStr(analysisType);

% Get checkbox flags (e.g. save, print, etc.)
b = getCheckVal(handles);

% Get files list
files = get(handles.([str 'FilesEdit']),'String');
if isempty(files)
    [fileInd ok] = listdlg('ListString',{expt.files.names{:}},'ListSize',[225 300],...
        'Name','Choose files for analysis');
else
    fileInd = str2num(files);
    ok = 1;
end

if ok
    if ~strcmp(analysisType,'overview') && ~b.chkVal.MUAanalysis
        % Get list of units to analyze
        if ~isfield(expt.sort,'unittags')
            disp('Can''t find any sorted units. Have you sorted any units yet?');
            return
        end
        unitList = expt.sort.unittags;
        % KR adds labels to units here
        showList={};
        for i=1:length(unitList)
            unitTag=unitList{i};
            [trodeNum unitInd] = readUnitTag(unitTag);
            if ~isempty(unitInd)
                label = getUnitLabel(expt,trodeNum,unitInd);
                showList{i}=[unitTag ' ' label];
            else
                showList{i}=unitTag;
            end
        end
        [unitInd, ok] = listdlg('ListString',showList,'ListSize',[175 300],...
            'Name','Choose units for analysis');
        unitList = unitList(unitInd);
    else
        if strcmp(analysisType,'overview')
            unitList = NaN;
        elseif b.chkVal.MUAanalysis
            disp('MUA analysis will take all spikes on the selected trodes.');
            % Get list of trodes to include in analysis
            trodeList={};
            for i=1:length(expt.sort.trode)
                trodeList{i}=expt.sort.trode(i).name;
            end 
            [trodeNums ok]=listdlg('ListString',trodeList,'ListSize',[100 200],...
                'Name','Choose trodes for analysis');
            unitList=trodeList(trodeNums);
        end
    end
    
    % Kim's analysis figs allow the combination of multiple units into
    % single fig
    % For non-"other" fig. types (e.g., orientation), separate units as
    % default, but for "other" fig. type, combine units as default
    separateUnits=1;
    sepUnits=get(handles.combineUnitsEdit,'String');
    if isequal(analysisType,'other')
        if isempty(expt.analysis.other.cond.values)
            % have not yet set expt.analysis.contrast.cond
            disp('Remember to save expt struct after entering analysis params and before analyzing!');
        return
        end
        if strcmp(sepUnits,'separate')
        %if isequal(expt.analysis.other.cond.type(end-13:end),'separate units')
            separateUnits=1;
        elseif strcmp(sepUnits,'combine')
            separateUnits=0;
        else
            if ~b.chkVal.MUAanalysis
                disp('As default, will combine all selected units into single figure.');
                separateUnits=0;
            end
        end
    end
    
    saveTag = get(handles.saveTag,'String');
    % Combine all spikes from these units into single spikes struct
    if separateUnits==0 || b.chkVal.MUAanalysis
        % Only makes "other" figs
        figType=get(handles.figureTypeEdit,'String');
        otherFig(expt,unitList,fileInd,b.chkVal,saveTag,figType);
    else % Make multiple figures, one for each unit
        for i = 1:length(unitList)
            switch analysisType
                
                case 'overview'
                    fileInd = makeFileIndStruct(handles);
                    overviewFig(expt,fileInd,b.chkVal,handles);
                case 'orientation'
                    unitTag = unitList{i};
                    [trodeNum unitInd] = readUnitTag(unitTag);
                    label = getUnitLabel(expt,trodeNum,unitInd);
                    
                    if ~strcmp('garbage',label)
                        saveTag = get(handles.saveTag,'String');
                        orientationFig(expt,unitList{i},fileInd,b.chkVal,saveTag)
                    end
                    
                case 'contrast'
                    contrastFig(expt,unitList{i},fileInd,b.chkVal)
                case 'srf'
                    srfFig(expt,unitList{i},fileInd,b.chkVal)
                    
                case 'other'
                    figType=get(handles.figureTypeEdit,'String');
                    otherFig(expt,unitList{i},fileInd,b.chkVal,saveTag,figType);
            end
        end
    end
end

function b = getCheckVal(handles)
chkNames = {'pause','close','save','print','MUAanalysis'};
for i = 1:length(chkNames)
    b.chkVal.(chkNames{i}) = get(handles.([chkNames{i} 'Chk']),'Value');
end
function str = setAnalysisTypeStr(analysisType)
switch analysisType
    case 'overview'
        str = 'ov';
    case 'orientation'
        str = 'ori';
    case 'contrast'
        str = 'rg';
    case 'srf'
        str = 'srf';
    case 'other'
        str = 'other';
end

function analysisType = getAnalysisType(handles)
% aType = {'overview','orientation','contrast','srf','other'};
% tagLabels = {'ov','ori','rg','srf','other'};
% for i = 5:length(tagLabels)
%     tempstr = [tagLabels{i} 'FigToggle'];
%     val(i) = get(handles.(tempstr),'Value');
% end
% 
% val = logical(val);
% 
% if ~any(val)
%    analysisType = '';
% else
%     analysisType = aType{val};
% end
analysisType='other';

function analysisTypeToggle(handles)


function f = makeFileIndStruct(handles)
analysisType = {'overview','orientation','contrast','srf','other'};

for i = 1:length(analysisType)
    str = setAnalysisTypeStr(analysisType{i});
    f.(analysisType{i}) = str2num(get(handles.([str 'FilesEdit']),'String'));
end


function saveExptButton_Callback(hObject, eventdata, handles)
h = handles;
h = getEditVal(h);
expt = getappdata(h.analysisGui,'expt');
assignin('base','expt',expt)
saveExpt = h.saveExptFcn{1};
ExptViewerHandles = guidata(h.saveExptFcn{2});
saveExpt(h.saveExptFcn{2},[],ExptViewerHandles);

guidata(hObject,h)

% Edit boxes
function ovTypeEdit_Callback(hObject, eventdata, handles)
function ovTagsEdit_Callback(hObject, eventdata, handles)
function ovValuesEdit_Callback(hObject, eventdata, handles)
function otherTypeEdit_Callback(hObject, eventdata, handles)
function otherTagsEdit_Callback(hObject, eventdata, handles)
function otherValuesEdit_Callback(hObject, eventdata, handles)
function srfTypeEdit_Callback(hObject, eventdata, handles)
function srfTagsEdit_Callback(hObject, eventdata, handles)
function srfValuesEdit_Callback(hObject, eventdata, handles)
function rgTypeEdit_Callback(hObject, eventdata, handles)
function rgTagsEdit_Callback(hObject, eventdata, handles)
function rgValuesEdit_Callback(hObject, eventdata, handles)
function oriValuesEdit_Callback(hObject, eventdata, handles)
function oriTagsEdit_Callback(hObject, eventdata, handles)
function oriTypeEdit_Callback(hObject, eventdata, handles)
function ovFilesEdit_Callback(hObject, eventdata, handles)
function otherFilesEdit_Callback(hObject, eventdata, handles)
function srfFilesEdit_Callback(hObject, eventdata, handles)
function oriFilesEdit_Callback(hObject, eventdata, handles)
function rgFilesEdit_Callback(hObject, eventdata, handles)

% Checkboxes
function pauseChk_Callback(hObject, eventdata, handles)
function printChk_Callback(hObject, eventdata, handles)
function closeChk_Callback(hObject, eventdata, handles)
function saveChk_Callback(hObject, eventdata, handles)

% Make figure toggles
function ovFigToggle_Callback(hObject, eventdata, handles)
function rgFigToggle_Callback(hObject, eventdata, handles)
function oriFigToggle_Callback(hObject, eventdata, handles)
function srfFigToggle_Callback(hObject, eventdata, handles)
function otherFigToggle_Callback(hObject, eventdata, handles)

function exptTableFigButton_Callback(hObject, eventdata, handles)
expt = evalin('base','expt');
bsave = get(handles.saveChk,'Value');
exptTableFig(expt,bsave)

function pushbutton10_Callback(hObject, eventdata, handles)

% --- Create function --- %
function otherTagsEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function otherTypeEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ovValuesEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ovTagsEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ovTypeEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function ovFilesEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function otherValuesEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function srfFilesEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function srfTypeEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function srfTagsEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function srfValuesEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function rgFilesEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function rgTypeEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function oriValuesEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function rgValuesEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function rgTagsEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function oriTagsEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function oriTypeEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function oriFilesEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function otherFilesEdit_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function saveTag_Callback(hObject, eventdata, handles)
% hObject    handle to saveTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of saveTag as text
%        str2double(get(hObject,'String')) returns contents of saveTag as a double


% --- Executes during object creation, after setting all properties.
function saveTag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to saveTag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% KR - show list of possible figure types and have user choose
% which figure to make

% Possible figure types
figureTypes={'KR_basicFiringRateOverTime','KR_basicFR_6layout','KR_orientation','KR_contrast','KR_separateByAllStimuli','KR_acrossChannelsFR','KR_acrossChannelsFR_withoutWindowsAnalysis','KR_acrossChannels_multiCountSpikes','SRO_overview','SRO_orientation','SRO_contrast','SRO_srf','SRO_other'};

[figType ok] = listdlg('ListString',figureTypes,'ListSize',[225 300],...
    'Name','Choose type of figure to make.');
if ok
    tag = get(hObject,'Tag');
%     temp = findstr(tag,'FilesButton');
%     tempField = [tag(1:temp-1) 'FilesEdit']; 
    %set(handles.figureTypeEdit,'String',mat2str(figType));
    set(handles.figureTypeEdit,'String',figureTypes{figType});
end
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function combineUnitsEdit_Callback(hObject, eventdata, handles)
% hObject    handle to combineUnitsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of combineUnitsEdit as text
%        str2double(get(hObject,'String')) returns contents of combineUnitsEdit as a double


% --- Executes during object creation, after setting all properties.
function combineUnitsEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to combineUnitsEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function figureTypeEdit_Callback(hObject, eventdata, handles)
% hObject    handle to figureTypeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of figureTypeEdit as text
%        str2double(get(hObject,'String')) returns contents of figureTypeEdit as a double


% --- Executes during object creation, after setting all properties.
function figureTypeEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figureTypeEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in MUAanalysisChk.
function MUAanalysisChk_Callback(hObject, eventdata, handles)
% hObject    handle to MUAanalysisChk (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of MUAanalysisChk
