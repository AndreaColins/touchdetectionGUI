function varargout = touch_gui(varargin)
% TOUCH_GUI MATLAB code for touch_gui.fig
% GUI for viewing and editing potential whisker-object touches in tracked
% whisker videos
%
% M.Evans 18.09.15
%
%      This GUI expects whisker tracking to have already been done,
%      resulting in either a .tr or _clean.mat file. Also pole
%      tracking/distance measurement should have been done already, though
%      this can be re-run from withing touch_gui
%
%
%
%      TOUCH_GUI, by itself, creates a new TOUCH_GUI or raises the existing
%      singleton*.
%
%      H = TOUCH_GUI returns the handle to a new TOUCH_GUI or the handle to
%      the existing singleton*.
%
%      TOUCH_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in TOUCH_GUI.M with the given input arguments.
%
%      TOUCH_GUI('Property','Value',...) creates a new TOUCH_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before touch_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to touch_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help touch_gui

% Last Modified by GUIDE v2.5 24-Jun-2016 13:50:13

%% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @touch_gui_OpeningFcn, ...
    'gui_OutputFcn',  @touch_gui_OutputFcn, ...
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

% --- Executes just before touch_gui is made visible.
function touch_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to touch_gui (see VARARGIN)

% Choose default command line output for touch_gui
handles.output = hObject;

% Default params for touch detection
handles.dist_default = 20; % One number per pole position/protraction or retraction, but just one initially.
handles.dist_thresh = handles.dist_default;
handles.angle_default = [10]; % plus or minus of angle from base to pole
handles.angle_thresh = handles.angle_default;
handles.kappa_default = 0.001; % Can ultimately be specific to each pole position/protraction or retraction.
handles.kappa_thresh = handles.kappa_default;
handles.trialtype_default = 1; % Read from .xlsx good_trials file
handles.trialtype = handles.trialtype_default;

% TO DO: Add in option to load working defaults from previous file

handles.dist_checkbox_val = 1;
handles.angle_checkbox_val = 0;
handles.kappa_checkbox_val = 0;
handles.trialtype_checkbox_val = 1;

% Default is for dist slider and trialtype slider to be checked
set(handles.dist_checkbox,'Value',1);
set(handles.trialtype_checkbox,'Value',1);

set(handles.dist_text,'String',num2str(handles.dist_default));
% set(handles.angle_text,'String',num2str(handles.angle_default));
% set(handles.kappa_text,'String',num2str(handles.kappa_default));
set(handles.trialtype_text,'String',num2str(handles.trialtype_default));
set(handles.choice_text,'String',num2str(handles.trialtype_default));

% Set slider values to defaults
set(handles.dist_slider,'Value',handles.dist_default)
% set(handles.angle_slider,'Value',handles.angle_default)
% set(handles.kappa_slider,'Value',handles.kappa_default)

% Initialise slider values
handles.dist_slider_val = handles.dist_default;

set(handles.text1,'BackgroundColor',[0.3,0.3,0.3])
set(handles.text1,'String','Starting up')


% Set zoom button default
set(handles.zoom_button,'Value',1)
handles.zoom = 1;

% Set plot whisker button default
set(handles.plot_whisker_tracking,'Value',0)
handles.plot_whisker = 0;

set(handles.threshold_text,'BackgroundColor',[0.3,0.3,0.3])
set(handles.threshold_text,'String','DEFAULT THRESHOLD')


% Set up key for params file
handles.touch_params_key = {'File ID';'Trial';'Trial type';'Touches detected';'Pole dist measured';'Dist thresh';'Angle thresh';'Kappa thresh';'Dist checkbox';'Angle checkbox';'Kappa checkbox';'Trial type checkbox';'Autotracked Y/N';'Choice';'Protraction 1st touch'}



% Update handles structure
guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

%% FILE MENU STUFF --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
    ['Close ' get(handles.figure1,'Name') '...'],...
    'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)

% --- Outputs from this function are returned to the command line.
function varargout = touch_gui_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
varargout{2} = handles;


%% CUSTOM FUNCTIONS START HERE

%%%% LOAD FILES FOR VIEWING
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end

ff = [dir('*.dat')];%,dir('*.avi')];% 31.03.16 limiting to .dat files, as only they work properly in matlab, and to save duplication
% ;dir('*.avi')]; % TO DO - List all files but only allow tracking of _clean.mat files


nr_files=size(ff,1);
string_list=cell(nr_files+1,1);
string_list{1}='Choose File';
for i=1:nr_files
    string_list{i+1}=ff(i).name;
end

% Load up excel file to determine trial type for all files
if ismac
    if exist('good_trials.xls','file');
        xlsfile = 'good_trials.xls';
    else
        error('xlsx file must be converted to Excel 98 xls file when working on Mac OSX')
    end
else
    xlsfile = 'good_trials.xlsx'; % Must be saved as .xls Excel 98 format for compatability
    
end
xls_info = xlsread(xlsfile);
handles.xls_info = xls_info;

% LOAD/CREATE touch_params.xls file
if exist('touch_params.csv','file');
    touchcsv = 'touch_params.csv';
    touch_params = csvread(touchcsv);
    
else
    % TO DO: NEED A WAY TO ONLY STORE PARAMS FOR VALID FILES. (?) e.g Specify only unique names (i.e. one per video)
    touch_params = zeros(nr_files,14); 
    
    
    for i = 2:numel(string_list)
        trial = str2double(string_list{i}(end-9:end-4));
        trial_num = find(xls_info(:,1) == trial);
        trialtype = xls_info(trial_num,204);
        choice = xls_info(trial_num,205);
        
        
        yes_touch_var = 0;
        yes_touch_file = 0;
        % Check for _touch.mat file
        fn = [string_list{i}(1:end-4)];
        if exist([fn,'_touch.mat'],'file')
            yes_touch_file = 1
            % Check for 'touches' variable
            touch_vars = whos('-file',[fn,'_touch.mat']);
            for j = 1:numel(touch_vars)
                if strcmp(touch_vars(j).name,'touches');
                    yes_touch_var = 1;
                end
            end
            
        else
            yes_touch_file = 0;
        end
        
        % Update touch_params
        touch_params(i-1,1) = i-1;
        touch_params(i-1,2) = trial;
        if trialtype == 1|2|3
            touch_params(i-1,3) = trialtype;
        else
            touch_params(i-1,3) = 0;
        end
        
        touch_params(i-1,4) = yes_touch_var;
        touch_params(i-1,5) = yes_touch_file;
        if choice == 1|2|3
            touch_params(i-1,14) = choice;
        else
            touch_params(i-1,14) = 0;
        end
        
    end
    
    disp('Writing params file to disk');
    dlmwrite('touch_params.csv',touch_params,'precision',6);
end

% Modify colour + font weight based on trial type and whether files
% has been tracked or not.

% Loop to format string list into HTML with appropriate colour/weight tags
html_string_list{1} = ['<HTML><FONT COLOR="black"<b>',string_list{1},'</b></HTML>'];
string_list;

handles.ttt = zeros(1,3); % number of each trial type that has been tracked.
handles.ttu = zeros(1,3); % Total number for each trial type

handles.ttu = hist(touch_params(:,3),linspace(0.5,3.5,3));
tracked = find(touch_params(:,4));
handles.ttt = hist(touch_params(tracked,3),linspace(0.5,3.5,3));

for i = 2:numel(string_list)
    yes_touch_var = touch_params(i-1,4);
    trialtype = touch_params(i-1,3);
    % If there is an existing touches variable, make text bold
    if yes_touch_var
        
        if trialtype == 1
            html_string_list{i} = ['<HTML><FONT COLOR="lime"<b>',string_list{i},'</b></HTML>'];
            %             handles.ttt(1) = handles.ttt(1) + 1;
            %             handles.ttu(1) = handles.ttu(1) + 1;
        elseif trialtype == 2
            html_string_list{i} = ['<HTML><FONT COLOR="red"<b>',string_list{i},'</b></HTML>'];
            %             handles.ttt(2) = handles.ttt(2) + 1;
            %             handles.ttu(2) = handles.ttu(2) + 1;
        elseif trialtype == 3
            html_string_list{i} = ['<HTML><FONT COLOR="black"<b>',string_list{i},'</b></HTML>'];
            %             handles.ttt(3) = handles.ttt(3) + 1;
            %             handles.ttu(3) = handles.ttu(3) + 1;
        else
            disp('Unknown trial type')
            html_string_list{i} = ['<HTML><FONT COLOR="white"<b>',string_list{i},'</b></HTML>'];
        end
        
        % Otherwise use regular text
    else
        if trialtype == 1
            html_string_list{i} = ['<HTML><FONT COLOR="lime">',string_list{i},'</HTML>'];
            %             handles.ttu(1) = handles.ttu(1) + 1;
        elseif trialtype == 2
            html_string_list{i} = ['<HTML><FONT COLOR="red">',string_list{i},'</HTML>'];
            %             handles.ttu(2) = handles.ttu(2) + 1;
        elseif trialtype == 3
            html_string_list{i} = ['<HTML><FONT COLOR="black">',string_list{i},'</HTML>'];
            %             handles.ttu(3) = handles.ttu(3) + 1;
        else
            disp('Unknown trial type')
            html_string_list{i} = ['<HTML><FONT COLOR="white">',string_list{i},'</HTML>'];
        end
    end
end


% TO DO: add option to run contact_detection on all untracked files
handles.html_string_list = html_string_list;
handles.string_list = string_list;
set(hObject,'String',html_string_list);
guidata(hObject, handles);


% --- DISPLAY FRAME OF VIDEO + TRACKED OUTPUTS
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1

val = get(hObject,'Value');
if(val~=1)
    %     string_list = get(hObject,'String');
    %     handles.fname = string_list{val};   % avi file or dat file.
    handles.fname = handles.string_list{val};
else
    error('')
end
handles.fid = val-1;

handles = update_ttt(handles);


handles = initialise_new_video(handles);

guidata(hObject, handles);



%% Button press functions

% --- Executes on button press in reset_defaults.
function reset_defaults_Callback(hObject, eventdata, handles)
% hObject    handle to reset_defaults (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = reset_everything(handles);

guidata(hObject, handles);

% --- Executes on button press in save_output_button.
function save_output_button_Callback(hObject, eventdata, handles)
% hObject    handle to save_output_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = save_everything(handles);

guidata(hObject, handles);


% --- Executes on button press in plot_whisker_tracking.
function plot_whisker_tracking_Callback(hObject, eventdata, handles)
% hObject    handle to plot_whisker_tracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.plot_whisker = get(hObject,'Value');
handles = update_frame(handles);

guidata(hObject,handles);



% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% NOTE THIS DOESN'T WORK QUITE AS NICELY AS IT COULD/SHOULD 30.09.15
%
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

vpos = get(handles.output_figure,'position'); % origin is upper left
pt = get(hObject,'CurrentPoint');    % mouse click location wrt figure

% location wrt output_figure axes:
tpt(1) = pt(1)-vpos(1);
tpt(2) = vpos(4)-(pt(2)-vpos(2));

% if
if tpt(1) >=0 && tpt(1)<=1 && tpt(2) >=0 && tpt(2)<=1
    %     datacursormode on
    dcm = datacursormode;
    cursor_pos = getCursorInfo(dcm);
    %     set(dcm,'UpdateFcn',@cursorposition,'SnapToDataVertex','on');
    
    % cursor_pos = get(0,'userdata')
    cursor_pos.Position(1)
    handles.CursorPos = cursor_pos.Position(1);
    
    handles.frameidx = handles.CursorPos;
    
    handles = update_frame(handles);
    datacursormode off
end
guidata(hObject,handles);



% --- PLOT TOUCHES BASED ON SELECTIONS BELOW
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

set(handles.text1,'BackgroundColor',[1,0,0])
set(handles.text1,'String','FINDING TOUCHES')

handles = find_touches(handles);

% Detect whether first touch is protraction/retraction
handles = protraction_detector(handles);

set(handles.text1,'BackgroundColor',[0,1,0])
set(handles.text1,'String','Touches found')

handles = update_frame(handles);

handles = plot_touches(handles);



guidata(hObject, handles);





%% Control video frame playback

% --- Executes on button press in back_1_frame.
function back_1_frame_Callback(hObject, eventdata, handles)
% hObject    handle to back_1_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.frameidx = max([handles.frameidx-1,1]);
% TO DO ADD SKIP TO NEXT TOUCH
handles = update_frame(handles);
handles = plot_touches(handles);

guidata(hObject,handles)

% --- Executes on button press in forward_1_frame.
function forward_1_frame_Callback(hObject, eventdata, handles)
% hObject    handle to forward_1_frame (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.frameidx = min([handles.frameidx+1,handles.nframes]);
% TO DO ADD SKIP TO NEXT TOUCH
handles = update_frame(handles);
handles = plot_touches(handles);

guidata(hObject,handles)




% --- Executes on button press in touch_yes_no.
function touch_yes_no_Callback(hObject, eventdata, handles)
% hObject    handle to touch_yes_no (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = toggle_touch(handles);

guidata(hObject,handles)


%% CHECKBOX AND SLIDER STUFF
% --- Executes on button press in dist_checkbox.
function dist_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to dist_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dist_checkbox
handles.dist_checkbox_val = get(hObject,'Value');
handles = plot_touches(handles);

guidata(hObject, handles);

% --- Executes on slider movement.
function dist_slider_Callback(hObject, eventdata, handles)
% hObject    handle to dist_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
handles.dist_slider_val = get(hObject,'Value');

if handles.dist_checkbox_val;
    handles.dist_thresh = handles.dist_slider_val;
else
    handles.dist_thresh = handles.dist_default;
end
set(handles.threshold_text,'BackgroundColor',[1,0.3,0.3])
set(handles.threshold_text,'String','USING CUSTOM THRESHOLD')

% set(handles.dist_text,'BackgroundColor',[0.5,0.5,0.5]);
set(handles.dist_text,'String',num2str(handles.dist_slider_val));
handles = plot_touches(handles);
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function dist_slider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dist_slider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
guidata(hObject,handles);


% --- Executes on button press in angle_checkbox.
function angle_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to angle_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of angle_checkbox
handles.angle_checkbox_val = get(hObject,'Value');
handles = plot_touches(handles);
guidata(hObject,handles);


% --- Executes on button press in data_driven.
function data_driven_Callback(hObject, eventdata, handles)
% hObject    handle to data_driven (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles = data_driven_function(handles);

guidata(hObject,handles);

% Function version to allow calling from other functions.
function handles = data_driven_function(handles)
% V1: set threshold to max of currently extracted touches
handles.dist_thresh = max(handles.pole_dist(find(handles.touches)));
% TO DO add thresholds for other variables

% TO DO add different threshold for protraction/retraction
set(handles.threshold_text,'BackgroundColor',[1,0.3,0.3])
set(handles.threshold_text,'String','USING DATA DERIVED THRESH')
set(handles.text1,'BackgroundColor',[1,0,0])
set(handles.text1,'String','FINDING TOUCHES')
handles = find_touches(handles);
set(handles.text1,'BackgroundColor',[0,1,0])
set(handles.text1,'String','Touches found')

% Detect whether first touch is protraction/retraction
handles = protraction_detector(handles);

handles = plot_touches(handles);
% Replot image just in case the touch toggle changes
handles = update_frame(handles);


% set(handles.dist_text,'BackgroundColor',[0.5,0.5,0.5]);
set(handles.dist_text,'String',num2str(handles.dist_thresh));

% --- Executes on button press in kappa_checkbox.
function kappa_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to kappa_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of kappa_checkbox
handles.kappa_checkbox_val = get(hObject,'Value');
handles = plot_touches(handles);
guidata(hObject,handles);

% --- Executes on button press in trialtype_checkbox.
function trialtype_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to trialtype_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of trialtype_checkbox

handles.trialtype_checkbox_val = get(hObject,'Value');
handles = plot_touches(handles);
guidata(hObject,handles);

% % --- Executes on slider movement.
% function angle_slider_Callback(hObject, eventdata, handles)
% % hObject    handle to angle_slider (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'Value') returns position of slider
% %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% handles.angle_slider_val = get(hObject,'Value');
% 
% if handles.angle_checkbox_val;
%     handles.angle_thresh = handles.angle_slider_val;
% else
%     handles.angle_thresh = handles.angle_default;
% end
% set(handles.threshold_text,'BackgroundColor',[1,0.3,0.3])
% set(handles.threshold_text,'String','USING CUSTOM THRESHOLD')

% % set(handles.dist_text,'BackgroundColor',[0.5,0.5,0.5]);
% % set(handles.angle_text,'String',num2str(handles.angle_slider_val));
% handles = plot_touches(handles);
% guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
% function angle_slider_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to angle_slider (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: slider controls usually have a light gray background.
% if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor',[.9 .9 .9]);
% end
% guidata(hObject,handles);


% % --- Executes on slider movement.
% function kappa_slider_Callback(hObject, eventdata, handles)
% % hObject    handle to kappa_slider (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% 
% % Hints: get(hObject,'Value') returns position of slider
% %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% handles = plot_touches(handles);
% guidata(hObject, handles);

% % --- Executes during object creation, after setting all properties.
% function kappa_slider_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to kappa_slider (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: slider controls usually have a light gray background.
% if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor',[.9 .9 .9]);
% end
% guidata(hObject,handles);


%% PRIVATE FUNCTIONS
function handles = initialise_new_video(handles)

set(handles.text1,'BackgroundColor',[0.5,0.5,0])
set(handles.text1,'String','LOADING VIDEO')

% Function version of video file initialisation
handles = open_video(handles);


% Check autotracked status.
touchcsv = 'touch_params.csv';
touch_params = csvread(touchcsv);

handles.autotracked = touch_params(handles.fid,13);

% Load choice and trial type
handles.trialtype = touch_params(handles.fid,3);
handles.choice = touch_params(handles.fid,14);

% Plot only or do contact detection stuff??
handles.plot_only = 0;

% Read in first frame
handles.frameidx = 1;
axes(handles.video_figure);
p = get(gca,'position');
cla(gca);
handles.frame = load_frame(handles.video,handles.frameidx);
imagesc(handles.frame);
set(gca,'XTick',[],'YTick',[])
title(['frame: ',num2str(handles.frameidx)]);

% LOAD POLE TRACKING OUTPUT
fn = [handles.fname(1:end-4),'_touch.mat']
if exist(fn,'file')
    load(fn,'start_frame','trigger_frame','pole_dist','gof','radius','barPos','closest_w')
    handles.start_frame = start_frame;
    handles.trigger_frame = trigger_frame;
    handles.pole_dist = pole_dist;
    handles.gof = gof;
    handles.radius = radius;
    handles.barPos = barPos;
    handles.closest_w = closest_w;
else
    %     error(['No _touch.mat file available. Run contact_detector first'])
    choice = questdlg('No _touch.mat file for this movie. Run contact detector here (it may be faster to run it separately) ?',...
        'Contact detector yay or nay','Yes please','No thanks','Try parallel version (experimental)','Yes please');
    switch choice
        case 'Yes please'
            set(handles.text1,'BackgroundColor',[1,1,0]);
            set(handles.text1,'String','DETECTING CONTACTS');
            contact_detector(handles.fname,14,0);
            
            set(handles.text1,'BackgroundColor',[0.5,1,0]);
            set(handles.text1,'String','LOADING POLE TRACKING OUTPUT');
            
            load(fn,'start_frame','trigger_frame','pole_dist','gof','radius','barPos','closest_w');
            handles.start_frame = start_frame;
            handles.trigger_frame = trigger_frame;
            handles.pole_dist = pole_dist;
            handles.gof = gof;
            handles.radius = radius;
            handles.barPos = barPos;
            handles.closest_w = closest_w;
            
            set(handles.text1,'BackgroundColor',[0,1,0]);
            set(handles.text1,'String','DONE');
            
            % Re-initialise video
            handles = open_video(handles);
            
        case 'No thanks'
            disp(['No _touch.mat file available. Touch detection requires you run contact_detector first, or choose another movie']);
            disp(['Plotting kappa and theta only for inspection']);
            handles.plot_only = 1;
            handles.touches = zeros(1,handles.nframes);
            
        case 'Try parallel version (experimental)'
            set(handles.text1,'BackgroundColor',[0.3,0,0]);
            set(handles.text1,'String','DETECTING CONTACTS (PARALLEL)');
            contact_detector_parallel(handles.fname,14);
            
            set(handles.text1,'BackgroundColor',[0.5,1,0]);
            set(handles.text1,'String','LOADING POLE TRACKING OUTPUT');
            
            load(fn,'start_frame','trigger_frame','pole_dist','gof','radius','barPos','closest_w');
            handles.start_frame = start_frame;
            handles.trigger_frame = trigger_frame;
            handles.pole_dist = pole_dist;
            handles.gof = gof;
            handles.radius = radius;
            handles.barPos = barPos;
            handles.closest_w = closest_w;
            
            set(handles.text1,'BackgroundColor',[0,1,0]);
            set(handles.text1,'String','DONE');
            
            % Re-initialise video
            handles = open_video(handles);
    end
end

% LOAD WHISKER TRACKING OUTPUT, EITHER WHIKERMAN OR CHIMERA
set(handles.text1,'BackgroundColor',[0.5,1,0]);
set(handles.text1,'String','LOADING WHISKER TRACKING OUTPUT')

fn = [handles.fname(1:end-4)];
if exist([fn,'.tr4'],'file')
    load([fn,'.tr4'],'whisker','-mat');
    handles.r_base = whisker.r3all(:,1:2,:);
    handles.kappa = whisker.kappa_all;
    handles.theta = whisker.theta_all;
    handles.fp = whisker.fp3_all(:,1:2);
    
elseif exist([fn,'_clean.mat'],'file')
    load([fn,'_clean.mat'],'r_base','theta_w','kappa_w');
    handles.r_base = r_base;
    handles.kappa = kappa_w;
    handles.theta = theta_w;
    handles.fp = squeeze(r_base(:,:,1));
else
    error(['No .tr or _clean.mat file available. Track video with chimera.m or Whikerman before continuing'])
end

% Check if pole up times, trial type or touches exist in _touch.mat.
% Compute them if not
% TO DO: COMPUTE THIS IN A LOOP FOR ALL FILES DURING INITIALISATION,
% MAKING A NOTE OF TOUCH-DETECTED FILES FOR DISPLAY IN THE GUI

% SOME DONE ABOVE: MAY NEED FURTHER REFACTORING
if exist([fn,'_touch.mat']);
    touch_vars = whos('-file',[fn,'_touch.mat']);
    
    has_variable = zeros(1,7);
    for i = 1:numel(touch_vars)
        % LOAD POLE UP TIMES
        if strcmp(touch_vars(i).name,'pole_up_times');
            has_variable(1) = 1;
            load([fn,'_touch.mat'],'pole_up_times');
            handles.pole_up_times = pole_up_times;
        end
        % LOAD TOUCHES
        if strcmp(touch_vars(i).name,'touches');
            has_variable(2) = 1;
            load([fn,'_touch.mat'],'touches');
            handles.touches = touches;
            set(handles.threshold_text,'BackgroundColor',[1,0.3,0.3])
            set(handles.threshold_text,'String','LOADED TOUCHES FROM FILE')
            
            % Re-compute touch period onsets
            handles = touch_periods(handles);
            if handles.touch_onset
                handles.frameidx = handles.touch_onset(1);
%                 handles = update_frame(handles);
            end
        end
        % LOAD TRIALTYPE
        if strcmp(touch_vars(i).name,'trialtype');
            has_variable(3) = 1;
            load([fn,'_touch.mat'],'trialtype');
            handles.trialtype = trialtype;
        end
        
        % LOAD Protraction/Retraction touch
        if strcmp(touch_vars(i).name,'protraction_touch');
            has_variable(4) = 1;
            load([fn,'_touch.mat'],'protraction_touch');
            handles.protraction_touch = protraction_touch;
        end
        
        % LOAD first touch touch
        if strcmp(touch_vars(i).name,'first_touch');
            has_variable(5) = 1;
            load([fn,'_touch.mat'],'first_touch');
            handles.first_touch = first_touch;
        end
        
        % LOAD POLE UPDOWN
        if strcmp(touch_vars(i).name,'pole_updown');
            has_variable(6) = 1;
            load([fn,'_touch.mat'],'pole_updown');
            handles.pole_updown = pole_updown;
            
%             handles.pole_updown(1) = pole_up_times(1);
%             handles.pole_updown(2) = pole_up_times(end);
        end
        
        % CHECK EXISTENCE OF PRO_RET
        if strcmp(touch_vars(i).name,'pro_ret');
            has_variable(7) = 1;
            load([fn,'_touch.mat'],'pro_ret');
            handles.pro_ret = pro_ret;
            
%             handles.pole_updown(1) = pole_up_times(1);
%             handles.pole_updown(2) = pole_up_times(end);
        end

        %     % LOAD TOUCH PARAMS
        %     if strcmp(touch_vars(i).name,'touch_params');
        %         has_variable(4) = 1;
        %         load([fn,'_touch.mat'],'touch_params');
        %         handles.touch_params = touch_params;
        %         % TO DO. UPDATE PARAM VALUES IN GUI
        %     end
    end
    
    % If any variables are missing, compute them
    display(['has_variable = ',num2str(has_variable)])
    
    % Compute pole up times
    if has_variable(1) == 0;
        handles = pole_uptimes(handles);
    end
    % Create empty touches array
    if has_variable(2) == 0;
        handles.touches = zeros(1,numel(handles.gof));
    end
    

    % Reset some values here to prevent errors
    handles.touch_thisframe = handles.touches(handles.frameidx);
    handles.CursorPos = 1;
    
    % Load trial type from excel metadata file
    if has_variable(3) == 0;
        handles = get_trial_type(handles);
    end

    if has_variable(4) == 0;
    % Compute whether first touch is protraction/retraction
    handles = pro_ret_calculator(handles);
    handles = protraction_detector(handles);
    end
    
    if has_variable(5) == 0;
    % Compute first touch times
    handles = pro_ret_calculator(handles);
    handles = protraction_detector(handles);
    end
    
    if has_variable(6) == 0;
        % Compute times when pole is up and down
        handles = pole_uptimes(handles);
    end
    
    if has_variable(7) == 0;
    % Compute first touch times
    handles = pro_ret_calculator(handles);
    handles = protraction_detector(handles);
    end
    
    if handles.first_touch >=1
        handles.frameidx = handles.first_touch;
    end
    
    handles = update_frame(handles);
    handles = pro_ret_text_set(handles);
    

%     if has_variable(4) == 0; % 
%         handles.pro_ret = zeros(1,numel(handles.gof));
%     end
    
    % Set touch_params to defaults
    % if has_variable(4) == 0;
    %     handles = reset_everything(handles);
    %
    %     handles.touch_params = {handles.dist_thresh, ...
    %                             handles.angle_thresh, ...
    %                             handles.kappa_thresh, ...
    %                             handles.trialtype, ...
    %                             handles.dist_checkbox_val, ...
    %                             handles.angle_checkbox_val, ...
    %                             handles.kappa_checkbox_val, ...
    %                             handles.trialtype_checkbox_val};
    % end
    
    
end

clear fn


% Signify trial type in gui
set(handles.trialtype_text,'String',num2str(handles.trialtype));
tt_color = ([1,1,1; 0,1,0; 1,0,0; 0,0,0]);
set(handles.trialtype_text,'ForegroundColor',tt_color(handles.trialtype + 1,:));

% Signify choice in gui
handles.choice
set(handles.choice_text,'String',num2str(handles.choice));
set(handles.choice_text,'ForegroundColor',tt_color(handles.choice + 1,:));
% PLOT DISTANCE ETC
handles = plot_touches(handles);



set(handles.text1,'BackgroundColor',[0,1,0]);
set(handles.text1,'String','DONE')


% LOAD_FRAME
function frame = load_frame(videopmtrs,frameidx)

switch videopmtrs.type
    case 'avi'
        video = read(videopmtrs.vObj, frameidx);
        frame = video(:,:,:,1);
    case 'dat'
        offset = videopmtrs.header.imagesize * (frameidx-1) + videopmtrs.offset;
        fseek(videopmtrs.fid,offset,-1);
        tmp = fread(videopmtrs.fid,videopmtrs.header.imagesize-24,'uint8=>uint8');
        tmp = reshape([tmp; zeros(24,1)],videopmtrs.width,videopmtrs.height)';
        frame = uint8(zeros(videopmtrs.height,videopmtrs.width,3));
        frame(:,:,1) = tmp;
        frame(:,:,2) = tmp;
        frame(:,:,3) = tmp;
        clear tmp
    otherwise
        error('Unhandled video file type')
end

% WORK OUT TIMES WHEN POLE IS UP FROM GOODNESS-OF-FIT
% Also computes angle to the pole
% TO DO: save pole up times to _touch.mat
function handles = pole_uptimes(handles)
pole_up_window = handles.trigger_frame:handles.trigger_frame + 700;
pole_up_window = mod(pole_up_window,numel(handles.gof)) + 1;

% Set gof threshold to mean(gof)-4.*std(gof)
%     gof_peak_thresh = nanmean(handles.gof(pole_up_window))+3.*nanstd(handles.gof(pole_up_window));

% New method using max of gof
gof_peak_thresh = max(handles.gof)-1; % May occasionally not find the correct peaks. Increasing this threshold can help. Usually max - 4, can go up to max - 0.5 etc


%%%%%%%%%%%preprocess gof
gof2=timeseries(handles.gof,[0:1:numel(handles.gof)-1]./1000);
gof3=idealfilter(gof2,[0,5],'pass');
% Find 2 peaks in gof corresponding to pole-in-focal-plane times
[~,pole_updown,peakwidth] = findpeaks(squeeze(gof3.Data),'MinPeakWidth',500,'Npeaks',1)
% handles
if numel(pole_updown) == 1;
    disp(['file ',handles.fname(1:end-4),' only has one gof peak']);
     pole_updown(1) = pole_updown(1)-200;
    pole_updown(2) = pole_updown(1)+floor(peakwidth);
elseif numel(pole_updown) > 2;
    disp(['file ',handles.fname(1:end-4),' has more than 2 gof peaks']);
    pole_updown = pole_updown(1:2);
elseif numel(pole_updown) == 0;
    disp(['No gof peaks found in file ',handles.fname(1:end-4)]);
    pole_updown = [handles.start_frame - 1,handles.trigger_frame];
end

% Work out which order the peaks are (depends on timing of trigger in
% the video

if numel(find(pole_updown>=handles.trigger_frame)) == 2||0;
    pole_up_times = pole_updown(1):pole_updown(2);
else
    pole_up_times = [pole_updown(2):numel(handles.pole_dist),1:pole_updown(1)];
end

handles.pole_up_window = pole_up_window;
handles.gof_peak_thresh = gof_peak_thresh;
handles.pole_updown = pole_updown;
handles.pole_up_times = pole_up_times;

% Compute angle to the pole
basePos = handles.r_base(:,:,1);
touch_angle_Vec = basePos - handles.closest_w;
handles.touch_angle = atan2(touch_angle_Vec(:,1),touch_angle_Vec(:,2)) *180./pi;

% Angle of touch from pole centre. > 180 = protraction touch
cont_angle_Vec = handles.barPos - handles.closest_w;
handles.cont_angle = atan2(cont_angle_Vec(:,1),cont_angle_Vec(:,2))*180./pi;

%%clear variables
clear gof2 gof3
% Protraction vs retraction touch now done with a function
handles = pro_ret_calculator(handles);



% Extract trial type from excel file
% TO DO: FIND WAY OF DOING THIS FOR ALL FILES IN A DIRECTORY
% This should probably be happening in contact_detector.
% Save trial type to _touch.mat here too, then load later if available
function handles = get_trial_type(handles)
% Check computer type/file existance
if ismac
    if exist('good_trials.xls','file');
        xlsfile = 'good_trials.xls';
    else
        error('xlsx file must be converted to Excel 98 xls file when working on Mac OSX')
    end
else
    xlsfile = 'good_trials.xlsx'; % Must be saved as .xls Excel 98 format for compatability
    
end
xls_info = xlsread(xlsfile);
trial = str2double(handles.fname(end-9:end-4));
trial_num = find(xls_info(:,1) == trial);
handles.trialtype = xls_info(trial_num,204);
handles.trial = trial;
trialtype = handles.trialtype;
save([handles.fname(1:end-4),'_touch.mat'],'trialtype','-append');

%% Update video frame image
function handles = update_frame(handles);

axes(handles.video_figure);
p = get(gca,'position');
cla(gca);
handles.frame = load_frame(handles.video,handles.frameidx);
imagesc(handles.frame);
set(gca,'XTick',[],'YTick',[]);

if handles.plot_only == 0;
    if ismember(handles.frameidx,find(handles.touches))
        handles.touch_thisframe = handles.touches(handles.frameidx);
        title(['frame: ',num2str(handles.frameidx)],'color',[0,0,1]);
    else
        handles.touch_thisframe = 0;
        title(['frame: ',num2str(handles.frameidx)],'color',[0,0,0]);
    end
    
    if handles.zoom
        xlim([handles.barPos(handles.frameidx,1)-40,handles.barPos(handles.frameidx,1)+40])
        ylim([handles.barPos(handles.frameidx,2)-20,handles.barPos(handles.frameidx,2)+40])
    end
    
else
    title(['frame: ',num2str(handles.frameidx)],'color',[0,0,0]);
end

if handles.plot_whisker
    temp = WhikerMan_functions({squeeze(1+handles.r_base(handles.frameidx,:,:)),[0:0.01:20]},8);
    b_tmp = temp{:};
    hold on
    plot(b_tmp(1,:),b_tmp(2,:),'m','Linewidth',2);
    plot(handles.barPos(handles.frameidx,1),handles.barPos(handles.frameidx,2),'w.')
    plot(handles.closest_w(handles.frameidx,1),handles.closest_w(handles.frameidx,2),'c.')
    
end
handles.r_base(handles.frameidx,:,:);



%% Reset everything. Allows calling by button or from inside another
% function
% TO DO. PUT SOME SETUP STUFF HERE TO BE CALLED DURING SETUP...
function handles = reset_everything(handles)
set(handles.dist_checkbox,'value',1)
set(handles.angle_checkbox,'value',0);
set(handles.kappa_checkbox,'value',0);
set(handles.trialtype_checkbox,'value',0);

handles.dist_thresh = handles.dist_default;
set(handles.dist_text,'String',num2str(handles.dist_thresh));
set(handles.dist_slider,'value',0);

% handles.angle_thresh = handles.angle_default;
% set(handles.angle_text,'String',num2str(handles.angle_thresh));
% set(handles.angle_slider,'value',0);

% handles.kappa_thresh = handles.kappa_default;
% set(handles.kappa_text,'String',num2str(handles.kappa_thresh));
% set(handles.kappa_slider,'value',0);

set(handles.threshold_text,'BackgroundColor',[0.3,0.3,0.3])
set(handles.threshold_text,'String','DEFAULT THRESHOLD')

% Find touches based on current criteria.
function handles = find_touches(handles)

handles.autotracked = 0;

% TO DO. Use touch-type appropriate distance threshold
touches = find(handles.pole_dist(handles.pole_up_times)<= handles.dist_thresh);
handles.touches = zeros(1,numel(handles.gof));
handles.touches(handles.pole_up_times(touches)) = 1;

% TO DO. Add switches for different checkboxes
if handles.dist_checkbox_val
end
if handles.angle_checkbox_val
end
if handles.kappa_checkbox_val
end

handles.dist_thresh
handles.angle_thresh
handles.kappa_thresh

% Re-compute touch period edges
handles = touch_periods(handles);

% Set frame index to first touch
if handles.first_touch >=1
    handles.frameidx = handles.first_touch;
end





% PLOT_TOUCHES
function handles = plot_touches(handles)
% Trying to plot in a way that updates the plot instead of re-drawing it
axes(handles.output_figure);
cla(gca);
if handles.plot_only == 0;
    plot(handles.output_figure,1:numel(handles.pole_dist),handles.pole_dist,'color',[0.2,0.2,0.2]);
    hold on
    plot(1:numel(handles.gof),handles.gof,'color',[0,0.4470,0.7410]); % 2015 default blue
    plot(handles.pole_up_times,handles.gof(handles.pole_up_times),'r.');
    handles.legend = {'pole distance','goodness of fit','pole up'};
%     handles.legend = {'Pole distance';'goodness of fit';'pole up';'touches'};
    % TO DO. Add switches for different checkboxes
    
%     keyboard
    if handles.dist_checkbox_val
    end
    if handles.angle_checkbox_val
        % SET PLOT COLOR BY PROTRACTION/RETRACTION
        theta = 100.*zscore(handles.theta);
        plot(theta,'m');
%         plot(find(handles.pro_ret),100.*(handles.theta(find(handles.pro_ret)) - nanmean(handles.theta))./nanstd(handles.theta),'g.');
        plot(find(handles.pro_ret),theta(find(handles.pro_ret)),'g.');
        
        handles.legend = {handles.legend{:}, 'theta','protraction'};
        
%         plot(handles.cont_angle,'b')    
%         plot(zeros(size(handles.touch_angle)),'k')
        %         handles.legend = {handles.legend{:};'theta'};
    end
    if handles.kappa_checkbox_val
        plot(100.*zscore(handles.kappa),'b');
        handles.legend = {handles.legend{:},'kappa'};
        %         handles.legend = {handles.legend{:};'kappa'};
    end
    
    plot(find(handles.touches),handles.pole_dist(find(handles.touches)),'k.');
    %     handles.legend = {'Pole distance';'goodness of fit';'pole up';'theta';'kappa';'touches'};
    handles.legend = {handles.legend{:},'touches'};
    % First touch
    %     handles.first_touch
%     handles
    if numel(find(handles.first_touch)) >= 1
        plot(handles.first_touch,handles.pole_dist(handles.first_touch),'c.','markersize',10);
        %     handles.legend = {'Pole distance';'goodness of fit';'pole up';'theta';'kappa';'touches'};
        handles.legend = {handles.legend{:},'first touch'};
    end
    % TO DO: add plot of first touch
else
    plot(100.*zscore(handles.kappa),'b');
    hold on
    plot(100.*zscore(handles.theta),'r');
    handles.legend = {'kappa','theta'};
    
end
% Add vertucal line indicating current frame
plot([handles.frameidx,handles.frameidx],[-80,80],'color',[1,0.5,0],'linewidth',2)
handles.legend = {handles.legend{:},'current frame'};
%legend(handles.legend);



% Toggle touch yes/no. Can be done with spacebar or by pressing the button
function handles = toggle_touch(handles)

handles.autotracked = 0 ;

if handles.touch_thisframe
    handles.touch_thisframe = 0;
else
    handles.touch_thisframe = 1;
end

% handles.touch_thisframe = mod(handles.touch_thisframe,1);
handles.touches(handles.frameidx) = handles.touch_thisframe;

handles = update_frame(handles);

handles = plot_touches(handles);

% Re-compute touch period onsets
handles = touch_periods(handles);

% Re-compute protraction/retraction first touch
handles = protraction_detector(handles);

set(handles.threshold_text,'BackgroundColor',[1,1,0]);
set(handles.threshold_text,'String','MANUALLY CURATED TOUCHES');

set(handles.text1,'BackgroundColor',[1,0.3,0]);
set(handles.text1,'String','UNSAVED MANUALLY CURATED TOUCHES');

% Determine edges of touch periods for faster navigation
function handles = touch_periods(handles);
handles.touch_onset = find(diff(handles.touches)>=1) + 1;
handles.touch_offset = find(diff(handles.touches)<=-1);
handles.touch_on_off = sort([handles.touch_onset,handles.touch_offset]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in zoom_button.
function zoom_button_Callback(hObject, eventdata, handles)
% hObject    handle to zoom_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.zoom = get(hObject,'Value');
axes(handles.video_figure);
if handles.zoom
    xlim([handles.barPos(handles.frameidx,1)-40,handles.barPos(handles.frameidx,1)+40]);
    ylim([handles.barPos(handles.frameidx,2)-20,handles.barPos(handles.frameidx,2)+40]);
else
    roi = size(handles.frame);
    xlim([0.5,roi(2)+0.5]);
    ylim([0.5,roi(1)+0.5]);
end
guidata(hObject,handles);


% --- Executes on button press in load_thresh.
function load_thresh_Callback(hObject, eventdata, handles)
% hObject    handle to load_thresh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% TO DO. ADD CODE TO SAVE THRESHOLDS TO METADATA FILE, THEN LOAD THEM HERE
guidata(hObject,handles);

% --- Executes on button press in contact_detection.
function contact_detection_Callback(hObject, eventdata, handles)
% hObject    handle to contact_detection (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.text1,'BackgroundColor',[1,1,0]);
set(handles.text1,'String','DETECTING CONTACTS')

choice = questdlg('Which pole radius should I use ?',...
    'Pole radius Q&A, for better gof peak identification.','3','14 (Default)','23','14 (Default)');
switch choice
    case '3'
        handles.radius = 3;
    case '14 (Default)'
        handles.radius = 14;
    case '23'
        handles.radius = 23;
   
end

choice = questdlg('Contact detection for all files in the directory or just this one?',...
    'Contact detection of all files or just one','All','Just this one (Default)','Just this one (Default)');
switch choice
    case 'All'
        for val = 2:numel(handles.string_list)
            handles.fname = handles.string_list{val};
            contact_detector(handles.fname,handles.radius,0);
        end
    case 'Just this one (Default)'
        handles.fname
        contact_detector(handles.fname,handles.radius,0); % Now using whole file name
end
set(handles.text1,'BackgroundColor',[1,1,0]);
set(handles.text1,'String','Idle.')
guidata(hObject,handles);


% --- Executes on key press with focus on figure1 and none of its controls.
function figure1_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.FIGURE)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

eventdata.Key

frameidx = handles.frameidx;

change_frame = 0;
toggle = 0;
pro_ret_toggle = 0;

if length(eventdata.Modifier) == 0
    switch eventdata.Key
        case 'leftarrow'
            frameidx = max([frameidx-1,1]); change_frame = 1;
            handles.frameidx = frameidx;
        case 'rightarrow'
            frameidx = min([frameidx+1,handles.nframes]); change_frame = 1;
            handles.frameidx = frameidx;
            
            % Can toggle touch with up, down, or space
        case 'uparrow'
            toggle = 1;
            handles.frameidx = frameidx;
        case 'downarrow'
            toggle = 1;
            handles.frameidx = frameidx;
        case 'space'
            toggle = 1;
            handles.frameidx = frameidx;
    end
    % SKIP TO NEXT TOUCH IF SHIFT IS HELD
    % TO DO: MAKE INDEX START AFTER POLE UP
    
elseif strmatch(eventdata.Modifier{1},'shift')
    next_tch = find(handles.touch_on_off > handles.frameidx);
    if next_tch
        next_tch = handles.touch_on_off(next_tch(1));
    else
        next_tch = handles.frameidx
        disp('Already at last touch')
    end
    
    prev_tch = find(handles.touch_on_off < handles.frameidx);
    if prev_tch
        prev_tch = handles.touch_on_off(prev_tch(end));
    else
        prev_touch = frameidx
        disp('Already at first touch')
    end
    
    switch eventdata.Key
        case 'leftarrow'
            frameidx = max([prev_tch,1]); change_frame = 1;
            handles.frameidx = frameidx;
        case 'rightarrow'
            frameidx = min([next_tch,handles.nframes]); change_frame = 1;
            handles.frameidx = frameidx;
        case 'uparrow'
            toggle = 1;
            handles.frameidx = frameidx;
        case 'downarrow'
            toggle = 1;
            handles.frameidx = frameidx;
        case 'space'
            toggle = 1;
            handles.frameidx = frameidx;
        case 's' % Save results
            handles = save_everything(handles);
        case 'd' % Derive from data
            handles = data_driven_function(handles);
    end
    
    % TO DO: Add keyboard shortcut for protraction/retraction
end

% TO DO: ADD MANUAL PROTRACTION/RETRACTION TOGGLE

if toggle
    handles = toggle_touch(handles);
end

if pro_ret_toggle
    handles = toggle_pro_ret(handles);
end


if change_frame
    handles = update_frame(handles);
    handles = plot_touches(handles);
end

guidata(hObject,handles);


% DONE: Add function that counts tracked files of each type and displays
% them. Can also update the html strings after each save.

% Function to SAVE EVERYTHING
function handles = save_everything(handles)

% TO DO: SAVE THESE INTO SEPARATE PARAMS TABLE INSTEAD

touchcsv = 'touch_params.csv';
touch_params = csvread(touchcsv);

% Update touch_params
touch_params(handles.fid,4) = 1; % Touches found
touch_params(handles.fid,5) = 1; % Pole distance
touch_params(handles.fid,6) = handles.dist_thresh;
touch_params(handles.fid,7) = handles.angle_thresh;
touch_params(handles.fid,8) = handles.kappa_thresh;
touch_params(handles.fid,9) = handles.dist_checkbox_val;
touch_params(handles.fid,10) = handles.angle_checkbox_val;
touch_params(handles.fid,11) = handles.kappa_checkbox_val;
touch_params(handles.fid,12) = handles.trialtype_checkbox_val;
touch_params(handles.fid,14) = handles.choice; 
% Adding 'tracked by machine' column to touch_params
touch_params(handles.fid,13) = handles.autotracked;
% Adding 'protraction 1st touch' column
touch_params(handles.fid,15) = handles.protraction_touch;

disp('Writing params file to disk');
dlmwrite('touch_params.csv',touch_params,'precision',6)

touches = handles.touches;

disp(['Saving file ',handles.fname(1:end-4),'_touch.mat']);
if exist('touches','var') % This check might be redundant now we are defining touches at startup
    save([handles.fname(1:end-4),'_touch.mat'],'touches','-append');
else
    error('No touches to save! Find touches before trying to save')
end

% TO DO: Add saving of pole-up-times and protraction/retraction 
protraction_touch = handles.protraction_touch;
pro_ret = handles.pro_ret;
first_touch = handles.first_touch;
pole_up_times = handles.pole_up_times;
pole_updown = handles.pole_updown;

save([handles.fname(1:end-4),'_touch.mat'],'protraction_touch','first_touch','pro_ret','pole_up_times','pole_updown','-append');

set(handles.text1,'BackgroundColor',[0.3,1,0.3])
set(handles.text1,'String','Touches saved')

set(handles.threshold_text,'BackgroundColor',[0.3,1,0.3])
set(handles.threshold_text,'String','THRESHOLD SAVED')

% Update popupmenu1 string list to include bold type
if handles.trialtype == 1
    handles.html_string_list{handles.fid+1} = ['<HTML><FONT COLOR="lime"><b>',handles.string_list{handles.fid+1},'</b></HTML>'];
elseif handles.trialtype == 2
    handles.html_string_list{handles.fid+1} = ['<HTML><FONT COLOR="red"><b>',handles.string_list{handles.fid+1},'</b></HTML>'];
elseif handles.trialtype == 3
    handles.html_string_list{handles.fid+1} = ['<HTML><FONT COLOR="black"><b>',handles.string_list{handles.fid+1},'</b></HTML>'];
else
    disp('Unknown trial type')
    handles.html_string_list{handles.fid+1} = ['<HTML><FONT COLOR="white"><b>',handles.string_list{handles.fid+1},'</b></HTML>'];
end

% Update tracked file counter in the GUI
tracked = find(touch_params(:,4));
handles.ttt = hist(touch_params(tracked,3),linspace(0.5,3.5,3));
handles = update_ttt(handles);

set(handles.popupmenu1,'String',handles.html_string_list)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Executes on button press in batch_tracking.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function batch_tracking_Callback(hObject, eventdata, handles)
% hObject    handle to batch_tracking (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% GOAL: Get the user to manually curate a few files, then the params for
% these files are applied to the rest of the data in the folder

% On push button the user will get a dialog box:
% (a) Apply current params to whole session? (b) Default params
% (c) Derive new params from tracked data

% If (c), we check how many files of each type have been tracked. This
% shoud be noted in the gui, and used as a guide about chances batch
% processing will work.

% After (c), or if the user chose (a) or (b), the code then just finds
% touches in all files based on specified params

touchcsv = 'touch_params.csv';
touch_params = csvread(touchcsv);

% DETERMINE NUMBER OF/IDENTITY OF FILES TO TRACK
% MAYBE JUST IGNORE .AVI AS FUNCTION FOR LOADING IS SO BRITTLE
ff = [dir('*_touch.mat')]; %dir('*.dat'); %

[all_files,file_id,~] = unique(touch_params(:,2));
tracked = find(touch_params(file_id,4));
autotracked = find(touch_params(file_id,13));
tracked(ismember(tracked,autotracked)) = []; % Re-track autotracked files
file_id(tracked) = [];
untracked = file_id;

% keyboard

% Determine correspondence between untracked and ff (for folders with 2
% videos per filename

% this_video = [];
% for i = 1:numel(untracked)
%     for j = 1:numel(ff)
%         [str2double(ff(j).name(end-15:end-10)),touch_params(untracked(i),2)]
%         if str2double(ff(j).name(end-15:end-10)) == touch_params(untracked(i),2) % end-9:end-4 for .dat or .avi
%             this_video(i) = j;
% %             j = numel(ff);
%         end
%     end
% end

DD = 0;

choice = questdlg('Which thesholds to use on the rest of the session?',...
    'Batch touch detection','Current','Default','Derive from tracked files','Default');
switch choice
    case 'Current' % Just use current thresholds/criteria
        set(handles.threshold_text,'BackgroundColor',[0.3,0.3,0.3])
        set(handles.threshold_text,'String','CURRENT THRESHOLD')
        
    case 'Default' % Set thresholds/criteria to defaults
        handles = reset_everything(handles);
        set(handles.threshold_text,'BackgroundColor',[0.3,0.3,0.3])
        set(handles.threshold_text,'String','DEFAULT THRESHOLD')
        
    case 'Derive from tracked files' % Derive thresholds/criteria from tracked data
        DD = 1;
        
        set(handles.dist_checkbox,'value',round(median(touch_params(tracked,9))));
        
        set(handles.angle_checkbox,'value',round(median(touch_params(tracked,10))));
        set(handles.kappa_checkbox,'value',round(median(touch_params(tracked,11))));
        set(handles.trialtype_checkbox,'value',round(median(touch_params(tracked,9))));
        
        % NEED TO WORK OUT HOW TO SET THIS IN A TRIAL TYPE SPECIFIC WAY
        % V1 is below. This is placeholder texts
        handles.dist_thresh = max(touch_params(tracked,6));
        set(handles.dist_text,'String',num2str(handles.dist_thresh));
        set(handles.dist_slider,'value',0);
        
        % NEXT TWO VALUES ARE PLACEHOLDERS
        handles.angle_thresh = max(touch_params(tracked,7));
%         set(handles.angle_text,'String',num2str(handles.angle_thresh));
%         set(handles.angle_slider,'value',0);
        
        handles.kappa_thresh = max(touch_params(tracked,8));
%         set(handles.kappa_text,'String',num2str(handles.kappa_thresh));
%         set(handles.kappa_slider,'value',0);
        
        set(handles.threshold_text,'BackgroundColor',[0.3,0.3,0.3])
        set(handles.threshold_text,'String','BATCH DATA THRESHOLD')
        
        
end

set(handles.text1,'BackgroundColor',[1,1,0]);
set(handles.text1,'String','DETECTING CONTACTS IN BATCH MODE');

% LOOP TO FIND CONTACTS IN THE REST OF THE SESSION
fname_backbone = ff(1).name(1:end-16);

for i = 1:numel(untracked)
    % LOAD DATA
    handles.fid = untracked(i);
    handles.fname = [fname_backbone,num2str(touch_params(untracked(i),2)),'.dat'];
    handles.trialtype = touch_params(untracked(i),3);
    fn = [handles.fname(1:end-4)];
    if exist([fn,'.tr4'],'file')
        load([fn,'.tr4'],'whisker','-mat');   
        handles.r_base = whisker.r3all(:,1:2,:);
        handles.kappa = whisker.kappa_all;
        handles.theta = whisker.theta_all;
        handles.fp = whisker.fp3_all(:,1:2);
        
    elseif exist([fn,'_clean.mat'],'file')
        load([fn,'_clean.mat'],'r_base','theta_w','kappa_w');
        handles.r_base = r_base;
        handles.kappa = kappa_w;
        handles.theta = theta_w;
        handles.fp = squeeze(r_base(:,:,1));
    else
        error(['No .tr or _clean.mat file available. Track video with chimera.m or Whikerman before continuing'])
    end
    
    
    load([fn,'_touch.mat'],'start_frame','trigger_frame','pole_dist','gof','radius','barPos','closest_w');
    handles.start_frame = start_frame;
    handles.trigger_frame = trigger_frame;
    handles.pole_dist = pole_dist;
    handles.gof = gof;
    handles.radius = radius;
    handles.barPos = barPos;
    handles.closest_w = closest_w;
    
    % Check if pole_up_times has been computed
    
    % Else, work out pole_up_times + protraction/retraction
    handles = pole_uptimes(handles);
    
    
    if DD % If using data-driven thesholds, specify threshold on a trial by trial (more accurately a trialtype by trialtype) way
        tt = touch_params(untracked(i),3);
        same_tt = find(touch_params(tracked,3) == tt);
        handles.dist_thresh = max(touch_params(same_tt,6));
    end
    % DETECT CONTACTS
    handles = find_touches(handles);
    
    % DETECT WHETHER FIRST TOUCH IS PROTRACTION/RETRACTION
    handles = protraction_detector(handles);
    
    set(handles.text1,'BackgroundColor',[0,1,0])
    set(handles.text1,'String','Touches found')
    
    % SAVE OUTPUTS and UPDATE touch_params
    handles.autotracked = 1; % Note that this video isn't being tracked by a human
  
    handles = save_everything(handles);
    
    % PLOT CONTACTS
    figure(2);
    clf;
    plot(1:numel(handles.pole_dist),handles.pole_dist,'color',[0.2,0.2,0.2]);
    hold on
    plot(1:numel(handles.gof),handles.gof,'color',[0,0.4470,0.7410]); % 2015 default blue
    plot(handles.pole_up_times,handles.gof(handles.pole_up_times),'r.');
    plot(find(handles.touches),handles.pole_dist(find(handles.touches)),'k.');
    title([handles.fname,' TT: ',num2str(handles.trialtype)])
    %legend({'Pole distance';'goodness of fit';'pole up';'touches'});
    % TO DO. Add switches for different checkboxes
    
    drawnow;
    %     pause

end


% Finally, plot kappa/theta space, dots coloured by trial type, then either darker or circled to indicate touch
handles = pollock_plot(handles);



% DELETE DATA
guidata(hObject,handles);

function handles = update_ttt(handles);
% Display number tracked of each trial type

set(handles.tracked_text1,'BackgroundColor',[0.5,0.5,0.5])
set(handles.tracked_text1,'ForegroundColor',[0,1,0])
set(handles.tracked_text1,'String',[num2str(handles.ttt(1)),'/',num2str(handles.ttu(1))])

set(handles.tracked_text2,'BackgroundColor',[0.5,0.5,0.5])
set(handles.tracked_text2,'ForegroundColor',[1,0,0])
set(handles.tracked_text2,'String',[num2str(handles.ttt(2)),'/',num2str(handles.ttu(2))])

set(handles.tracked_text3,'BackgroundColor',[0.5,0.5,0.5])
set(handles.tracked_text3,'ForegroundColor',[0,0,0])
set(handles.tracked_text3,'String',[num2str(handles.ttt(3)),'/',num2str(handles.ttu(3))])



% --- Executes on button press in session_summary.
function session_summary_Callback(hObject, eventdata, handles)
% hObject    handle to session_summary (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = pollock_plot(handles);
guidata(hObject,handles);


% Generate 'Jackson Pollock' theta/kappa plot coloured by contact
% TO DO: colour points by trial type
function handles = pollock_plot(handles);

% RE-LOAD ALL DATA
% TO DO: use touch_params for this, and do separate plots/colours for
% autotracked files.
ff = dir('*_touch.mat');
kappa_touch = [];
theta_touch = [];
kappa_nontouch = [];
theta_nontouch = [];

for i = 1:numel(ff);
    disp(['Loading data from trial ',ff(i).name(1:end-10)])
    if exist([ff(i).name(1:end-10),'.tr4'],'file')
        load([ff(i).name(1:end-10),'.tr4'],'whisker','-mat');
        kappa = whisker.kappa_all;
        theta = whisker.theta_all;
    elseif exist([ff(i).name(1:end-10),'_clean.mat'],'file')
        load([ff(i).name(1:end-10),'_clean.mat'],'theta_w','kappa_w');
        kappa = kappa_w;
        theta = theta_w;
    else
        error(['No .tr or _clean.mat file available. Track video with chimera.m or Whikerman before continuing'])
    end
    
    load([ff(i).name],'touches')
    
    % CREATE LONG VECTORS OF ALL DATA
    %         plot(touches); title(ff(i).name);drawnow; % debugging
    %     keyboard
    if exist('touches','var')
 %       if numel(find(touches))
            kappa_touch = [kappa_touch,kappa(find(touches))];
            kappa(find(touches)) = [];
            kappa_nontouch = [kappa_nontouch,kappa];
            
            theta_touch = [theta_touch,theta(find(touches))];
            theta(find(touches)) = [];
            theta_nontouch = [theta_nontouch,theta];
%         end
    end
    clear touches
end

% PLOT
figure;
plot(conv(theta_nontouch,gausswin(10,1),'same')./10,conv(kappa_nontouch,gausswin(10,1),'same')./10,'k.','markersize',1)
hold all
plot(conv(theta_touch,gausswin(10,1),'same')./10,conv(kappa_touch,gausswin(10,1),'same')./10,'r.','markersize',1)
set(gca,'Xdir','reverse')
xlabel('Whisker angle')
ylabel('Whisker curvature')

% --- Executes on button press in pole_up_fix.
function pole_up_fix_Callback(hObject, eventdata, handles)
% hObject    handle to pole_up_fix (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

prompt = {'Enter Pole-Up time:','Enter Pole-down time:'};
dlg_title = 'Manual pole-up time';
num_lines = 1;
defaultans = {num2str(handles.pole_updown(1)),num2str(handles.pole_updown(2))};
answer = inputdlg(prompt,dlg_title,num_lines,defaultans);

pole_updown = [str2num(answer{1}),str2num(answer{2})]

% Work out which order the peaks are (depends on timing of trigger in
% the video
if numel(find(pole_updown>=handles.trigger_frame)) == 2||0;
    pole_up_times = pole_updown(1):pole_updown(2);
else
    pole_up_times = [pole_updown(2):numel(handles.pole_dist),1:pole_updown(1)];
end

handles.pole_updown = pole_updown;
handles.pole_up_times = pole_up_times;
handles = plot_touches(handles);

guidata(hObject,handles);

% Function to calculate whether a given moment is protraction is retraction
function handles = pro_ret_calculator(handles);
theta_ts = timeseries(handles.theta,(1:numel(handles.theta))./1000);
bandpass = [6,30];
theta_filt = idealfilter(theta_ts,bandpass,'pass');
H = hilbert(theta_filt.data);

pro = find(angle(H)<=0);

handles.pro_ret = zeros(size(handles.theta));
handles.pro_ret(pro) = 1;


% Function to determine if first touch is protraction or retraction
function handles = protraction_detector(handles)
touches = circshift(handles.touches',[-handles.start_frame,0]);
pro_ret = circshift(handles.pro_ret',[-handles.start_frame,0]);
pro = find(pro_ret);

% % Debug plotting
% figure(5);clf;
% subplot(2,1,1);
% theta = circshift(handles.theta',[-handles.start_frame,0]);
% plot(theta);
% hold all
% plot(pro,theta(pro),'g.')
% plot(find(touches),theta(find(touches)),'k.')
% first_touch = find(touches,1,'first');
% plot(first_touch,theta(first_touch),'m.')
% 
% dummy_array = zeros(size(theta));
% dummy_array(first_touch) = 1;
% first_touch_array = circshift(dummy_array,[handles.start_frame,0]);
% handles.first_touch = find(first_touch_array);
% 
% subplot(2,1,2);
% plot(handles.theta);
% hold all
% plot(find(handles.pro_ret),handles.theta(find(handles.pro_ret)),'g.')
% plot(find(handles.touches),handles.theta(find(handles.touches)),'k.')
% plot(handles.first_touch,handles.theta(handles.first_touch),'m.')

handles.first_touch = 0;
handles.protraction_touch = 0;

if numel(find(touches)) >= 1
    first_touch = find(touches,1,'first');
    handles.first_touch = mod(first_touch + handles.start_frame, numel(handles.gof));
    
    
    if ismember(first_touch,pro)
        handles.protraction_touch = 1;
        
        % Compute first touch frame in un-circshifted video
        dummy_array = zeros(size(pro_ret));
        dummy_array(first_touch) = 1;
        first_touch_array = circshift(dummy_array,[handles.start_frame,0]);
        handles.first_touch = find(first_touch_array);
        
    else
        handles.protraction_touch = 2;
    end
end

handles.first_touch

% keyboard
handles = pro_ret_text_set(handles);

    
    
    
% --- Executes on button press in pro_ret_button.
function pro_ret_button_Callback(hObject, eventdata, handles)
% hObject    handle to pro_ret_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pro_ret_button
handles = pro_ret_toggle(handles);


guidata(hObject,handles);

% Update Protraction/Retraction text colour/string
function handles = pro_ret_text_set(handles)
% Protraction
handles.protraction_touch
if handles.protraction_touch == 1
    set(handles.pro_ret_text,'BackgroundColor',[0.5,0.5,0.5])
    set(handles.pro_ret_text,'ForegroundColor','g')
    set(handles.pro_ret_text,'String','Protraction')
elseif handles.protraction_touch == 2
    set(handles.pro_ret_text,'BackgroundColor',[0.5,0.5,0.5])
    set(handles.pro_ret_text,'ForegroundColor','m')
    set(handles.pro_ret_text,'String','Retraction')
else
    set(handles.pro_ret_text,'BackgroundColor',[0.5,0.5,0.5])
    set(handles.pro_ret_text,'ForegroundColor','k')
    set(handles.pro_ret_text,'String','No Touches')
end

% Toggle protraction/retraction first touch
function handles = pro_ret_toggle(handles)

if handles.protraction_touch == 1
    handles.protraction_touch = 2;
elseif handles.protraction_touch == 2;
    handles.protraction_touch = 1;
end
handles = pro_ret_text_set(handles);

% --- Open video file
function handles = open_video(handles);
switch handles.fname(end-2:end)
    
    case 'avi'
        video.type = 'avi';
        video.vObj = VideoReader(handles.fname);
        handles.nframes = video.vObj.NumberOfFrames;
        video.width = video.vObj.Width;
        video.height = video.vObj.Height;
    case 'dat'
        video.type = 'dat';
%         fclose('all');
        video.fid = fopen(handles.fname,'r')
        video.header = read_mikrotron_datfile_header(video.fid);
        handles.nframes = video.header.nframes;
        video.width = video.header.width;
        video.height = video.header.height;
        video.offset = 8192;
        % set file position to start of first frame
        fseek(video.fid,8192,-1);
        
    otherwise
        fprintf('Please choose a video file to display instead of an analysis output file \n \n')
        msg = ['Unhandled file type: ',handles.fname(end-3:end)];
        error(msg)
        
end
handles.video = video;
    

        
