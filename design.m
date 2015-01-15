%Light Study Data Automation GUI
%Neesirg Patel, nmp52(at)pitt.edu
%Yates Lab, University of Pittsburgh
%Last Updated: 1/14/15 @ 11:00AM by Neesirg Patel;

function varargout = design(varargin)
% DESIGN M-file for design.fig
%      DESIGN, by itself, creates a new DESIGN or raises the existing
%      singleton*.
%
%      H = DESIGN returns the handle to a new DESIGN or the handle to
%      the existing singleton*.
%
%      DESIGN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DESIGN.M with the given input arguments.
%
%      DESIGN('Property','Value',...) creates a new DESIGN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before design_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to design_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help design

% Last Modified by GUIDE v2.5 14-Jan-2015 11:02:00

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @design_OpeningFcn, ...
                   'gui_OutputFcn',  @design_OutputFcn, ...
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


% --- Executes just before design is made visible.
function design_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for design
handles.output = hObject;
movegui('center');
% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = design_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Data Input Selection
function pushbutton1_Callback(hObject, eventdata, handles)
handles.in = uigetdir('C:\Users\', 'Select the data input path');
if isequal(handles.in, 0)
   guidata(hObject, handles);
else
   set(findobj('Tag','text5'),'String',['Input Path: ' handles.in]);
   set(handles.pushbutton2, 'Visible', 'on');
   guidata(hObject, handles);
end


% --- Date Selection
function pushbutton2_Callback(hObject, eventdata, handles)
uicalendar('Weekend', [1 0 0 0 0 0 1], 'DestinationUI', {handles.pushbutton2, {'String'}});
waitfor(handles.pushbutton2, 'String');
if or(strcmp(get(handles.pushbutton2, 'String'),'Select Date'), strcmp(get(handles.pushbutton2, 'String'),''))
    set(handles.pushbutton2, 'String', 'Select Date');
    guidata(hObject, handles);
else
    handles.date = get(handles.pushbutton2, 'String');
    set(handles.pushbutton3, 'Enable', 'on');
    guidata(hObject, handles);
end


% --- Data Output Selection
function pushbutton3_Callback(hObject, eventdata, handles)
[handles.outFileName, handles.outPathName] = uigetfile({'*.xls;*.xlsx','All Excel files';}, 'Select the data output file');
if isequal(handles.outFileName,0)
   guidata(hObject, handles);
else
   handles.out = [handles.outPathName, handles.outFileName]; 
   set(findobj('Tag','text4'),'String',['Output Path: ' fullfile(handles.outPathName, handles.outFileName)]);
   set(handles.pushbutton4, 'Visible', 'on');
   guidata(hObject, handles);
end


% --- Callback for 'Analyze' Button
function pushbutton4_Callback(hObject, eventdata, handles)
dataMatrix = analyzeRVLM(handles.in); % call automation function for RVLM study
a = xlsread(handles.out);
nRows = size(a,1);
set(handles.pushbutton4, 'Enable', 'off');
set(handles.pushbutton5, 'Visible', 'on');
if (size(dataMatrix) == [0,0])
    msgbox('No data was written to the output file', 'Input Error');
    guidata(hObject, handles);
else
    xlswrite(handles.out, cellstr(handles.date), 'Sheet1', strcat('A', num2str(nRows+2), ':', 'A', num2str(nRows+1+size(dataMatrix,1))));
    xlswrite(handles.out, dataMatrix, 'Sheet1', strcat('B', num2str(nRows+2)));
    guidata(hObject, handles);
end


% --- Callback for 'Start Over' Button
function pushbutton5_Callback(hObject, eventdata, handles)
close(gcbf);
design;
