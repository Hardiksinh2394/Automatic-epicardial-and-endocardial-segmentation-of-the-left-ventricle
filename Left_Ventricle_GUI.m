function varargout = Left_Ventricle_GUI(varargin)
% LEFT_VENTRICLE_GUI MATLAB code for Left_Ventricle_GUI.fig
%      LEFT_VENTRICLE_GUI, by itself, creates a new LEFT_VENTRICLE_GUI or raises the existing
%      singleton*.
%
%      H = LEFT_VENTRICLE_GUI returns the handle to a new LEFT_VENTRICLE_GUI or the handle to
%      the existing singleton*.
%
%      LEFT_VENTRICLE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LEFT_VENTRICLE_GUI.M with the given input arguments.
%
%      LEFT_VENTRICLE_GUI('Property','Value',...) creates a new LEFT_VENTRICLE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Left_Ventricle_GUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Left_Ventricle_GUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Left_Ventricle_GUI

% Last Modified by GUIDE v2.5 23-May-2019 12:27:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Left_Ventricle_GUI_OpeningFcn, ...
                   'gui_OutputFcn',  @Left_Ventricle_GUI_OutputFcn, ...
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


% --- Executes just before Left_Ventricle_GUI is made visible.
function Left_Ventricle_GUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Left_Ventricle_GUI (see VARARGIN)

% Choose default command line output for Left_Ventricle_GUI
handles.output = hObject;
index = 0;
handles.index = index;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Left_Ventricle_GUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Left_Ventricle_GUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in LOAD_NIFTI.
function LOAD_NIFTI_Callback(hObject, eventdata, handles)
% hObject    handle to LOAD_NIFTI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename pathname] = uigetfile('*.nii','select nifti file');
completename = strcat(pathname,filename);

Img_4d = niftiread(completename);

[A, B, C, D] = size(Img_4d);

[FTimg, DC_AVG, Harmonic_1] = HeartROI_2(Img_4d, A, B, C, D);


axes(handles.Slice_AVG);
imshow(DC_AVG(:,:,1),[]);
handles.size_3 = C;
handles.average = DC_AVG;
handles.H1 = Harmonic_1;

Mask_ROI = zeros(size(DC_AVG));
Epi = zeros(size(DC_AVG));
Endo = zeros(size(DC_AVG));

handles.EPI = Epi;
handles.Endo = Endo;
handles.ROI_disp = Mask_ROI;
handles.filename = filename;
guidata(hObject,handles)


% --- Executes on slider movement.
function Switch_Callback(hObject, eventdata, handles)
% hObject    handle to Switch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
DC_AVG = handles.average;

Epi = handles.EPI;
Endo = handles.Endo;
ROI = handles.ROI_disp;

minimum = handles.index;
maximum = handles.size_3;
set(handles.Switch,'Min',minimum + 1 );
set(handles.Switch,'Max',maximum);
set(handles.Switch, 'SliderStep', [1/(maximum-1) , 10/(maximum-1) ]);
A = get(hObject,'Value');
A = round(A);
axes(handles.Slice_AVG);
% https://fr.mathworks.com/matlabcentral/answers/405856-i-have-a-gui-slider-that-calls-a-value-from-an-array-that-gets-the-error-index-in-position-2-is-inv
% useful link
tmp = DC_AVG(:,:,A);
axes(handles.Slice_AVG);
imshow(tmp,[]);

handles.A = A;

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
axes(handles.Segmented_display);
imshow(tmp,[]);
hold on;
visboundaries(Endo(:,:,A), 'LineWidth', 0.3, 'Color', 'g');
hold on;
visboundaries(Epi(:,:,A), 'LineWidth', 0.3, 'Color', 'r');
hold on;
visboundaries(ROI(:,:,A), 'LineWidth', 0.6, 'Color', 'y');

I = getframe(handles.Segmented_display);
handles.I = I;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function Switch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Switch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton20.
function pushbutton20_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in LV_Segmentation.
function LV_Segmentation_Callback(hObject, eventdata, handles)
% hObject    handle to LV_Segmentation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Harmonic_1 = handles.H1;
DC_AVG = handles.average;
C = handles.size_3;
A = handles.A;
[ROI_slices, H1_ROI, H1_post, Mask_ROI] = Heart_detector_3(Harmonic_1, DC_AVG);
%[ROI_slices_2, H1_ROI_2, H1_post_2] = Heart_detector_part2(H1_ROI, ROI_slices);

[refined_mask,Extracted_ROI,H1_Extracted] = ROI_extractor(ROI_slices,H1_ROI, Mask_ROI);


Pol_ROI = Centroid_pol(Extracted_ROI, DC_AVG);
%Pol_ROI = zeros(size(DC_AVG));
% Do polar transform wrt weighted centroid from H1_Extracted

ROI_padd = Niveleur(Extracted_ROI, Pol_ROI);
ROI_padd2 = Nivelage(ROI_slices, Pol_ROI, H1_ROI);

BW = zeros(size(ROI_padd));
BW2 = BW;

for i = 1:C
    [counts,x] = imhist(uint8(ROI_padd(:,:,i)), 16);
    stem(x,counts)
    T1 = otsuthresh(counts);
    BW(:,:,i) = imbinarize(uint8(ROI_padd(:,:,i)),T1);
end


for i = 1:C
    [counts,x] = imhist(histeq(uint8(ROI_padd2(:,:,i)), 16));
    stem(x,counts)
    T2 = otsuthresh(counts);
    BW2(:,:,i) = imbinarize(uint8(ROI_padd2(:,:,i)), T2);
    if BW2(:,:,i) == 0
        return
    end
end


Endo = Endocardium(BW);
Epi = Epicardium(BW2);


handles.ROI = Mask_ROI;
handles.EPI = Epi;
handles.Endo = Endo;
guidata(hObject,handles)







% --- Executes on button press in ROI_Button.
function ROI_Button_Callback(hObject, eventdata, handles)
% hObject    handle to ROI_Button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Mask_ROI = handles.ROI;

handles.ROI_disp = Mask_ROI;
guidata(hObject,handles)


% --- Executes on button press in JPEG_Converter.
function JPEG_Converter_Callback(hObject, eventdata, handles)
% hObject    handle to JPEG_Converter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
size3 = handles.size_3;
filename = handles.filename;
Slice_num = handles.A;
Slice_chosen = handles.I;

folder = strcat('JPEG_images', filename);

if isfolder(folder)
   rmdir(folder,'s');
   mkdir(folder);
else
   mkdir(folder);
end        
ite = Slice_num; 

    if ite < 10
        s = sprintf(strcat(folder,'/segmented_slice000%d'), ite);
    else
        s = sprintf(strcat(folder,'/segmented_slice00%d'), ite);
    end
            
    capture = frame2im(Slice_chosen);
    imwrite(capture, strcat(s,'.jpg'));

guidata(hObject,handles)


% --- Executes on button press in Clear_figure.
function Clear_figure_Callback(hObject, eventdata, handles)
% hObject    handle to Clear_figure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cla(handles.Segmented_display);
cla(handles.Slice_AVG);

