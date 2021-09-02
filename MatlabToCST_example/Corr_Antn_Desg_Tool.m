function varargout = Corr_Antn_Desg_Tool(varargin)
% CORR_ANTN_DESG_TOOL MATLAB code for Corr_Antn_Desg_Tool.fig
%      CORR_ANTN_DESG_TOOL, by itself, creates a new CORR_ANTN_DESG_TOOL or raises the existing
%      singleton*.
%
%      H = CORR_ANTN_DESG_TOOL returns the handle to a new CORR_ANTN_DESG_TOOL or the handle to
%      the existing singleton*.
%
%      CORR_ANTN_DESG_TOOL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CORR_ANTN_DESG_TOOL.M with the given input arguments.
%
%      CORR_ANTN_DESG_TOOL('Property','Value',...) creates a new CORR_ANTN_DESG_TOOL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Corr_Antn_Desg_Tool_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Corr_Antn_Desg_Tool_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Corr_Antn_Desg_Tool

% Last Modified by GUIDE v2.5 24-Feb-2020 00:01:32

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Corr_Antn_Desg_Tool_OpeningFcn, ...
                   'gui_OutputFcn',  @Corr_Antn_Desg_Tool_OutputFcn, ...
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


% --- Executes just before Corr_Antn_Desg_Tool is made visible.
function Corr_Antn_Desg_Tool_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Corr_Antn_Desg_Tool (see VARARGIN)

% Choose default command line output for Corr_Antn_Desg_Tool
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Corr_Antn_Desg_Tool wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Corr_Antn_Desg_Tool_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles) %edit1 es la fmin
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles) %edit 2 es la f_max
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit7_Callback(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit7 as text
%        str2double(get(hObject,'String')) returns contents of edit7 as a double


% --- Executes during object creation, after setting all properties.
function edit7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);
fmin = str2double(get(handles.edit1, 'String')); % in ghz
handles.fmin=fmin;
fmax = str2double(get(handles.edit2,'string')); % in ghz
handles.fmax=fmax;
scl=fmax/fmin;

if scl < 1.4
    fc=sqrt(fmin*fmax);
    fo=1.02*fc;
else
    fc=1.2*fmin;
    fo=1.1*fc;
end

lam_c=physconst('LightSpeed')/(fc*10^9);
lam_c_mm = lam_c * 1000;
lam_o=physconst('LightSpeed')/(fo*10^9);
lam_o_mm = lam_o * 1000;

ai=3*lam_c_mm/(2*pi); %in mm
%from fig 3
ao=str2double(get(handles.edit3,'string'))*lam_c_mm; %in mm
%pitch between 5 (narrowband) - 10 (broadband)
p=lam_c_mm/str2double(get(handles.edit4,'string'));
% pitch to width ratio (affect XPD)
delta = str2double(get(handles.edit5,'string'));
sigma = str2double(get(handles.edit6,'string')); % sigma 0.4 - 0.5
kc = (2*pi)/lam_c_mm;            % Wave number at center frequency
ko = (2*pi)/lam_o_mm;            % Wave number at output frequency
wgl = str2double(get(handles.edit7,'string')); % input waveguide length

%Number of total corrugations
N = str2double(get(handles.edit8,'string'));
len = N*p;
%Number of mode converters
n_mc = str2double(get(handles.edit9,'string'));
z=0:p:len;

% delta2 is only used for case 2, ring loaded slot mode converter
delta2 = str2double(get(handles.edit10,'string'));

% delta_min is only used for case 3,variable pitch to width slot mode converter
delta_min = str2double(get(handles.edit11,'string')); % Should be greater than 0.125 and less than delta

%1=LINEAR, 2=SINUSOID, 3=ASYMMETRIC SINE-SQUARED, 4=TANGENTIAL,
% 5=x.rho, 6=EXPONENTIAL, 7=HYPERBOLIC, 8=POLYNOMIAL

contents = get(handles.popupmenu1,'String');
popupmenu1value = contents{get(handles.popupmenu1,'Value')};

switch popupmenu1value
    case 'Linear'
        horn_profile = 1;
    case 'Sinusoidal'
        horn_profile = 2;
    case 'Asymmetric Sine-Squared'
        horn_profile = 3;
    case 'Tangential'
        horn_profile = 4;
    case 'x.rho'
        horn_profile = 5;
    case 'Exponential'
        horn_profile = 6;
    case 'Hyperbolic'
        horn_profile = 7;
    case 'Poynomial'
        horn_profile = 8;
end

%1=VARIABLE SLOT DEPTH MC,
% 2=RING LOADED SLOTS MC, 3=VARIABLE PITCH TO WIDTH SLOT MC
contents = get(handles.popupmenu2,'String');
popupmenu2value = contents{get(handles.popupmenu2,'Value')};
switch popupmenu2value
    case 'Variable Slot Depth'
        mode_converter_type = 1;
    case 'Ring Loaded Slots'
        mode_converter_type = 2;
    case 'Variable p/w Slot'
        mode_converter_type = 3;
end

switch(horn_profile)
    case 1
        %%% Linear profile %%%
        a = ai+(ao-ai)*z/len;
        
        plot(z, a);
        %set(gca, "linewidth",2, "fontsize", 14 )
        axis equal;
        xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
        ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
        title( 'Linear Horn Profile', 'FontSize', 16 );
        
    case 2
        %%% Sinusoid profile %%%
        A = 1;     % Amplitude factor 'A' should be between 0 and 1
        rho = 2;     % rho should be between 0.5 and 5, default is 2
        a = ai+(ao-ai)*((1-A)*(z/len)+A*power(sin((pi*z)/(2*len)),rho));
        
        plot(z, a);
        %set(gca, "linewidth",2, "fontsize", 14 )
        axis equal;
        xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
        ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
        title( 'Sinusoidal Horn Profile', 'FontSize', 16 );
        
    case 3
        %%% Asymmetric Sine-Squared profile %%%
        L1 = len/3;    % Choose a value for L1, must be less than the horn length
        L2 = len-L1;   % L2 is the length between L1 and the end of the horn
        gamma = L2/L1;
        idx = find(z <= L1);     % Find the index of z corresponding to L1
        zelements = size(z,2);   % Total number of points in z axis of the horn
        
        za = z(1: max(idx));
        aa = ai+((2*(ao-ai))/(1+gamma))*sin((pi*za)/(4*L1)).^2;
        
        zb = z(max(idx)+1 : zelements);
        ab = ai+((2*(ao-ai))/(1+gamma))*(gamma*sin(((pi*(zb+L2-L1))/(4*L2))).^2+((1-gamma)/2));
        
        a = [aa,ab];
        z = [za,zb];
        
        plot(z, a);
        %set(gca, "linewidth",2, "fontsize", 14 )
        axis equal;
        xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
        ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
        title( 'Asymmetric Sine Squared Horn Profile', 'FontSize', 16 );
        
    case 4
        %%% Tangential profile %%%
        A = 1;
        rho = 2;
        a = ai+(ao-ai)*((1-A)*(z/len)+A*power(tan((pi*z)/(4*len)),rho));
        
        plot(z, a);
        %set(gca, "linewidth",2, "fontsize", 14 )
        axis equal;
        xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
        ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
        title( 'Tangential Horn Profile', 'FontSize', 16 );
        
    case 5
        %%% x.rho profile %%%
        A = 1;
        rho = 2;
        a = ai+(ao-ai)*((1-A)*(z/len)+A*power(z/len,rho));
        
        plot(z, a);
        %%set(gca, "linewidth",2, "fontsize", 14 )
        axis equal;
        xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
        ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
        title( 'xp Horn Profile', 'FontSize', 16 );
        
    case 6
        %%% Exponential profile %%%
        a=ai*exp(log(ao/ai)*(z/len));
        
        plot(z, a);
        %%set(gca, "linewidth",2, "fontsize", 14 )
        axis equal;
        xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
        ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
        title( 'Exponential Horn Profile', 'FontSize', 16 );
        
    case 7
        %%% Hyperbolic profile %%%
        a = sqrt(ai^2 + (power(z,2) * (ao^2-ai^2) / len^2));
        
        plot(z, a);
        %set(gca, "linewidth",2, "fontsize", 14 )
        axis equal;
        xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
        ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
        title( 'Hyperbolic Horn Profile', 'FontSize', 16 );
        
    case 8
        %%% POLYNOMIAL Profile %%%
        rho = 3;
        a=ai+(rho+1)*(ao-ai)*(1-((rho*z)/((rho+1)*len))).*power(z/len,rho);
        
        plot(z, a);
        %set(gca, "linewidth",2, "fontsize", 14 )
        axis equal;
        xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
        ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
        title( 'Polynomial Horn Profile', 'FontSize', 16 );
        
end

handles.a=a;
switch(mode_converter_type)
    
    case 1                            % case 1=VARIABLE SLOT DEPTH MODE CONVERTER
        % Mode Converter depths for element j
        ajmc = a(1:n_mc);                     % Index range for mode converter
        idx = 1:n_mc;
        djmc = (sigma-((idx-1)./n_mc).*(sigma-(0.25.*exp(1./(2.114.*(kc*ajmc).^1.134)))))*lam_c_mm;
        % Depth of remaining corrugations
        aj = a(n_mc+1:end);
        idx = n_mc+1:N+1;
        dj = ((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-n_mc-1)/(N-n_mc-1)).*((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lam_o_mm/4).*exp(1./(2.114.*(ko*ao).^1.134)));
        d = [djmc, dj];       % Combining the mode converter and horn depth values
        
        % Generate z,y coordinates as len and rad vector
        n = 0;
        lent(1) = 0;
        lent(2) = 0;
        for i = 1:N
            rad(i+n) = a(i);
            rad(i+n+1) = a(i)+d(i);
            rad(i+n+2) = a(i)+d(i);
            rad(i+n+3) = a(i+1);
            rad(i+n+4) = a(i+1);
            lent(i+n+2) = lent(i+n)+delta*p;
            lent(i+n+3) = lent(i+n+2);
            lent(i+n+4) = lent(i+n+3)+(1-delta)*p;
            lent(i+n+5) = lent(i+n+4);
            n = n+3;
        end
        
        z_number = (N*4)+1; % Number of coordinate points for corrugated length of horn
        lent = lent(1:z_number); % Truncate z axis data points to equal rad vector length
        
%         figure
%         plot(lent,rad);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Mode Converter and Corrugation Profile', 'FontSize', 16 );
%         axis equal;   % Scale axis equally for aspect ratio 1:1
        
    case 2                             % case 2=RING LOADED SLOT MODE CONVERTER
        % Mode Converter depths for element j
        ajmc = a(1:n_mc);                 % Index range for mode converter
        idx = 1:n_mc;
        djmc = (lam_c_mm/4).*exp(1./(2.114.*(kc*ajmc).^1.134));
        % Width of bjth slot for mode converter
        bj = (0.1+(idx-1).*((delta2-0.1)./n_mc)).*p;
        % Height of hjth slot for mode converter
        hj = (2/3).*djmc;
        
        % Depth of remaining corrugations
        aj = a(n_mc+1:end);
        idx = n_mc+1:N+1;
        dj = ((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-n_mc-1)/(N-n_mc-1)).*((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lam_o_mm/4).*exp(1./(2.114.*(ko*ao).^1.134)));
        d = [djmc, dj];       % Combining the mode converter and horn depth values
        
        % Generate z,y coordinates as lent,rad vector
        n = 5;
        lent = [0, 0, (-delta*p)+bj(1), (-delta*p)+bj(1), bj(1), bj(1)];
        rad = [a(1), a(1)+d(1)-hj(1), a(1)+d(1)-hj(1), a(1)+d(1), a(1)+d(1), a(2)];
        for i = 2:n_mc
            rad(i+n) = a(i);
            rad(i+n+1) = a(i)+d(i)-hj(i);
            rad(i+n+2) = a(i)+d(i)-hj(i);
            rad(i+n+3) = a(i)+d(i);
            rad(i+n+4) = a(i)+d(i);
            rad(i+n+5) = a(i+1);
            lent(i+n) = i*p-p;
            lent(i+n+1) = i*p-p;
            lent(i+n+2) = lent(i+n+1)-(delta*p)+bj(i);
            lent(i+n+3) = lent(i+n+1)-(delta*p)+bj(i);
            lent(i+n+4) = lent(i+n+3)+(delta*p)+bj(i);
            lent(i+n+5) = lent(i+n+3)+(delta*p)+bj(i);
            n = n+5;
        end
        % Add extra coordinate points before remaining corrugations
        lent(n_mc*(n_mc+1)+1) = lent(n_mc*(n_mc+1))+(1-delta)*p;
        rad(n_mc*(n_mc+1)+1) = a(n_mc+1);
        
%         figure
%         plot(lent,rad);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Ring Loaded Mode Converter Profile', 'FontSize', 16 );
%         axis equal;   % Scale axis equally for aspect ratio 1:1
        
        n = n+n_mc+1;
        for i = n_mc+1:N
            rad(n) = a(i);
            rad(n+1) = a(i)+d(i);
            rad(n+2) = a(i)+d(i);
            rad(n+3) = a(i+1);
            rad(n+4) = a(i+1);
            lent(n+1) = lent(n);
            lent(n+2) = lent(n+1)+delta*p;
            lent(n+3) = lent(n+2);
            lent(n+4) = lent(n+3)+(1-delta)*p;
            n = n+4;
        end
        
        z_number = (n_mc*2)+(N*4)+1; % Number of coordinate points for corrugated length of horn
        
%         figure
%         plot(lent,rad);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Ring Loaded Mode Converter Horn Profile', 'FontSize', 16 );
%         axis equal;   % Scale axis equally for aspect ratio 1:1
        
    case 3                 % case 3=VARIABLE PITCH TO WIDTH SLOT MODE CONVERTER
        % Mode Converter depths for element j
        ajmc = a(1:n_mc);                     % First indexes for mode converter
        idx = 1:n_mc;
        djmc = (sigma*(lam_c_mm/1.15)+((idx-1)./(n_mc-1)).*(lam_c_mm/4-(sigma*lam_c_mm/1.15))).*exp(1./(2.114.*(kc*ajmc).^1.134));
        % Depth of remaining corrugations
        aj = a(n_mc+1:end);
        idx = n_mc+1:N+1;
        dj = ((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-n_mc-1)/(N-n_mc-1)).*((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lam_o_mm/4).*exp(1./(2.114.*(ko*ao).^1.134)));
        d = [djmc, dj];       % Combining the mode converter and horn depth values
        
        % Generate z,y coordinates as lent,rad vector
        n = 0;
        lent(1) = 0;
        lent(2) = 0;
        for i = 1:n_mc
            rad(i+n) = a(i);
            rad(i+n+1) = a(i)+d(i);
            rad(i+n+2) = a(i)+d(i);
            rad(i+n+3) = a(i+1);
            rad(i+n+4) = a(i+1);
            lent(i+n+2) = lent(i+n)+(delta_min+(((i-1)./(n_mc-1))*(delta-delta_min)))*p;
            lent(i+n+3) = lent(i+n+2);
            lent(i+n+4) = lent(i+n+3)+(1-(delta_min+(((i-1)./(n_mc-1))*(delta-delta_min))))*p;
            lent(i+n+5) = lent(i+n+4);
            n = n+3;
        end
        
%         figure
%         plot(lent(1:end-1),rad);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Variable Pitch To Width Mode Converter Profile', 'FontSize', 16 );
%         axis equal;   % Scale axis equally for aspect ratio 1:1
        
        for i = n_mc+1:N
            rad(i+n) = a(i);
            rad(i+n+1) = a(i)+d(i);
            rad(i+n+2) = a(i)+d(i);
            rad(i+n+3) = a(i+1);
            rad(i+n+4) = a(i+1);
            lent(i+n+2) = lent(i+n)+delta*p;
            lent(i+n+3) = lent(i+n+2);
            lent(i+n+4) = lent(i+n+3)+(1-delta)*p;
            lent(i+n+5) = lent(i+n+4);
            n = n+3;
        end
        
        z_number = (N*4)+1; % Number of coordinate points for corrugated length of horn
        lent = lent(1:z_number); % Truncate z axis data points to equal rad vector length
        
%         figure
%         plot(lent,rad);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Variable Pitch To Width Mode Converter Horn Profile', 'FontSize', 16 );
%         axis equal;   % Scale axis equally for aspect ratio 1:1
        
end

% Add the rest of the geometry to create a closed path
% a_offset is the inner horn profile shifted up to give the horn a thickness
a_offset = a+(lam_c_mm/2+2);
%figure;                    % Uncomment these three lines for debugging
%plot(z, a_offset);
%axis equal;

% Add vertical surface at horn aperture
lent = [lent, lent(z_number)];
rad = [rad, a_offset(N)];
radmsh=rad;                 % radmesh to fix mesh lines to corrugations
%figure;                    % Uncomment these three lines for debugging
%plot(lent, rad);
%axis equal;

% Flip outer surface profile so that widest horn dimensions comes next in the outline coordinates
outer_surface = fliplr(a_offset);
z_flip = fliplr(z);
extent = lent(end);  % Fudge to make horn aperture planar for ring loaded slot MC
z_flip(1) = extent; % Fudge to make horn aperture planar for ring loaded slot MC
% Add outer profile and circular waveguide to horn
lent = [lent, z_flip, -wgl, -wgl, 0];
rad = [rad, outer_surface, ai+(lam_c_mm/2+2), ai, ai];

axes(handles.axes1);
plot(lent,rad);
xlim([-wgl-10 N*p+10])
%set(gca, "linewidth",2, "fontsize", 14 )
xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
title( 'Complete Corrugated Horn Profile', 'FontSize', 16 );
grid on


axis equal;   % Scale axis equally for aspect ratio 1:1
guidata(hObject, handles);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles = guidata(hObject);

addpath(genpath('C:\Users\seymur\Google Drive\CST\cst-matlab\CST-MATLAB-API-master'));

cst = actxserver('CSTStudio.application');

mws = cst.invoke('NewMWS');

%CstDefaultUnits(mws)
fmin=handles.fmin;
fmax=handles.fmax;
a=handles.a;
%CstDefineFrequencyRange(mws,fmin,fmax)

%CstMeshInitiator(mws)

Xmin='expanded open';
Xmax='expanded open';
Ymin='expanded open';
Ymax='expanded open';
Zmin='expanded open';
Zmax='expanded open';

minfrequency = fmin;
CstDefineOpenBoundary(mws,minfrequency,Xmin,Xmax,Ymin,Ymax,Zmin,Zmax)

XminSpace = 0;
XmaxSpace = 0;
YminSpace = 0;
YmaxSpace = 0;
ZminSpace = 0;
ZmaxSpace = 0;

CstDefineBackroundMaterial(mws,XminSpace,XmaxSpace,YminSpace,YmaxSpace,ZminSpace,ZmaxSpace)

CstCopperAnnealedLossy(mws)


scl=fmax/fmin;

if scl < 1.4
    fc=sqrt(fmin*fmax);
    fo=1.02*fc;
else
    fc=1.2*fmin;
    fo=1.1*fc;
end

lam_c=physconst('LightSpeed')/(fc*10^9);
lam_c_mm = lam_c * 1000;
lam_o=physconst('LightSpeed')/(fo*10^9);
lam_o_mm = lam_o * 1000;

ai=3*lam_c_mm/(2*pi); %in mm
%from fig 3
ao=str2double(get(handles.edit3,'string'))*lam_c_mm; %in mm
%pitch between 5 (narrowband) - 10 (broadband)
p=lam_c_mm/str2double(get(handles.edit4,'string'));
% pitch to width ratio (affect XPD)
delta = str2double(get(handles.edit5,'string'));
sigma = str2double(get(handles.edit6,'string')); % sigma 0.4 - 0.5
kc = (2*pi)/lam_c_mm;            % Wave number at center frequency
ko = (2*pi)/lam_o_mm;            % Wave number at output frequency
wgl = str2double(get(handles.edit7,'string')); % input waveguide length

%Number of total corrugations
N = str2double(get(handles.edit8,'string'));
len = N*p;
%Number of mode converters
n_mc = str2double(get(handles.edit9,'string'));
z=0:p:len;

% delta2 is only used for case 2, ring loaded slot mode converter
delta2 = str2double(get(handles.edit10,'string'));

% delta_min is only used for case 3,variable pitch to width slot mode converter
delta_min = str2double(get(handles.edit11,'string')); % Should be greater than 0.125 and less than delta


%1=LINEAR, 2=SINUSOID, 3=ASYMMETRIC SINE-SQUARED, 4=TANGENTIAL,
% 5=x.rho, 6=EXPONENTIAL, 7=HYPERBOLIC, 8=POLYNOMIAL

contents = get(handles.popupmenu1,'String');
popupmenu1value = contents{get(handles.popupmenu1,'Value')};

switch popupmenu1value
    case 'Linear'
        horn_profile = 1;
    case 'Sinusoidal'
        horn_profile = 2;
    case 'Asymmetric Sine-Squared'
        horn_profile = 3;
    case 'Tangential'
        horn_profile = 4;
    case 'x.rho'
        horn_profile = 5;
    case 'Exponential'
        horn_profile = 6;
    case 'Hyperbolic'
        horn_profile = 7;
    case 'Poynomial'
        horn_profile = 8;
end

%1=VARIABLE SLOT DEPTH MC,
% 2=RING LOADED SLOTS MC, 3=VARIABLE PITCH TO WIDTH SLOT MC
contents = get(handles.popupmenu2,'String');
popupmenu2value = contents{get(handles.popupmenu2,'Value')};
switch popupmenu2value
    case 'Variable Slot Depth'
        mode_converter_type = 1;
    case 'Ring Loaded Slots'
        mode_converter_type = 2;
    case 'Variable p/w Slot'
        mode_converter_type = 3;
end

switch(horn_profile)
    case 1
        %%% Linear profile %%%
        a = ai+(ao-ai)*z/len;
        
%         plot(z, a);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         axis equal;
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Linear Horn Profile', 'FontSize', 16 );
        
    case 2
        %%% Sinusoid profile %%%
        A = 1;     % Amplitude factor 'A' should be between 0 and 1
        rho = 2;     % rho should be between 0.5 and 5, default is 2
        a = ai+(ao-ai)*((1-A)*(z/len)+A*power(sin((pi*z)/(2*len)),rho));
%         
%         plot(z, a);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         axis equal;
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Sinusoidal Horn Profile', 'FontSize', 16 );
%         
    case 3
        %%% Asymmetric Sine-Squared profile %%%
        L1 = len/3;    % Choose a value for L1, must be less than the horn length
        L2 = len-L1;   % L2 is the length between L1 and the end of the horn
        gamma = L2/L1;
        idx = find(z <= L1);     % Find the index of z corresponding to L1
        zelements = size(z,2);   % Total number of points in z axis of the horn
        
        za = z(1: max(idx));
        aa = ai+((2*(ao-ai))/(1+gamma))*sin((pi*za)/(4*L1)).^2;
        
        zb = z(max(idx)+1 : zelements);
        ab = ai+((2*(ao-ai))/(1+gamma))*(gamma*sin(((pi*(zb+L2-L1))/(4*L2))).^2+((1-gamma)/2));
        
        a = [aa,ab];
        z = [za,zb];
        
%         plot(z, a);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         axis equal;
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Asymmetric Sine Squared Horn Profile', 'FontSize', 16 );
        
    case 4
        %%% Tangential profile %%%
        A = 1;
        rho = 2;
        a = ai+(ao-ai)*((1-A)*(z/len)+A*power(tan((pi*z)/(4*len)),rho));
        
%         plot(z, a);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         axis equal;
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Tangential Horn Profile', 'FontSize', 16 );
        
    case 5
        %%% x.rho profile %%%
        A = 1;
        rho = 2;
        a = ai+(ao-ai)*((1-A)*(z/len)+A*power(z/len,rho));
        
%         plot(z, a);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         axis equal;
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'xp Horn Profile', 'FontSize', 16 );
        
    case 6
        %%% Exponential profile %%%
        a=ai*exp(log(ao/ai)*(z/len));
        
%         plot(z, a);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         axis equal;
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Exponential Horn Profile', 'FontSize', 16 );
        
    case 7
        %%% Hyperbolic profile %%%
        a = sqrt(ai^2 + (power(z,2) * (ao^2-ai^2) / len^2));
        
%         plot(z, a);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         axis equal;
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Hyperbolic Horn Profile', 'FontSize', 16 );
        
    case 8
        %%% POLYNOMIAL Profile %%%
        rho = 3;
        a=ai+(rho+1)*(ao-ai)*(1-((rho*z)/((rho+1)*len))).*power(z/len,rho);
        
%         plot(z, a);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         axis equal;
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Polynomial Horn Profile', 'FontSize', 16 );
        
end

switch(mode_converter_type)
    
    case 1                            % case 1=VARIABLE SLOT DEPTH MODE CONVERTER
        % Mode Converter depths for element j
        ajmc = a(1:n_mc);                     % Index range for mode converter
        idx = 1:n_mc;
        djmc = (sigma-((idx-1)./n_mc).*(sigma-(0.25.*exp(1./(2.114.*(kc*ajmc).^1.134)))))*lam_c_mm;
        % Depth of remaining corrugations
        aj = a(n_mc+1:end);
        idx = n_mc+1:N+1;
        dj = ((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-n_mc-1)/(N-n_mc-1)).*((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lam_o_mm/4).*exp(1./(2.114.*(ko*ao).^1.134)));
        d = [djmc, dj];       % Combining the mode converter and horn depth values
        
        % Generate z,y coordinates as len and rad vector
        n = 0;
        lent(1) = 0;
        lent(2) = 0;
        for i = 1:N
            rad(i+n) = a(i);
            rad(i+n+1) = a(i)+d(i);
            rad(i+n+2) = a(i)+d(i);
            rad(i+n+3) = a(i+1);
            rad(i+n+4) = a(i+1);
            lent(i+n+2) = lent(i+n)+delta*p;
            lent(i+n+3) = lent(i+n+2);
            lent(i+n+4) = lent(i+n+3)+(1-delta)*p;
            lent(i+n+5) = lent(i+n+4);
            n = n+3;
        end
        
        z_number = (N*4)+1; % Number of coordinate points for corrugated length of horn
        lent = lent(1:z_number); % Truncate z axis data points to equal rad vector length
        
%         figure
%         plot(lent,rad);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Mode Converter and Corrugation Profile', 'FontSize', 16 );
%         axis equal;   % Scale axis equally for aspect ratio 1:1
        
    case 2                             % case 2=RING LOADED SLOT MODE CONVERTER
        % Mode Converter depths for element j
        ajmc = a(1:n_mc);                 % Index range for mode converter
        idx = 1:n_mc;
        djmc = (lam_c_mm/4).*exp(1./(2.114.*(kc*ajmc).^1.134));
        % Width of bjth slot for mode converter
        bj = (0.1+(idx-1).*((delta2-0.1)./n_mc)).*p;
        % Height of hjth slot for mode converter
        hj = (2/3).*djmc;
        
        % Depth of remaining corrugations
        aj = a(n_mc+1:end);
        idx = n_mc+1:N+1;
        dj = ((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-n_mc-1)/(N-n_mc-1)).*((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lam_o_mm/4).*exp(1./(2.114.*(ko*ao).^1.134)));
        d = [djmc, dj];       % Combining the mode converter and horn depth values
        
        % Generate z,y coordinates as lent,rad vector
        n = 5;
        lent = [0, 0, (-delta*p)+bj(1), (-delta*p)+bj(1), bj(1), bj(1)];
        rad = [a(1), a(1)+d(1)-hj(1), a(1)+d(1)-hj(1), a(1)+d(1), a(1)+d(1), a(2)];
        for i = 2:n_mc
            rad(i+n) = a(i);
            rad(i+n+1) = a(i)+d(i)-hj(i);
            rad(i+n+2) = a(i)+d(i)-hj(i);
            rad(i+n+3) = a(i)+d(i);
            rad(i+n+4) = a(i)+d(i);
            rad(i+n+5) = a(i+1);
            lent(i+n) = i*p-p;
            lent(i+n+1) = i*p-p;
            lent(i+n+2) = lent(i+n+1)-(delta*p)+bj(i);
            lent(i+n+3) = lent(i+n+1)-(delta*p)+bj(i);
            lent(i+n+4) = lent(i+n+3)+(delta*p)+bj(i);
            lent(i+n+5) = lent(i+n+3)+(delta*p)+bj(i);
            n = n+5;
        end
        % Add extra coordinate points before remaining corrugations
        lent(n_mc*(n_mc+1)+1) = lent(n_mc*(n_mc+1))+(1-delta)*p;
        rad(n_mc*(n_mc+1)+1) = a(n_mc+1);
%         
%         figure
%         plot(lent,rad);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Ring Loaded Mode Converter Profile', 'FontSize', 16 );
%         axis equal;   % Scale axis equally for aspect ratio 1:1
        
        n = n+n_mc+1;
        for i = n_mc+1:N
            rad(n) = a(i);
            rad(n+1) = a(i)+d(i);
            rad(n+2) = a(i)+d(i);
            rad(n+3) = a(i+1);
            rad(n+4) = a(i+1);
            lent(n+1) = lent(n);
            lent(n+2) = lent(n+1)+delta*p;
            lent(n+3) = lent(n+2);
            lent(n+4) = lent(n+3)+(1-delta)*p;
            n = n+4;
        end
        
        z_number = (n_mc*2)+(N*4)+1; % Number of coordinate points for corrugated length of horn
%         
%         figure
%         plot(lent,rad);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Ring Loaded Mode Converter Horn Profile', 'FontSize', 16 );
%         axis equal;   % Scale axis equally for aspect ratio 1:1
        
    case 3                 % case 3=VARIABLE PITCH TO WIDTH SLOT MODE CONVERTER
        % Mode Converter depths for element j
        ajmc = a(1:n_mc);                     % First indexes for mode converter
        idx = 1:n_mc;
        djmc = (sigma*(lam_c_mm/1.15)+((idx-1)./(n_mc-1)).*(lam_c_mm/4-(sigma*lam_c_mm/1.15))).*exp(1./(2.114.*(kc*ajmc).^1.134));
        % Depth of remaining corrugations
        aj = a(n_mc+1:end);
        idx = n_mc+1:N+1;
        dj = ((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-n_mc-1)/(N-n_mc-1)).*((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lam_o_mm/4).*exp(1./(2.114.*(ko*ao).^1.134)));
        d = [djmc, dj];       % Combining the mode converter and horn depth values
        
        % Generate z,y coordinates as lent,rad vector
        n = 0;
        lent(1) = 0;
        lent(2) = 0;
        for i = 1:n_mc
            rad(i+n) = a(i);
            rad(i+n+1) = a(i)+d(i);
            rad(i+n+2) = a(i)+d(i);
            rad(i+n+3) = a(i+1);
            rad(i+n+4) = a(i+1);
            lent(i+n+2) = lent(i+n)+(delta_min+(((i-1)./(n_mc-1))*(delta-delta_min)))*p;
            lent(i+n+3) = lent(i+n+2);
            lent(i+n+4) = lent(i+n+3)+(1-(delta_min+(((i-1)./(n_mc-1))*(delta-delta_min))))*p;
            lent(i+n+5) = lent(i+n+4);
            n = n+3;
        end
        
%         figure
%         plot(lent(1:end-1),rad);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Variable Pitch To Width Mode Converter Profile', 'FontSize', 16 );
%         axis equal;   % Scale axis equally for aspect ratio 1:1
        
        for i = n_mc+1:N
            rad(i+n) = a(i);
            rad(i+n+1) = a(i)+d(i);
            rad(i+n+2) = a(i)+d(i);
            rad(i+n+3) = a(i+1);
            rad(i+n+4) = a(i+1);
            lent(i+n+2) = lent(i+n)+delta*p;
            lent(i+n+3) = lent(i+n+2);
            lent(i+n+4) = lent(i+n+3)+(1-delta)*p;
            lent(i+n+5) = lent(i+n+4);
            n = n+3;
        end
        
        z_number = (N*4)+1; % Number of coordinate points for corrugated length of horn
        lent = lent(1:z_number); % Truncate z axis data points to equal rad vector length
        
%         figure
%         plot(lent,rad);
%         set(gca, "linewidth",2, "fontsize", 14 )
%         xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
%         ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
%         title( 'Variable Pitch To Width Mode Converter Horn Profile', 'FontSize', 16 );
%         axis equal;   % Scale axis equally for aspect ratio 1:1
        
end

% Add the rest of the geometry to create a closed path
% a_offset is the inner horn profile shifted up to give the horn a thickness
a_offset = a+(lam_c_mm/2+2);
%figure;                    % Uncomment these three lines for debugging
%plot(z, a_offset);
%axis equal;

% Add vertical surface at horn aperture
lent = [lent, lent(z_number)];
rad = [rad, a_offset(N)];
radmsh=rad;                 % radmesh to fix mesh lines to corrugations
%figure;                    % Uncomment these three lines for debugging
%plot(lent, rad);
%axis equal;

% Flip outer surface profile so that widest horn dimensions comes next in the outline coordinates
outer_surface = fliplr(a_offset);
z_flip = fliplr(z);
extent = lent(end);  % Fudge to make horn aperture planar for ring loaded slot MC
z_flip(1) = extent; % Fudge to make horn aperture planar for ring loaded slot MC
% Add outer profile and circular waveguide to horn
lent = [lent, z_flip, -wgl, -wgl, 0];
rad = [rad, outer_surface, ai+(lam_c_mm/2+2), ai, ai];
% 
% figure
% plot(lent,rad);
% set(gca, "linewidth",2, "fontsize", 14 )
% xlabel( 'Dimension in z Direction (mm)', 'FontSize', 14 );
% ylabel( 'Dimension in y Direction (mm)', 'FontSize', 14 );
% title( 'Complete Corrugated Horn Profile', 'FontSize', 16 );
% axis equal;   % Scale axis equally for aspect ratio 1:1

material = 'Copper (annealed)';
Name = 'Waveguide';
OuterRadius = ao+(lam_c_mm/2+2);
InnerRadius = ai;
Xcenter = 0;
Ycenter = 0;
Zrange = [-wgl 0];
Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)

switch(mode_converter_type)
    
    case 1
        
        n = 0;
        lent(1) = 0;
        lent(2) = 0;
        for i=1:N
            rad(i+n) = a(i);
            rad(i+n+1) = a(i)+d(i);
            rad(i+n+2) = a(i)+d(i);
            rad(i+n+3) = a(i+1);
            rad(i+n+4) = a(i+1);
            lent(i+n+2) = lent(i+n)+delta*p;
            lent(i+n+3) = lent(i+n+2);
            lent(i+n+4) = lent(i+n+3)+(1-delta)*p;
            lent(i+n+5) = lent(i+n+4);
            
            material = 'Copper (annealed)';
            Name = ['Coor' num2str(i) 'b'];
            OuterRadius = ao+(lam_c_mm/2+2);
            InnerRadius = rad(i+n+1);
            Xcenter = 0;
            Ycenter = 0;
            Zrange = [lent(i+n) lent(i+n+2)];
            Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
            
            material = 'Copper (annealed)';
            Name = ['Coor' num2str(i) 'a'];
            OuterRadius = ao+(lam_c_mm/2+2);
            InnerRadius = rad(i+n+3);
            Xcenter = 0;
            Ycenter = 0;
            Zrange = [lent(i+n+3) lent(i+n+5)];
            Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
            
            n = n+3;
        end
        
    case 2                             % case 2=RING LOADED SLOT MODE CONVERTER
        % Mode Converter depths for element j
        ajmc = a(1:n_mc);                 % Index range for mode converter
        idx = 1:n_mc;
        djmc = (lam_c_mm/4).*exp(1./(2.114.*(kc*ajmc).^1.134));
        % Width of bjth slot for mode converter
        bj = (0.1+(idx-1).*((delta2-0.1)./n_mc)).*p;
        % Height of hjth slot for mode converter
        hj = (2/3).*djmc;
        
        % Depth of remaining corrugations
        aj = a(n_mc+1:end);
        idx = n_mc+1:N+1;
        dj = ((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134)))-((idx-n_mc-1)/(N-n_mc-1)).*((lam_c_mm/4).*exp(1./(2.114.*(kc*aj).^1.134))-(lam_o_mm/4).*exp(1./(2.114.*(ko*ao).^1.134)));
        d = [djmc, dj];       % Combining the mode converter and horn depth values
        
        % Generate z,y coordinates as lent,rad vector
        n = 5;
        lent = [0, 0, (-delta*p)+bj(1), (-delta*p)+bj(1), bj(1), bj(1)];
        rad = [a(1), a(1)+d(1)-hj(1), a(1)+d(1)-hj(1), a(1)+d(1), a(1)+d(1), a(2)];
        i=1;
        material = 'Copper (annealed)';
        Name = ['Coor' num2str(1) 'c'];
        OuterRadius = rad(i+1);
        InnerRadius = rad(i);
        Xcenter = 0;
        Ycenter = 0;
        Zrange = [lent(i+1) lent(i+2)];
        Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
        material = 'Copper (annealed)';
        
        Name = ['Coor' num2str(1) 'd'];
        OuterRadius = ao+(lam_c_mm/2+2);
        InnerRadius = rad(i+3);
        Xcenter = 0;
        Ycenter = 0;
        Zrange = [lent(i+3) lent(i+4)];
        Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
        
        material = 'Copper (annealed)';
        Name = ['Coor' num2str(1) 'e'];
        OuterRadius = ao+(lam_c_mm/2+2);
        InnerRadius = rad(i+5);
        Xcenter = 0;
        Ycenter = 0;
        Zrange = [lent(i+4) lent(i+4)+p-(delta*p)];
        Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
        
        material = 'Copper (annealed)';
        Name = ['Coor' num2str(1) 'f'];
        OuterRadius = rad(i+3);
        InnerRadius = rad(i+1);
        Xcenter = 0;
        Ycenter = 0;
        Zrange = [lent(i+1) lent(i+2)];
        Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
        
        component1 = 'component1:Waveguide';
        component2 = 'component1:Coor1f';
        CstSubtract(mws,component1,component2)
        
        for i = 2:n_mc
            rad(i+n) = a(i);
            rad(i+n+1) = a(i)+d(i)-hj(i);
            rad(i+n+2) = a(i)+d(i)-hj(i);
            rad(i+n+3) = a(i)+d(i);
            rad(i+n+4) = a(i)+d(i);
            rad(i+n+5) = a(i+1);
            lent(i+n) = i*p-p;
            lent(i+n+1) = i*p-p;
            lent(i+n+2) = lent(i+n+1)-(delta*p)+bj(i);
            lent(i+n+3) = lent(i+n+1)-(delta*p)+bj(i);
            lent(i+n+4) = lent(i+n+3)+(delta*p)+bj(i);
            lent(i+n+5) = lent(i+n+3)+(delta*p)+bj(i);
            
            material = 'Copper (annealed)';
            Name = ['Coor' num2str(i) 'c'];
            OuterRadius = rad(i+n+1);
            InnerRadius = rad(i+n);
            Xcenter = 0;
            Ycenter = 0;
            Zrange = [lent(i+n+1) lent(i+n+2)];
            Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
            material = 'Copper (annealed)';
            
            Name = ['Coor' num2str(i) 'd'];
            OuterRadius = ao+(lam_c_mm/2+2);
            InnerRadius = rad(i+n+3);
            Xcenter = 0;
            Ycenter = 0;
            Zrange = [lent(i+n+3) lent(i+n+4)];
            Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
            
            material = 'Copper (annealed)';
            Name = ['Coor' num2str(i) 'e'];
            OuterRadius = ao+(lam_c_mm/2+2);
            InnerRadius = rad(i+n+5);
            Xcenter = 0;
            Ycenter = 0;
            Zrange = [lent(i+n+4) lent(i+n+4)+(1-delta)*p];
            Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
            
            n = n+5;
        end
        % Add extra coordinate points before remaining corrugations
        lent(n_mc*(n_mc+1)+1) = lent(n_mc*(n_mc+1))+(1-delta)*p;
        rad(n_mc*(n_mc+1)+1) = a(n_mc+1);
        
        n = n+n_mc+1;
        for i = n_mc+1:N
            rad(n) = a(i);
            rad(n+1) = a(i)+d(i);
            rad(n+2) = a(i)+d(i);
            rad(n+3) = a(i+1);
            rad(n+4) = a(i+1);
            lent(n+1) = lent(n);
            lent(n+2) = lent(n+1)+delta*p;
            lent(n+3) = lent(n+2);
            lent(n+4) = lent(n+3)+(1-delta)*p;
            
            material = 'Copper (annealed)';
            Name = ['Coor' num2str(i) 'b'];
            OuterRadius = ao+(lam_c_mm/2+2);
            InnerRadius = rad(n+1);
            Xcenter = 0;
            Ycenter = 0;
            Zrange = [lent(n) lent(n+2)];
            Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
            
            material = 'Copper (annealed)';
            Name = ['Coor' num2str(i) 'a'];
            OuterRadius = ao+(lam_c_mm/2+2);
            InnerRadius = rad(n+3);
            Xcenter = 0;
            Ycenter = 0;
            Zrange = [lent(n+2) lent(n+4)];
            Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
            n = n+4;
        end
        
        
    case 3
        n = 0;
        lent(1) = 0;
        lent(2) = 0;
        for i = 1:n_mc
            rad(i+n) = a(i);
            rad(i+n+1) = a(i)+d(i);
            rad(i+n+2) = a(i)+d(i);
            rad(i+n+3) = a(i+1);
            rad(i+n+4) = a(i+1);
            lent(i+n+2) = lent(i+n)+(delta_min+(((i-1)./(n_mc-1))*(delta-delta_min)))*p;
            lent(i+n+3) = lent(i+n+2);
            lent(i+n+4) = lent(i+n+3)+(1-(delta_min+(((i-1)./(n_mc-1))*(delta-delta_min))))*p;
            lent(i+n+5) = lent(i+n+4);
            
            material = 'Copper (annealed)';
            Name = ['Coor' num2str(i) 'b'];
            OuterRadius = ao+(lam_c_mm/2+2);
            InnerRadius = rad(i+n+1);
            Xcenter = 0;
            Ycenter = 0;
            Zrange = [lent(i+n) lent(i+n+2)];
            Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
            
            material = 'Copper (annealed)';
            Name = ['Coor' num2str(i) 'a'];
            OuterRadius = ao+(lam_c_mm/2+2);
            InnerRadius = rad(i+n+3);
            Xcenter = 0;
            Ycenter = 0;
            Zrange = [lent(i+n+3) lent(i+n+5)];
            Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
            n = n+3;
        end
        
        for i = n_mc+1:N
            rad(i+n) = a(i);
            rad(i+n+1) = a(i)+d(i);
            rad(i+n+2) = a(i)+d(i);
            rad(i+n+3) = a(i+1);
            rad(i+n+4) = a(i+1);
            lent(i+n+2) = lent(i+n)+delta*p;
            lent(i+n+3) = lent(i+n+2);
            lent(i+n+4) = lent(i+n+3)+(1-delta)*p;
            lent(i+n+5) = lent(i+n+4);
            
            material = 'Copper (annealed)';
            Name = ['Coor' num2str(i) 'b'];
            OuterRadius = ao+(lam_c_mm/2+2);
            InnerRadius = rad(i+n+1);
            Xcenter = 0;
            Ycenter = 0;
            Zrange = [lent(i+n) lent(i+n+2)];
            Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
            
            material = 'Copper (annealed)';
            Name = ['Coor' num2str(i) 'a'];
            OuterRadius = ao+(lam_c_mm/2+2);
            InnerRadius = rad(i+n+3);
            Xcenter = 0;
            Ycenter = 0;
            Zrange = [lent(i+n+3) lent(i+n+5)];
            Cstcylinder(mws, Name, 'component1', material, 'Z', OuterRadius, InnerRadius, Xcenter, Ycenter, Zrange)
            n = n+3;
        end
end

guidata(hObject, handles);


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popupmenu2.
function popupmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu2


% --- Executes during object creation, after setting all properties.
function popupmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
