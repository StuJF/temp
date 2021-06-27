function handles = topo3d(Values,chanlocs,EEG3dset,varargin)
%topo3d plot 3-D topography of chanlocs/spatial weights
%   Functions used are from Brainstorm software: https://neuroimage.usc.edu/brainstorm   
%   handles = topo3d(Values,chanlocs,EEG3dset,varargin)
%   Input:
%      -Values(needed): []/ nchan X 1 vector if Values is
%               empty([]), then topo3d plots EEG channel locations; if Values
%               is nchan by 1 vector(where nchan is the number of channels,length(chanlocs)), 
%               then topo3d plots spatial weights topography;
%      -chanlocs(needed): Channel locations structure. Using EEGlab EEG.chanlocs or
%               location structure S containing:
%                   S.labels,S.theta,S.radiu,S.X,S.Y,S.Z
%      -EEG3dset(needed but the programme load it by default): 
%               Settings of topography, which has been provided as 'EEG3dset'.
%               topo3d will load it if exists, or you need to provide it
%   Options:
%      -CLim(default [0,1]):color limits of topography, 2 by 1 vector
%      -electrodes(default 'on'): 'off': don't show electrodes; 'on': show
%               electordes; 'numbers': show electrodes and channel serial
%               number; 'text': show electrodes and channel labels
%      -plotchans(default is 1:length(chanlocs)): channel numbers to plot,vector 
%      -headColor(default is [1 0.9 0.8]): color of cartoon head and ear; 3 by
%                1 vector.
%      -electrodeSize(default is 5): size of electrode. float
%      -electrodColor(default is [1 0.7 0.3]): color of electrode
%      -textSize(default is 10): size of electrode numbers/labels
%      -textColor(default is [0 0 0]): color of electrode numbers/labels.(see Options:electrodes)
%      -colormap(default is jet):colormap, 
%      -Axes: default is gca; plot topography on specific axis if provided.

if nargin<3||isempty(EEG3dset)
    try
        path2file = which('topo3d');
        path2set = fileparts(path2file);
        path2set = fullfile(path2set,'EEG3dset');
        load(path2set,'EEG3dset');
    catch
        error('No 3d setting found');
    end
end
% default setting
ELECTRODES = 'on';
headColor = [1 0.9 0.8];
electrodeColor = [1 0.7 0.3];
textColor = [0 0 0];
CLim = [0 1];
colorMap = CBar;
electrodeSize = 5;
plotchans = 1:length(chanlocs);
textSize = 10;
ax2plot = [];

nargs = nargin;
if nargs > 3
    for i = 1:2:length(varargin)
        Param = lower(varargin{i});
        Value = varargin{i+1};
        switch Param
            case 'electrodes'
                ELECTRODES = lower(Value);
            case 'headcolor'
                headColor = Value;
                if size(Value,2)~=3
                    error('Head color must be a 1 x 3 matrix')
                end
            case 'electrodecolor'
                electrodeColor = Value;
                if size(Value,2)~=3
                    error('Electrode color must be a 1 x 3 matrix')
                end
            case 'textcolor'
                textColor = Value;
                if size(Value,2)~=3
                    error('Head color must be a 1 x 3 matrix')
                end
            case 'colormap'
                colorMap = Value;
                if size(Value,2)~=3
                    error('Colormap must be a n x 3 matrix')
                end
            case 'electrodesize'
                electrodeSize = Value;
            case 'clim'
                CLim = Value;
                if numel(CLim)~=2
                    error('Color Limit must be 1 x 2 vector');
                end
            case 'plotchans'
                plotchans = Value;
            case 'textsize'
                textSize = Value;
            case 'axes'
                ax2plot = Value;
        end
    end
end

headvertices = EEG3dset.HeadModel.Vertices;
headfaces = EEG3dset.HeadModel.Faces;
chanmark = zeros(length(chanlocs),3);
rejChan = [];
for ic = 1:length(chanlocs)
    fname = chanlocs(ic).labels;
    for jc = 1:length(necChanloc)
        if strcmpi(EEG3dChannel.Name,fname)
            chanmark(ic,:) = EEG3dChannel.Loc;
            break;
        end
        rejChan = [rejChan,ic];
    end
end

% select channel to plot
plotchans = setdiff(plotchans,rejChan);
chanmark = chanmark(plotchans,:);
% project to head surface
chanmark = channel_project_scalp(headvertices,chanmark);
if isempty(ax2plot)
    handles.h = figure('color',[1 1 1]);
    handles.axes = axes('Parent',handles.h);
else
    handles.axes = ax2plot;
end
hold on
% plot 3d head
handles.hs = patch(handles.axes,...
    'Faces',            headfaces, ...
    'Vertices',         headvertices,...
    'FaceVertexCData',  [], ...
    'FaceColor',        headColor, ...
    'FaceAlpha',        1 , ...
    'AlphaDataMapping', 'none', ...
    'EdgeColor',        'none', ...
    'EdgeAlpha',        1, ...
    'BackfaceLighting', 'lit', ...
    'AmbientStrength',  0.5, ...
    'DiffuseStrength',  0.5, ...
    'SpecularStrength', 0.2, ...
    'SpecularExponent', 1, ...
    'SpecularColorReflectance', 0.5, ...
    'FaceLighting',     'gouraud', ...
    'EdgeLighting',     'gouraud',...
    'hittest','off');
material dull
set(handles.axes,'Xcolor',[1 1 1],'ycolor',[1 1 1],'zcolor',[1 1 1],'CameraPosition', [0 100 0],'CameraViewAngle',6);
view([ -90 90 ])

axis equal
axlim=[min(headvertices);max(headvertices)]*1.2;
axis(reshape(axlim,1,[]));

hr = rotate3d;
LightSource = camlight('headlight');set(LightSource,'style','infinite');
set(hr,'Enable','on','ActionPostCallback',{@move_light_source,LightSource});

if strcmp(ELECTRODES,'off')
else
    % electrodes
    [sx,sy,sz] = sphere(16);
    sscale = 1000/electrodeSize;
    for ichan = 1:size(chanmark,1)
        surf(sx/sscale+chanmark(ichan,1),sy/sscale+chanmark(ichan,2),sz/sscale+chanmark(ichan,3),'linestyle','none','facecolor',electrodeColor);
    end
    if strcmp(ELECTRODES,'labels')
        %labels
        X = 1.05*chanmark(:,1)+sign(chanmark(:,1))/sscale;
        Y = 1.1*chanmark(:,2);
        Z = 1.05*chanmark(:,3)+1/sscale;
        text(X, Y, Z, ...
            {chanlocs(plotchans).labels}, ...
            'Parent',              handles.axes, ...
            'HorizontalAlignment', 'center', ...
            'FontSize',            textSize, ...
            'FontUnits',           'points', ...
            'FontWeight',          'bold', ...
            'Tag',                 'SensorsLabels', ...
            'Color',               textColor, ...
            'hittest','off','HorizontalAlignment','center');
    elseif strcmp(ELECTRODES,'numbers')
        X = 1.05*chanmark(:,1)+sign(chanmark(:,1))/sscale;
        Y = 1.1*chanmark(:,2);
        Z = 1.05*chanmark(:,3)+1/sscale;
        text(X, Y, Z, ...
            {num2str(plotchans(:))}, ...
            'Parent',              handles.axes, ...
            'HorizontalAlignment', 'center', ...
            'FontSize',            textSize, ...
            'FontUnits',           'points', ...
            'FontWeight',          'bold', ...
            'Tag',                 'SensorsLabels', ...
            'Color',               textColor, ...
            'hittest','off','HorizontalAlignment','center');
    end
    
end


if ~isempty(Values)
    Values = Values(plotchans);
    valueVertices = chanmark;
    iZero = find(all(abs(valueVertices) < 1/1e4, 2));
    if ~isempty(iZero)
        valueVertices(iZero,:) = [];
    end
    % Tesselate sensor cap
    valueFaces = channel_tesselate(valueVertices, 1);
    % Remove some  pathological triangles
    valueFaces = tess_threshold(valueVertices, valueFaces, 5, 3, 170);
    % Refine mesh
    [valueVertices, valueFaces] = tess_refine(valueVertices, valueFaces, [], [], 1);
    [valueVertices, valueFaces] = tess_refine(valueVertices, valueFaces, [], [], 1);
    valueVertices = channel_project_scalp(headvertices*1.0265,valueVertices);
    
    WExtrap = bst_shepards(valueVertices, chanmark, 12, 0, 1);
    
    % spatial weights
    scalpMap = WExtrap*Values;
    scalpMap = (scalpMap - min(CLim))/abs(diff(CLim));
    
    colorMap = colorMap(floor(min(1,scalpMap)*(length(colorMap)-1)+1),:);
    handles.hSurf = patch(handles.axes,...
        'Faces',            valueFaces, ...
        'Vertices',         valueVertices, ...
        'EdgeColor',        'none', ...
        'FaceColor',        'interp', ...
        'FaceVertexCData',  colorMap, ...
        'BackfaceLighting', 'lit', ...
        'AmbientStrength',  0.95, ...
        'DiffuseStrength',  0, ...
        'SpecularStrength', 0, ...
        'FaceLighting',     'gouraud', ...
        'EdgeLighting',     'gouraud', ...
        'Parent',           handles.axes, ...
        'Tag',              'TopoSurface');
    
end
caxis(CLim);
end

function move_light_source(src,evt,LightSource)
LightSource = camlight(LightSource,'headlight');
end


%%% Project channel to head surface
function ChanLoc = channel_project_scalp(Vertices, ChanLoc)
% CHANNEL_ALIGN_MANUAL: Align manually an electrodes net on the scalp surface of the subject.
%
% INPUT:
%     - Vertices : [Mx3] positions of the scalp vertices
%     - ChanLoc  : [Nx3] positions of the EEG electrodes
%
% Authors: Francois Tadel, 2014

% Center the surface on its center of mass
center = mean(Vertices, 1);
Vertices = bsxfun(@minus, Vertices, center);
% Parametrize the surface
p   = .2;
th  = -pi-p   : 0.01 : pi+p;
phi = -pi/2-p : 0.01 : pi/2+p;
rVertices = tess_parametrize_new(Vertices, th, phi);

% Process each sensor
for iChan = 1:size(ChanLoc,1)
    % Get the closest surface from the point
    c = ChanLoc(iChan,:);
    % Center electrode
    c = c - center;
    % Convert in spherical coordinates
    [c_th,c_phi,c_r] = cart2sph(c(1), c(2), c(3));
    % Interpolate
    c_r = interp2(th, phi, rVertices, c_th, c_phi);
    % Project back in cartesian coordinates
    [c(1),c(2),c(3)] = sph2cart(c_th, c_phi, c_r);
    % Restore initial origin
    c = c + center;
    ChanLoc(iChan,:) = c;
end
end

function [fr,th_val,phi_val] = tess_parametrize_new(vert, th_val, phi_val, p)
% TESS_PARAMETRIZE_NEW: Get a parametric representation of a closed and non-overlapping surface.
%
% USAGE:  tess_parametrize_new(vert, th_val, phi_val, p)
%         tess_parametrize_new(vert, th_val, phi_val)
%
% INPUTS:
%    - vert    : x,y,z coordinates of the surface to parametrize
%    - th_val  : Row array of theta values
%    - phi_val : Row array of phi values
%    - p       : Pad the function of p values in both direction with circular values


% Parse inputs
if (nargin < 4) || isempty(p)
    p = 0;
end

% Convert surface to spherical coordinates
[th,phi,r] = cart2sph(vert(:,1), vert(:,2), vert(:,3));
% Replicate the values to get a full coverage at the edges
th  = [th;  th;        th;        th+2*pi;  th+2*pi;  th+2*pi;   th-2*pi;  th-2*pi;  th-2*pi  ];
phi = [phi; phi+2*pi;  phi-2*pi;  phi;      phi+2*pi; phi-2*pi;  phi;      phi+2*pi; phi-2*pi ];
r   = [r;   r;         r;         r;        r;        r;         r;        r;        r        ];

% Build grid of (theta,phi) values
[fth, fphi] = meshgrid(th_val, phi_val);
% Estimate radius values for all those points
fr = griddata(th, phi, r, fth, fphi);

% Pad function with circular values, to avoid edge irregularities
if (p > 0)
    fr_pad  = bst_pad(fr, [0  p], 'circular');
    fr_pad  = bst_pad(fr_pad, [p  0], 'mirror');
    fr = fr_pad;
    th_val  = [th_val(end-p+1:end)-2*pi,  th_val,  th_val(1:p)+2*pi];
    phi_val = [phi_val(end-p+1:end)-2*pi, phi_val, phi_val(1:p)+2*pi];
end
end


function [az,elev,r] = cart2sph(x,y,z)
%CART2SPH Transform Cartesian to spherical coordinates.
%   [TH,PHI,R] = CART2SPH(X,Y,Z) transforms corresponding elements of
%   data stored in Cartesian coordinates X,Y,Z to spherical
%   coordinates (azimuth TH, elevation PHI, and radius R).  The arrays
%   X,Y, and Z must be the same size (or any of them can be scalar).
%   TH and PHI are returned in radians.
%
%   TH is the counterclockwise angle in the xy plane measured from the
%   positive x axis.  PHI is the elevation angle from the xy plane.
%
%   Class support for inputs X,Y,Z:
%      float: double, single
%
%   See also CART2POL, SPH2CART, POL2CART.

%   Copyright 1984-2005 The MathWorks, Inc.

hypotxy = hypot(x,y);
r = hypot(hypotxy,z);
elev = atan2(z,hypotxy);
az = atan2(y,x);
end

function Apad = bst_pad(A, p, method)
% BST_PAD: Pad a 2D array with circular/mirrored values, or zeros
%
% USAGE:  Apad = bst_pad(A, p, method)
%
% INPUT:
%    - A      : Array to pad
%    - p      : Number of values to add before and after, in each direction
%               If one value: use the same for dimensions 1 and 2
%               If two values: pad dimensions 1 and 2 with a different amount of values
%    - method : {'zeros', 'mirror', 'circular'}


% Parse inputs
if (length(p) == 1)
    px = p;
    py = p;
elseif (length(p) == 2)
    px = p(1);
    py = p(2);
else
    error('Invalid call');
end

% Initialize returned array
Apad = zeros(size(A) + 2*[px py]);
Apad(px+1:end-px, py+1:end-py) = A;
% Padding method
switch lower(method)
    case 'zeros'
        % Nothing to do...
    case 'circular'
        if (px > size(A,1)) || (py > size(A,2))
            error('Paddind size exceeds limit for "circular" method.');
        end
        % Padding x
        Apad(1:px, py+1:end-py) = A(end-px+1:end, :);
        Apad(end-px+1:end, py+1:end-py) = A(1:px,:);
        % Padding y
        Apad(:, 1:py) = Apad(:, end-2*py+1:end-py);
        Apad(:, end-py+1:end) = Apad(:,py+1:2*py);
    case 'mirror'
        if (px > size(A,1)-1) || (py > size(A,2)-1)
            error('Paddind size exceeds limit for "mirror" method.');
        end
        % Padding x
        Apad(1:px, py+1:end-py) = A(px+1:-1:2, :);
        Apad(end-px+1:end, py+1:end-py) = A(end-1:-1:end-px,:);
        % Padding y
        Apad(:, 1:py) = Apad(:, 2*py+1:-1:py+2);
        Apad(:, end-py+1:end) = Apad(:, end-py-1:-1:end-2*py);
    otherwise
        error('Unknown method');
end
end

%%% Plot channels

function Faces = channel_tesselate( Vertices, isPerimThresh )
% CHANNEL_TESSELATE: Tesselate a set of EEG or MEG sensors, for display purpose only.
%
% USAGE:  Faces = channel_tesselate( Vertices, isPerimThresh=1 )
%
% INPUT:
%    - Vertices      : [Nx3], set of 3D points (MEG or EEG sensors)
%    - isPerimThresh : If 1, remove the Faces that are too big
% OUTPUT:
%    - Faces    : [Mx3], result of the tesselation

% @=============================================================================
% This function is part of the Brainstorm software:
% https://neuroimage.usc.edu/brainstorm


% Parse inputs
if (nargin < 2) || isempty(isPerimThresh)
    isPerimThresh = 1;
end

% === TESSELATE ===
% 2D Projection
[X,Y] = bst_project_2d(Vertices(:,1), Vertices(:,2), Vertices(:,3), '2dcap');
% Compute best fitting sphere
bfs_center = bst_bfs(Vertices)';
% Center Vertices on BFS center
coordC = bsxfun(@minus, Vertices, bfs_center);
% Normalize coordinates
coordC = bsxfun(@rdivide, coordC, sqrt(sum(coordC.^2,2)));
coordC = bsxfun(@rdivide, coordC, sqrt(sum(coordC.^2,2)));
% Tesselation of the sensor array
Faces = convhulln(coordC);


% === REMOVE UNNECESSARY TRIANGLES ===
% For instance: the holes for the ears on high-density EEG caps
if isPerimThresh
    % Get border of the representation
    border = convhull(X,Y);
    % Keep Faces inside the border
    iInside = find(~(ismember(Faces(:,1),border) & ismember(Faces(:,2),border)& ismember(Faces(:,3),border)));
    %Faces   = Faces(iInside, :);
    
    % Compute perimeter
    triPerimeter = tess_perimeter(Vertices, Faces);
    % Threshold values
    thresholdPerim = mean(triPerimeter(iInside)) + 6 * std(triPerimeter(iInside));
    % Apply threshold
    iFacesOk = intersect(find(triPerimeter <= thresholdPerim), iInside);
    % Find Vertices that are not in the Faces matrix
    iVertNotInFaces = setdiff(1:length(Vertices), unique(Faces(:)));
    if ~isempty(iVertNotInFaces)
        % disp(['CHANNEL_TESSELATE> WARNING: Some sensors are not in the Faces list: ' sprintf('%d ', iVertNotInFaces)]);
    end
    % Loop until all the Vertices are visible
    isMissing = 1;
    while isMissing
        % List all the Vertices ignored by the reduced mesh
        iVertOk = unique(reshape(Faces(iFacesOk,:),[],1));
        iVertMissing = setdiff(1:length(Vertices), iVertOk);
        iVertMissing = setdiff(iVertMissing, iVertNotInFaces);
        % If all the Vertices are included, next step
        if isempty(iVertMissing)
            isMissing = 0;
        else
            % Find Faces connected to the first missing vertex
            iFacesAdd = find(any(Faces == iVertMissing(1), 2));
            % From the potential candidate Faces, keep the one that has the smaller perimeter
            [minP, iMinP] = min(triPerimeter(iFacesAdd));
            % Add the smallest face to the list of Faces we keep
            iFacesOk(end+1) = iFacesAdd(iMinP);
        end
    end
    % Remove the Faces
    Faces = Faces(sort(iFacesOk),:);
end

end

function [X,Y] = bst_project_2d(x, y, z, Method)
% BST_PROJECT_2D: Project a set of 3D points (EEG or MEG sensors) on a 2D surface.
%
% USAGE:  [X,Y] = bst_project_2d(x, y, z, Method='2dcap');


% Parse inputs
if (nargin < 4) || isempty(Method)
    Method = '2dcap';
end

% Different projection methods
switch (Method)
    case '2dcap'
        % Spherical coordinates
        [TH,PHI,R] = cart2sph(x, y, z);
        % Flat projection
        R = 1 - PHI ./ pi*2;
        % Convert back to cartesian coordinates
        [X,Y] = pol2cart(TH,R);
        
    case '2dlayout'
        % Spherical coordinates
        z = z - max(z);
        [TH,PHI,R] = cart2sph(x, y, z);
        % Remove the too smal values for PHI
        PHI(PHI < 0.001) = 0.001;
        % Flat projection
        R2 = R ./ cos(PHI) .^ .2;
        [X,Y] = pol2cart(TH,R2);
        
    case 'circle'
        % Spherical coordinates
        [TH,PHI,R] = cart2sph(x, y, z);
        % Flat projection
        R = 1 - PHI ./ pi*2;
        % Convert back to cartesian coordinates
        [X,Y] = pol2cart(TH,R);
        
        % Convert back to cartesian coordinates
        [TH,R] = cart2pol(X,Y);
        % Convex hull
        facesBorder = convhull(X,Y);
        iBorder = unique(facesBorder);
        
        % Deformation field in radius computed from the border sensors, projected onto the circle
        Rcircle = 1;
        Dborder = Rcircle ./ R(iBorder);
        % Compute the radius deformation to apply to all the sensors
        funcTh = [TH(iBorder)-2*pi; TH(iBorder); TH(iBorder)+2*pi];
        funcD  = [Dborder; Dborder; Dborder];
        D = interp1(funcTh, funcD, TH, 'linear', 0);
        
        % Remove the possible zero values
        D(D == 0) = 1;
        % Compute new radius: the closer to the center, the less transformed
        R = min(R .* D, Rcircle);
        % Convert back to cartesian coordinates
        [X,Y] = pol2cart(TH, R);
        
end
end

function [ HeadCenter, Radius ] = bst_bfs( Vertices )
% BST_BFS: Compute best fitting sphere (BFS) for a set of points.
%
% USAGE:  [ HeadCenter, Radius ] = bst_bfs( Vertices )
%
% INPUT:
%     - Vertices  : [Nv,3] double, (x,y,z) coordinates of the points
% OUTPUT:
%     - HeadCenter: (x,y,z) coordinates of the BFS
%     - Radius    : Radius of the BFS

% Check matrices orientation
if (size(Vertices, 2) ~= 3)
    error('Vertices must have 3 columns (X,Y,Z).');
end
% Convert vertices into double
Vertices = double(Vertices);
% 500 points is more than enough to compute scalp's best fitting sphere
nscalp = size(Vertices,1);
if (nscalp > 500)
    Vertices = Vertices(unique(round(linspace(1,nscalp,500))),:);
end

% Center of mass of the scalp vertex locations
mass = mean(Vertices);
% Average distance between the center of mass and the scalp points
diffvert = bst_bsxfun(@minus, Vertices, mass);
R0 = mean(sqrt(sum(diffvert.^2, 2)));

% Optimization
vec0 = [mass,R0];
minn = fminsearch(@dist_sph, vec0, [], Vertices);

% Results : Center
HeadCenter = minn(1:end-1)'; % 3x1
% Largest radius (largest sphere radius)
Radius = minn(end);

end


% FMINS distance function used to minimize the fit to a sphere.
% Given center and list of sensor points, find the average distance from the center to these points
function d = dist_sph(vec,sensloc)
R = vec(end);
center = vec(1:end-1);
% Average distance between the center if mass and the electrodes
diffvert = bst_bsxfun(@minus, sensloc, center);
d = mean(abs(sqrt(sum(diffvert.^2,2)) - R));
end

function [Faces, iFacesRemove] = tess_threshold(Vertices, Faces, threshArea, threshRatio, threshAngle, threshEdge)
% TESS_THRESHOLD: Detect pathological triangles in a surface mesh.
%
% INPUTS:
%    - Vertices    : [Nvert x 3] surface vertices
%    - Faces       : [Nfaces x 3] surface triangles
%    - threshArea  : Detect large triangles (area > thresh * std)
%    - threshRatio : Detect asymetric triangles (ratio perimeter/area > thresh * std)
%    - threshAngle : Detect triangles with angles that are too open (angle in degrees > thresh)

%
% Authors: Francois Tadel, 2012-2016

% ===== PARSE INPUTS =====
if (nargin < 6) || isempty(threshEdge)
    threshEdge = [];
end
if (nargin < 5) || isempty(threshAngle)
    threshAngle = [];
end
if (nargin < 4) || isempty(threshRatio)
    threshRatio = [];
end
if (nargin < 3)
    threshArea = [];
end
iFacesRemoveArea  = [];
iFacesRemoveRatio = [];
iFacesRemoveAngle = [];
iFacesRemoveEdge  = [];

% ===== COMPUTE SURFACE STATISTICS =====
% Triangles area
triArea = tess_area(Vertices, Faces);
% Compute the vector of each edge
v1 = Vertices(Faces(:,1),:) - Vertices(Faces(:,2),:);
v2 = Vertices(Faces(:,1),:) - Vertices(Faces(:,3),:);
v3 = Vertices(Faces(:,2),:) - Vertices(Faces(:,3),:);

% ===== THRESHOLD: AREA =====
% Detect the faces that have an area above the threshold
if ~isempty(threshArea) && (threshArea > 0)
    iFacesRemoveArea = find(triArea - mean(triArea) > threshArea * std(triArea));
end

% ===== THRESHOLD: PERIMETER/AREA =====
if ~isempty(threshRatio) && (threshRatio > 0)
    % Compute perimeter again
    triPerimeter = tess_perimeter(Vertices, Faces);
    % Ratio perimeter / area
    ratio = (triPerimeter ./ triArea);
    % Detect the Faces that have an area above the threshold
    iFacesRemoveRatio = find(ratio - mean(ratio) > threshRatio * std(ratio));
end

% ===== THRESHOLD: ANGLE =====
if ~isempty(threshAngle) && (threshAngle > 0)
    % Compute the angle between all the vectors
    maxAngle = zeros(size(Vertices,1),1);
    for i = 1:size(v1,1)
        maxAngle(i) = max([atan2(norm(cross(v1(i,:),v2(i,:))), dot(v1(i,:),v2(i,:))), ...
            atan2(norm(cross(v1(i,:),v3(i,:))), dot(v1(i,:),v3(i,:))), ...
            atan2(norm(cross(v2(i,:),v3(i,:))), dot(v2(i,:),v3(i,:)))]);
    end
    % Convert to degrees
    maxAngle = maxAngle / 2 / pi * 360;
    % Detect the Faces that have an area above the threshold
    iFacesRemoveAngle = find(maxAngle > threshAngle);
end

% ===== THRESHOLD: EDGE LENGTH =====
if ~isempty(threshEdge) && (threshEdge > 0)
    % Compute the length of all the edges
    edgeLength = sqrt(v1.^2 + v2.^2 + v3.^2);
    % Split long edges
    iFacesRemoveEdge = find(edgeLength - mean(edgeLength) > threshEdge * std(edgeLength));
end

% List of faces to remove
iFacesRemove = [iFacesRemoveArea(:); iFacesRemoveRatio(:); iFacesRemoveAngle(:)];
% Keep only the good faces
if ~isempty(iFacesRemove)
    Faces(iFacesRemove,:) = [];
end

end

function [Vertices, Faces] = tess_refine(Vertices, Faces, threshArea, threshEdge, isThresh )
% TESS_REFINE: Refine a surface mesh.
%
% USAGE:  [Vertices, Faces] = tess_refine(Vertices, Faces, threshArea=[], threshEdge=[], isThresh=0);
%
% INPUT:
%     - Vertices   : Mx3 double matrix
%     - Faces      : Nx3 double matrix
%     - threshArea : Only split the faces that have an area above a certain threshold (edge > thresh * std)
%     - threshEdge : Only split the edges that have a length above the threshold (edge > thresh * std)
%     - isThresh   : If 1, call channel_tesselate with a isThresh = 1 (remove big triangles)
%
% DESCRIPTION: Each triangle is subdivided in 4 triangles.
%
%             /\1
%            /  \
%           /    \
%          /      \
%        4/--------\5
%        /  \    /  \
%       /    \6 /    \
%     2'--------------'3

%
% Authors: François Tadel, 2009-2016

% Parse inputs
if (nargin < 5) || isempty(isThresh)
    isThresh = 0;
end
if (nargin < 4) || isempty(threshEdge)
    threshEdge = [];
end
if (nargin < 3) || isempty(threshArea)
    threshArea = [];
end
% Check matrices orientation
if (size(Vertices, 2) ~= 3) || (size(Faces, 2) ~= 3)
    error('Faces and Vertices must have 3 columns (X,Y,Z).');
end

% ===== SPLIT EDGES =====
% Get all the edges of the surface
i1 = [Faces(:,1); Faces(:,1); Faces(:,2)];
i2 = [Faces(:,2); Faces(:,3); Faces(:,3)];
% List of edges to split in half
iSplit = [];
% Split long edges
if ~isempty(threshEdge) && (threshEdge > 0)
    % Compute the length of all the edges
    edgeLength = sqrt((Vertices(i1,1)-Vertices(i2,1)).^2 + (Vertices(i1,2)-Vertices(i2,2)).^2 + (Vertices(i1,3)-Vertices(i2,3)).^2);
    % Split long edges
    iSplit = [iSplit; find(edgeLength - mean(edgeLength) > threshEdge * std(edgeLength))];
end
% Split large faces
if ~isempty(threshArea) && (threshArea > 0)
    % Detect the faces that have an area above the threshold
    [tmp, iFacesSplit] = tess_threshold(Vertices, Faces, threshArea);
    % Split all the edges of the large surfaces
    iSplit = [iFacesSplit(:); iFacesSplit(:) + size(Faces,1); iFacesSplit(:) + 2*size(Faces,1)];
end
% Split all faces (if there are no other constraints)
if isempty(threshArea) && isempty(threshEdge)
    iSplit = (1:length(i1))';
end
% Nothing to split: return
if isempty(iSplit)
    return;
end

% ===== REFINE MESH =====
% New vertices
newVertices = unique((Vertices(i1(iSplit),:) + Vertices(i2(iSplit),:)) ./ 2, 'rows');
% Add to the existing vertices
Vertices = [Vertices; newVertices];

% ===== TESSELATE NEW SURFACE =====
% If 3D surface, use channel_tesselate.m
if ~all(Vertices(:,3) == Vertices(1,3))
    Faces = channel_tesselate(Vertices, isThresh);
    % Else: flat surface, use delaunay.m
else
    Faces = delaunay(Vertices(:,1), Vertices(:,2));
end

end

function [FaceArea, VertArea] = tess_area(Vertices, Faces)
% TESS_AREA: Compute the surface area associated with each face and each vertex.

% Compute the area of all the faces
r12 = Vertices(Faces(:,1),:);        % temporary holding
r13 = Vertices(Faces(:,3),:) - r12;  % negative of r31
r12 = Vertices(Faces(:,2),:) - r12;  % from 1 to 2
FaceArea = sqrt(sum(bst_cross(r12,r13,2).^2, 2)) / 2;

% Compute the triangle area only if needed
if (nargout >= 2)
    % Build vertex-face connectivity matrix, with the area information
    nFaces = size(Faces,1);
    rowno = double([Faces(:,1); Faces(:,2); Faces(:,3)]);
    colno = [1:nFaces, 1:nFaces, 1:nFaces]';
    data  = [FaceArea; FaceArea; FaceArea];
    VertFacesArea = sparse(rowno,colno,data);
    
    % Compute the vertex area: 1/3 of each triangle involving this vertex
    VertArea = 1/3 * full(sum(VertFacesArea,2));
end
end

function triPerimeter = tess_perimeter(Vertices, Faces)
% TESS_PERIMETER: Computes the perimeter of each face of the tesselation

%
% Authors: Francois Tadel, 2009-2011


my_norm = @(v)sqrt(sum(v .^ 2, 2));
% Get coordinates of Vertices for each face
vertFacesX = reshape(Vertices(reshape(Faces,1,[]), 1), size(Faces));
vertFacesY = reshape(Vertices(reshape(Faces,1,[]), 2), size(Faces));
vertFacesZ = reshape(Vertices(reshape(Faces,1,[]), 3), size(Faces));
% For each face : compute triangle perimeter
triSides = [my_norm([vertFacesX(:,1)-vertFacesX(:,2), vertFacesY(:,1)-vertFacesY(:,2), vertFacesZ(:,1)-vertFacesZ(:,2)]), ...
    my_norm([vertFacesX(:,1)-vertFacesX(:,3), vertFacesY(:,1)-vertFacesY(:,3), vertFacesZ(:,1)-vertFacesZ(:,3)]), ...
    my_norm([vertFacesX(:,2)-vertFacesX(:,3), vertFacesY(:,2)-vertFacesY(:,3), vertFacesZ(:,2)-vertFacesZ(:,3)])];
triPerimeter = sum(triSides, 2);

end

function Wmat = bst_shepards(destLoc, srcLoc, nbNeighbors, excludeParam, expDistance)
% BST_SHEPARDS: 3D nearest-neighbor interpolation using Shepard's weighting.
%
% USAGE:  Wmat = bst_shepards(destLoc, srcLoc, nbNeighbors=8, excludeParam=0, expDistance=2)
%
% INPUT:
%    - srcLoc       : Nx3 array of original locations, or tesselation structure (Faces,Vertices,VertConn)
%    - destLoc      : NNx3 array of locations onto original data will be interpolated, or tesselation structure (Faces,Vertices,VertConn)
%    - nbNeighbors  : Number of nearest neighbors to be considered in the interpolation (default is 8)
%    - excludeParam : If > 0, the source points that are two far away from the destination surface are ignored.
%                     Excluded points #i that have: (minDist(i) > mean(minDist) + excludeParam * std(minDist))
%                     where minDist represents the minimal distance between each source point and the destination surface
%                     If < 0, exclude the vertices that are further from the absolute distance excludeParam  (in millimeters)
%    - expDistance  : Distance exponent (if higher, influence of a value decreases faster)
%
% OUTPUT:
%    - Wmat : Interpolation matrix

%% ===== PARSE INPUTS =====
% Check number of arguments
if (nargin < 2)
    error('Usage: Wmat = bst_shepards(destLoc, srcLoc, nbNeighbors, excludeParam, expDistance)');
end
% Check matrices orientation
if ((size(destLoc, 2) ~= 3) || (size(srcLoc, 2) ~= 3)) && ((size(destLoc, 2) ~= 2) || (size(srcLoc, 2) ~= 2))
    error('destLoc and srcLoc must have 2 or 3 columns.');
end
% Argument: Number of neighbors for interpolation
if (nargin < 3) || isempty(nbNeighbors)
    nbNeighbors = 8;
end
% Argument: excludeParam
if (nargin < 4) || isempty(excludeParam)
    excludeParam = 0;
end
% Argument: expDistance
if (nargin < 5) || isempty(expDistance)
    expDistance = 2;
end

%% ===== SHEPARDS INTERPOLATION =====
% Allocate interpolation matrix
nDest = size(destLoc,1);
nSrc  = size(srcLoc,1);
% Maximum number of neighbors = number of electrodes
if (nbNeighbors > nSrc)
    nbNeighbors = nSrc;
end

% Find nearest neighbors
[I,dist] = bst_nearest(srcLoc, destLoc, nbNeighbors, 1);
% Square the distance matrix
dist = dist .^ 2;
% Eliminate zeros in distance matrix for stability
dist(dist == 0) = eps;

% One neighbor
if (nbNeighbors == 1)
    Wmat = sparse(1:nDest, I(:)', ones(1,nDest), nDest, nSrc);
    
    % More complicated cases
elseif (nbNeighbors > 1)
    % Interpolation weights from Shepards method
    W = (bst_bsxfun(@minus, dist(:,nbNeighbors), dist) ./ bst_bsxfun(@times, dist(:,nbNeighbors), dist)) .^ expDistance;
    sumW = sum(W(:,1:nbNeighbors-1),2);
    % Correct zero values: points overlap exactly => take only the first point
    iZeroW = find(sumW == 0);
    if ~isempty(iZeroW)
        sumW(iZeroW) = 1;
        W(iZeroW, 1:nbNeighbors-1) = ones(length(iZeroW),1) * [1,zeros(1,nbNeighbors-2)];
    end
    W = W(:,1:nbNeighbors-1) ./ (sumW * ones(1,nbNeighbors-1));
    % Create sparse matrix with those weights
    i = repmat((1:nDest)', nbNeighbors-1, 1);
    j = reshape(I(:, 1:nbNeighbors-1), [], 1);
    Wmat = sparse(i, j, W(:), nDest, nSrc);
end


%% ===== IGNORE VERTICES TOO FAR AWAY =====
% Set to zero the weights of the vertices that are too far away from the sources
% EEG: Distance relative to the mean distance between sensors
if (excludeParam > 0)
    % Find vertices that are too far from their nearest neighbors
    iTooFarVertices = (dist(:,1) > mean(dist(:,1)) + excludeParam * std(dist(:,1)));
    % Remove them from the interpolation matrix
    Wmat(iTooFarVertices, :) = 0;
    % SEEG/ECOG: Absolute distance
elseif (excludeParam < 0)
    % Find vertices that are too far from their nearest neighbors (in millimeters)
    iTooFarVertices = (sqrt(dist(:,1)) >  abs(excludeParam));
    % Remove them from the interpolation matrix
    Wmat(iTooFarVertices, :) = 0;
end

end

function cbars = CBar
cbars = [1 0 0;1 0.5 0;1 1 0;0.5 1 0.5;0 1 1;0 0.5 1;0.6 0.5 1];
cbars = lineInterp(cbars,256,8);
cbars = cbars(end:-1:1,:);
end