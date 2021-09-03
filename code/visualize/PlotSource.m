function handles = PlotSource(SourceData,BrainModel,varargin)
%PlotSource plot cortex domain activity
% handles = PlotSource(SourceData,BrainModel,varargin)
%   Input:
%      -Values(needed): []/ nchan X 1 vector. If Values is
%               empty([]), then PlotSource plots EEG channel locations; if Values
%               is nchan by 1 vector(where nchan is the number of channels,length(chanlocs)), 
%               then topo3d plots spatial weights topography;
%      -BrainModel(needed but the programme load it by default): 
%               Settings of cortex
%   Options:
%      -CLim(default [0,1]):color limits of topography, 2 by 1 vector
%      -roiIdx (default is 1:length(vertices)): region of interst to plot,vector 
%      -cortexcolor(default is [0.75 0.75 0.75]/[0.9 0.7 0.7]): color of cortex; 3 by
%                1 vector.

%      -colormap:colormap, 
%      -Axes: default is gca; plot topography on specific axis if provided.

if nargin<3||isempty(BrainModel)
    try
        path2file = which('PlotSource');
        path2set = fileparts(path2file);
        path2set = fullfile(path2set,'BrainModel');
        load(path2set,'BrainModel');
    catch
        error('No brain model found');
    end
end
if nargin<1
    SourceData = [];
end

loc=BrainModel.Vertices;tri=BrainModel.Faces;

% default setting
smooth = [];
thresh = 0;
roiIdx = 1:length(loc);
hemisphere = [];
ax2plot = gca;
CLim = [0,1];
colorMap = CBar;
cortexcolor=[.9 .7 .7];

nargs = nargin;
if nargs > 3
    for i = 1:2:length(varargin)
        Param = lower(varargin{i});
        Value = varargin{i+1};
        switch Param
            case 'smooth'
                smooth = Value;
                if Value<=0 || Value>100
                   warning('smooth must be in (0,100]');
                   smooth = 30;
                end
            case 'cortexcolor'
                cortexcolor = Value;
                if size(Value,2)~=3
                    error('Head color must be a 1 x 3 matrix')
                end
            case 'colormap'
                colorMap = Value;
                if size(Value,2)~=3
                    error('Colormap must be a n x 3 matrix')
                end
            case 'clim'
                CLim = Value;
                if numel(CLim)~=2
                    error('Color Limit must be 1 x 2 vector');
                end
            case 'roiidx'
                roiIdx = Value;
            case 'axes'
                ax2plot = Value;
        end
    end
end


if ~isempty(smooth)
    SmoothValue=smooth/100;
    VertConn=BrainModel.VertConn;
    if ~isempty(SmoothValue)
        iVertices=1:size(loc,1);
        % Smoothing factor
        SurfSmoothIterations = ceil(300 * SmoothValue * length(iVertices) / 100000);
        % Calculate smoothed vertices locations
        loc_sm=loc;
        loc_sm(iVertices,:) = tess_smooth(loc(iVertices,:), SmoothValue, SurfSmoothIterations, VertConn(iVertices,iVertices), 1);
        % Apply smoothed locations
        loc=loc_sm;
        
    end
end

if isempty(ax2plot)
    handles.h = figure('color',[1 1 1]);
    handles.axes = axes('Parent',handles.h);
else
    handles.axes = ax2plot;
end

handles.hp=patch(ax2plot,'vertices',loc,'faces',tri,...
    'FaceColor',cortexcolor,'edgecolor', 'none','facelighting','gouraud'...
    ,'specularstrength',0.2,'ambientstrength',0.5,'diffusestrength',0.5,...
    'BackfaceLighting','lit','AmbientStrength',0.5,'SpecularExponent',1,...
    'SpecularColorReflectance',0.5,'EdgeLighting','gouraud');
material dull
% camlight('headlight','infinite');
set(handles.axes,'Xcolor',[1 1 1],'ycolor',[1 1 1],'zcolor',[1 1 1],'CameraPosition', [0 100 0],'CameraViewAngle',6);
view([ -90 90 ])
axis equal
axlim=[min(loc);max(loc)]*1.2;
axis(reshape(axlim,1,[]));
hr = rotate3d;
LightSource = camlight('headlight');set(LightSource,'style','infinite');
set(hr,'Enable','on','ActionPostCallback',{@move_light_source,LightSource});


if ~isempty(SourceData)
    cortexcolor=[.75 .75 .75];
    SourceData = SourceData(roiIdx);
    cdata=repmat(cortexcolor,length(SourceData),1);
    
    if ~isempty(hemisphere)
        Valpha = ones(size(loc,1),1);
    else
        switch hemisphere
            case 'left'
                hemiIdx = BrainModel.Atlas(3).Scouts(1).Vertices;
                Valpha = zeros(size(loc,1),1);
                Valpha(hemiIdx) = 1;
                roiIdx = intersect(roiIdx,hemiIdx);
            case 'right'
                hemiIdx = BrainModel.Atlas(3).Scouts(2).Vertices;
                Valpha = zeros(size(loc,1),1);
                Valpha(hemiIdx) = 1;
                roiIdx = intersect(roiIdx,hemiIdx);
        end
    end
    
    s2plot = intersect(find(SourceData>thresh),roiIdx);
    cmin = min(CLim);cmax = max(CLim);
    cdata(s2plot,:)=colorMap(floor(min((SourceData(s2plot)-cmin)/(cmax-cmin),1)*(length(colorMap)-1))+1,:);
    
    
    set(handles.hp,'FaceAlpha','flat','AlphaDataMapping','none','FaceVertexAlphaData',Valpha, 'FaceVertexCData',cdata, 'facecolor', 'interp');
    caxis([cmin,cmax]);
end




end

function move_light_source(src,evt,LightSource)
LightSource = camlight(LightSource,'headlight');
end

function cbars = CBar
cbars = [1 0.3 0;1 1 0;1 1 1];
cbars = lineInterp(cbars,256,8);
cbars = cbars(end:-1:1,:);
end

function Cbar = lineInterp(cbar,cnum,pernum)
np = size(cbar,1);
nump = max(round(cnum/(np-1)),pernum);
for ic = 1:3
    tcolor = cbar(:,ic);
    ncolor = [];
    for itc = 1:np-1
        ncolor = [ncolor,linspace(tcolor(itc),tcolor(itc+1),nump)];
    end
    Cbar(:,ic) = ncolor';
end
while size(Cbar,1)>cnum
    Cbar(round(end/2),:)=[];
end
end