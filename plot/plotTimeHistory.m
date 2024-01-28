function plotTimeHistory(model, struData, reducedBasis, modes, tout, options)

% Set default options
iOpt = 0;
iOpt = iOpt+1; baseOpt.scale = 10;        descr{iOpt} = 'Scale factor to magnify displacements [1].';
iOpt = iOpt+1; baseOpt.dt = 0.001;        descr{iOpt} = 'Step size in seconds for the output calculation [0.001].';
iOpt = iOpt+1; baseOpt.step = 0.001;      descr{iOpt} = 'Step size in seconds between one frame and the other [0.01].';
iOpt = iOpt+1; baseOpt.xlim = [-0.2 1];   descr{iOpt} = 'xlim for the plot [-0.2 1].';
iOpt = iOpt+1; baseOpt.ylim = [-0.6 0.6]; descr{iOpt} = 'xlim for the plot [-0.6 0.6].';
iOpt = iOpt+1; baseOpt.zlim = [-0.2 0.5]; descr{iOpt} = 'ylim for the plot [-0.2 0.5].';

if nargin==0
    printOptionDescription(baseOpt, descr);
    return
end

% Process input options
options = setOptions(baseOpt, 'ignore', options);

%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

% First, we need to interpolate to have the variables at constant time
% steps
timeVect = seconds(tout-tout(1));
TT1 = timetable(timeVect.',modes.');
TT2 = retime(TT1,'regular','linear','TimeStep',seconds(options.dt));

Tgz = struData.Tgz;
U = reducedBasis.V;

[~, gridDOF, spointDOF] = getIDtable(model.Node.ID, model.Spoint.ID);

% Define parameters
outputVideoFile = 'outputVideo.mp4';
frameRate = 30; % Adjust the frame rate as needed

% Create VideoWriter object
outputVideo = VideoWriter(outputVideoFile, 'MPEG-4');
outputVideo.FrameRate = frameRate;

% Open the VideoWriter object
open(outputVideo);

step = options.step/options.dt;
for time = 25000:step:30000%time = 1:step:size(TT2,1)
    deformedModel = model;
    [MODES, ~] = getModeShapes(gridDOF, spointDOF, Tgz, U*TT2.Var1(time,:).');
    for i = 1:size(model.Node.Coord,2)
        deformedModel.Node.Coord(:,i) = model.Node.Coord(:,i)+MODES(1:3,i)*options.scale;
    end
    figHandle = plotNeoModel(deformedModel, [], 0, false);
    title ''
    axis off
    xlim(options.xlim)
    ylim(options.ylim)
    zlim(options.zlim)
    view([74, 24])
    set(gcf,'color','w')
    frame = getframe(figHandle); 
    writeVideo(outputVideo, frame);
    close(figHandle)
end

close(outputVideo);
end
