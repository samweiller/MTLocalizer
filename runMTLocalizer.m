function [MTLOC] = runMTLocalizer(sub, cbl, acq)
%% Start me up
clc
DOG.curDir = cd;
if isempty(sub); DOG.subID = input('\nPlease Enter Your Participant Code #: ', 's'); else DOG.subID = sub; end;
if isempty(cbl); DOG.run =  input('\nPlease Enter The Run #: ', 's'); else  DOG.run = cbl; end;
if isempty(acq); DOG.acq =  input('\nPlease Enter The Aquisition #: ', 's'); else  DOG.acq = acq; end;
PATH = fullfile(DOG.curDir, sprintf('DOG_S%d_R%d_A%d.mat', DOG.subID, DOG.run, DOG.acq));
save(PATH);
if ~exist(PATH);
    [Path, File] = uigetfile('*.mat', 'Select .MAT with RS');
    PATH = fullfile(Path, File);
end
load(PATH);

pause(.8)
disp('Dog Localizer')
pause(.7)
disp('  Version 0.50')
disp('  Mar. 20, 2014')
pause(.4)
disp('Script Written by Sam Weiller')
pause(3)
clc

%% Control Panel
designs = [...
    4 3 2 6 7 5 1 1 5 7 6 2 3 4;
    3 2 1 4 6 7 5 5 7 6 4 1 2 3;
    1 5 6 7 4 2 3 3 2 4 7 6 5 1;
    2 1 5 3 6 7 4 4 7 6 3 5 1 2;
    5 4 6 7 3 1 2 2 1 3 7 6 4 5;
    ];

numBlocks = size(designs, 2);



condition = designs(cbl, :);

%% MT Params
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
grey = ceil((white+black)/2);
ifi = Screen('GetFlipInterval', w);
[tw, th] = Screen('WindowSize', w);

% Scanning Parameters
initFixation = 32;            % in seconds
postFixation = initFixation;
% Session Parameters
blockDur = 20;               % should be a multiple of osc
% createBlockOrderMTloc;
% blockOrder = [2 3 2 3 2 3 2 3 2 3 2 3 2 3 2 3 2 3 2 3 2 3 2 3 2 3 2 3];

ibi = .5;                    % wait between blocks in seconds
% (very rough) setting of oscillation frequency
osc = .8;                    % Oscillating in seconds

nframes     = floor((blockDur-ibi)/ifi); % number of animation frames in loop
mon_width   = 39;   % horizontal dimension of viewable screen (cm)
v_dist      = 60;   % viewing distance (cm)
dot_speed   = 7;    % dot speed (deg/sec)
ndots       = 500; % number of dots
max_d       = 15;   % maximum radius of  annulus (degrees)
min_d       = 1;    % minumum
dot_w       = 0.6;  % width of dot (deg)
fix_r       = 0.15; % radius of fixation point (deg)
f_kill      = 0.05; % fraction of dots to kill each frame (limited lifetime)
differentcolors = 1; % Use a different color for each point if == 1. Use common color white if == 0.
differentsizes = 0; % Use different sizes for each point if >= 1. Use one common size if == 0.
waitframes = 1;     % Show new dot-images at each waitframes'th monitor refresh.

ppd = pi * (winRect(3)-winRect(1)) / atan(mon_width/v_dist/2) / 360;    % pixels per degree
pfs = dot_speed * ppd / (1/ifi);                            % dot speed (pixels/frame)
s = dot_w * ppd;                                        % dot size (pixels)
fix_cord = [[tw/2 th/2]-fix_r*ppd [tw/2 th/2]+fix_r*ppd];
rmax = max_d * ppd;	% maximum radius of annulus (pixels from center)
rmin = min_d * ppd; % minimum
r = rmax * sqrt(rand(ndots,1));	% r
r(r<rmin) = rmin;
t = 2*pi*rand(ndots,1);                     % theta polar coordinate
cs = [cos(t), sin(t)];
xy = [r r] .* cs;   % dot positions in Cartesian coordinates (pixels from center)

mdir = ones(ndots,1);%2 * floor(rand(ndots,1)+0.5) - 1;    % motion direction (in or out) for each dot
dr = pfs * mdir;                            % change in radius per frame (pixels)
dxdy = [dr dr] .* cs;                       % change in x and y per frame (pixels)

% Create a vector with different colors for each single dot, if
% requested:
if (differentcolors==1)
    % colvect = uint8(round(rand(3,ndots)*255));
    colvect = grey + (grey-1)* (2 * floor(rand(ndots,1)+0.5) - 1);
    colvect = [colvect colvect colvect]';
else
    colvect=white;
end;

% Create a vector with different point sizes for each single dot, if
% requested:
if (differentsizes>0)
    s=(1+rand(1, ndots)*(differentsizes-1))*s;
end;

%% Main Loop
experimentStartTime = GetSecs;
for blocks = 1:numBlocks
    fixStartLog = GetSecs;
    tempCov = max(max(designs)) + 1; % Defines a non-zero covariate number for fixation covariate files.
    
    timeLogger.block(blocks).startTime = GetSecs - experimentStartTime
    TIMER.block(blocks).start = GetSecs - TIMER.init;
    
    fixate(w);
    while GetSecs <= initFixation
        % Wait for Initial Fixation
    end;
    
    %         WaitSecs(fixationTime);
    TIMER.block(blocks).end = GetSecs - TIMER.init;
    TIMER.block(blocks).length = TIMER.block(blocks).end - TIMER.block(blocks).start;
    TIME_MAT(blocks, 3) = (GetSecs - TIMER.init) - TIME_MAT(blocks, 2);
    TIME_MAT(blocks, 4) = 1;
    if condition(blocks) == 0 % fixation
        
    elseif condition(blocks) == 6
        motionStartLog = GetSecs;
        motionEnd = motionStartLog + motionLength;
        vbl = Screen('Flip', w, 0, 1);
        % --------------
        % animation loop
        % --------------
        for i = 1:nframes
            if GetSecs <= motionEnd
                if i>1
                    Screen('FillOval', w, uint8(white), fix_cord);	% draw fixation dot (flip erases it)
                    Screen('DrawDots', w, xymatrix, s, colvect, [tw/2 th/2],1);  % change 1 to 0 to draw square dots
                    Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
                end
                
                mdir = ones(ndots,1) * mod(ceil(i/(nframes/4)),2)*2-1;
                dr = pfs * mdir;                            % change in radius per frame (pixels)
                dxdy = [dr dr] .* cs;                       % change in x and y per frame (pixels)
                
                xy = xy + dxdy;						% move dots
                r = r + dr;							% update polar coordinates too
                
                r_out = find(r > rmax | r < rmin | rand(ndots,1) < f_kill);	% dots to reposition
                nout = length(r_out);
                
                if nout
                    r(r_out) = rmax * sqrt(rand(nout,1));
                    r(r<rmin) = rmin;
                    t(r_out) = 2*pi*(rand(nout,1));
                    
                    cs(r_out,:) = [cos(t(r_out)), sin(t(r_out))];
                    xy(r_out,:) = [r(r_out) r(r_out)] .* cs(r_out,:);
                    
                    dxdy(r_out,:) = [dr(r_out) dr(r_out)] .* cs(r_out,:);
                end;
                xymatrix = transpose(xy);
                
                vbl=Screen('Flip', w, vbl + (waitframes-0.5)*ifi);
            else
                break;
            end;
        end;
        
        disp('Time:');
        disp(GetSecs - motionStartLog);
        TIMER.block(blocks).end = GetSecs - TIMER.init;
        TIMER.block(blocks).length = TIMER.block(blocks).end - TIMER.block(blocks).start;
        TIME_MAT(blocks, 3) = (GetSecs - TIMER.init) - TIME_MAT(blocks, 2);
        TIME_MAT(blocks, 4) = 1;
    elseif condition(blocks) == 7
        motionStartLog = GetSecs;
        motionEnd = motionStartLog + motionLength;
        vbl = Screen('Flip', w, 0, 1);
        for i = 1:nframes
            if GetSecs <= motionEnd
                if i>1 && rem(i,10)~=0
                    Screen('FillOval', w, uint8(white), fix_cord);	% draw fixation dot (flip erases it)
                    Screen('DrawDots', w, xymatrix, s, colvect, [tw/2 th/2],1);  % change 1 to 0 to draw square dots
                    Screen('DrawingFinished', w); % Tell PTB that no further drawing commands will follow before Screen('Flip')
                end
                
                r_out = find(r > rmax | r < rmin | rand(ndots,1) < f_kill/100);	% dots to reposition
                nout = length(r_out);
                
                if nout && ~mod(i/(nframes/4),2)
                    r(r_out) = rmax * sqrt(rand(nout,1));
                    r(r<rmin) = rmin;
                    t(r_out) = 2*pi*(rand(nout,1));
                    
                    cs(r_out,:) = [cos(t(r_out)), sin(t(r_out))];
                    xy(r_out,:) = [r(r_out) r(r_out)] .* cs(r_out,:);
                    
                    dxdy(r_out,:) = [dr(r_out) dr(r_out)] .* cs(r_out,:);
                end;
                xymatrix = transpose(xy);
                
                vbl=Screen('Flip', w, vbl + (waitframes-0.5)*ifi);
            else
                break;
            end;
        end
    end;
    %         disp('Time:');
    %         disp(GetSecs - motionStartLog);
    TIMER.block(blocks).end = GetSecs - TIMER.init;
    TIMER.block(blocks).length = TIMER.block(blocks).end - TIMER.block(blocks).start;
    TIME_MAT(blocks, 3) = (GetSecs - TIMER.init) - TIME_MAT(blocks, 2);
    TIME_MAT(blocks, 4) = 1;
end;

function [w rect xc yc] = startPTB(screenNumber, oGl, color)

if nargin == 0
    oGl = 0;
    color = [0 0 0];
elseif nargin == 1;
    color = [0 0 0];
end;

[w rect] = Screen('OpenWindow', screenNumber, color);
xc = rect(3)/2;
yc = rect(4)/2;

if oGl == 1
    AssertOpenGL;
    Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA, [1 1 1 1]);
end;

function fixate(w)
Screen('TextSize', w, 40);
DrawFormattedText(w, '+', 'center', 'center', [200 200 200]);
Screen('TextSize', w, 25);
Screen('Flip', w);

function flipNwait(w)
Screen('Flip', w);
WaitSecs(.3);
KbCheck(-1);

function trigger(triggerKey)

KbName('UnifyKeyNames');

touch = 0;
go = 0;
while go == 0
    [touch, secs, keyCode] = KbCheck(-1);
    WaitSecs(.0001);
    if touch && keyCode(triggerKey)
        go = 1;
    end;
end;