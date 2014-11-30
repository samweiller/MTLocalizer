function [MTLOC] = runMTLocalizer_long(sub, cbl, acq)
%% Start me up
clc
MTLOC.curDir = cd;
if isempty(sub); MTLOC.subID = input('\nPlease Enter Your Participant Code #: ', 's'); else MTLOC.subID = sub; end;
if isempty(cbl); MTLOC.run =  input('\nPlease Enter The Run #: ', 's'); else  MTLOC.run = cbl; end;
if isempty(acq); MTLOC.acq =  input('\nPlease Enter The Aquisition #: ', 's'); else  MTLOC.acq = acq; end;
PATH = fullfile(MTLOC.curDir, sprintf('MTLOC_S%d_C%d_A%d.mat', MTLOC.subID, MTLOC.run, MTLOC.acq));
save(PATH);
if ~exist(PATH);
    [Path, File] = uigetfile('*.mat', 'Select .MAT with RS');
    PATH = fullfile(Path, File);
end
load(PATH);

pause(.8)
disp('MT Localizer')
pause(.7)
disp('  Version 1.00')
disp('  Nov. 30, 2014')
pause(.4)
disp('Script Written by Sam Weiller')
pause(3)
clc

%% Control Panel
designs = [...
    3 1 2 1 2 1 2 1 2 1 2 1 2 1 2 3;
    3 2 1 2 1 2 1 2 1 2 1 2 1 2 1 3;
    ];

triggerKey = KbName('t');

numBlocks = size(designs, 2);
condition = designs(cbl, :);

%% PTB Setup
screenNumber = 0;
Screen('Preference', 'SkipSyncTests', 2);
[w winRect xMid yMid] = startPTB(screenNumber, 1, [128 128 128]);
HideCursor;

%% MT Params
white = WhiteIndex(screenNumber);
black = BlackIndex(screenNumber);
grey = ceil((white+black)/2);
ifi = Screen('GetFlipInterval', w);
[tw, th] = Screen('WindowSize', w);

% Scanning Parameters
initFixation = 32;            % in seconds

% Session Parameters
blockDur = 16;               % should be a multiple of osc

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
Screen('TextSize', w, 20);
DrawFormattedText(w, 'Waiting for trigger...', 'center', 'center', 0);
Screen('Flip', w);
trigger(triggerKey);

experimentStartTime = GetSecs;
for blocks = 1:numBlocks
    timeLogger.block(blocks).startTime = GetSecs - experimentStartTime;
    timeLogger.block(blocks).condition = condition(blocks);
    
    if condition(blocks) == 3 % fixation
        fixationEndTime = GetSecs + blockDur;
        
        fixate(w);
        while GetSecs <= fixationEndTime
            % Wait for Initial Fixation
        end;
    elseif condition(blocks) == 1
        motionStartLog = GetSecs;
        motionEnd = motionStartLog + blockDur;
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
    elseif condition(blocks) == 2
        motionStartLog = GetSecs;
        motionEnd = motionStartLog + blockDur;
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
    timeLogger.block(blocks).endTime = GetSecs - experimentStartTime;
    timeLogger.block(blocks).length  = timeLogger.block(blocks).endTime - timeLogger.block(blocks).startTime;
end;

%% Logging & Cleanup
MTLOC.ANSMAT = ANSMAT;
save(PATH, 'MTLOC', 'timeLogger');

cov1Filename = sprintf('MTLOC%02d_CBL%02d_Acq%02d_Cov1_motion.txt', sub, cbl, acq);
cov2Filename = sprintf('MTLOC%02d_CBL%02d_Acq%02d_Cov2_flicker.txt', sub, cbl, acq);
cov3Filename = sprintf('MTLOC%02d_CBL%02d_Acq%02d_Cov3_fixation.txt', sub, cbl, acq);

for block = 1:numBlocks
    temp = [round(timeLogger.block(block).startTime), round(timeLogger.block(block).endTime), round(timeLogger.block(block).length)];
    
    eval(sprintf('dlmwrite(cov%dFilename, temp, ''delimiter'', ''\t'', ''-append'');', timeLogger.block(block).condition));
end;

%% Shutdown Procedures
ShowCursor;
clear screen;

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