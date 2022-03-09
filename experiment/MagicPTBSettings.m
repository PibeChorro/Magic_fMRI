function ptb = MagicPTBSettings()
% Turn initial start screen from white to black
Screen('Preference', 'VisualDebugLevel', 1);

% Decide which set up to use
SetUp = 2;

% Here we call some default settings for setting up Psychtoolbox
PsychDefaultSetup(2);
% ListenChar(2);

% Keyboard pre-setup
KbName('UnifyKeyNames');

% Set a new random seed.
ptb.seed = rng('shuffle');

% Get the screen numbers. This gives us a number for each of the screens
% attached to our computer.
ptb.screens = Screen('Screens');

% To draw we select the maximum of these numbers. So in a situation where we
% have two screens attached to our monitor we will draw to the external
% screen.
ptb.screenNumber = max(ptb.screens);

% Define black and white (white will be 1 and black 0). This is because
% in general luminace values are defined between 0 and 1 with 255 steps in
% between. All values in Psychtoolbox are defined between 0 and 1
ptb.white = WhiteIndex(ptb.screenNumber);
ptb.black = BlackIndex(ptb.screenNumber);

% Do a simply calculation to calculate the luminance value for grey. This
% will be half     the luminace values for white
ptb.grey = ptb.white / 2;

% Open an on screen window and color it grey/or RGB. This function returns a
% number that identifies the window we have opened "window" and a vector
% "windowRect".
% "windowRect" is a vector of numbers: the first is the X coordinate
% representing the far left of our screen, the second the Y coordinate
% representing the top of our screen,
% the third the X coordinate representing
% the far right of our screen and finally the Y coordinate representing the
% bottom of our screen.

ptb.KeyList1 = zeros (256, 1); % initiate an empty array that later, when keys are defined, is filled with the allowed keys
ptb.KeyList2 = zeros (256, 1); % initiate an empty array that later, when keys are defined, is filled with the allowed keys
ptb.usbTrg = 1; % If 1 --> wait for scanner triggers & check USB inputs.


if SetUp == 0
    ptb.FontColor = [1 1 1];
    ptb.BackgroundColor = ptb.grey;
    [ptb.window, ptb.windowRect] = PsychImaging('OpenWindow', ptb.screenNumber, ptb.BackgroundColor, [0 0 1000 800]); % For troubleshooting
    ptb.widthMonitor = 153;
    ptb.heightMonitor = 123;
    %                                       The KeyList must be filled with doubles
    ptb.Keys.escape     = KbName ('ESCAPE');ptb.KeyList1(ptb.Keys.escape)= double(1);
    ptb.Keys.yes        = KbName ('y');     ptb.KeyList2(ptb.Keys.yes)   = double(1);
    ptb.Keys.no         = KbName ('n');     ptb.KeyList2(ptb.Keys.no)    = double(1);
    ptb.Keys.one        = KbName ('1!');    ptb.KeyList2(ptb.Keys.one)  = double(1);
    ptb.Keys.two        = KbName ('2"');    ptb.KeyList2(ptb.Keys.two)   = double(1);
    ptb.Keys.three      = KbName ('3&');    ptb.KeyList2(ptb.Keys.three) = double(1);
    ptb.Keys.four       = KbName ('4$');    ptb.KeyList2(ptb.Keys.four)  = double(1);
    ptb.Keys.five       = KbName ('5%');    ptb.KeyList2(ptb.Keys.five)   = double(1);
    ptb.Keys.debug      = 1;
elseif SetUp == 1
    ptb.FontColor = [0 1 1];
    ptb.BackgroundColor = ptb.black;
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'UsePanelFitter', [1600 1080], 'Aspect'); 
    [ptb.window, ptb.windowRect] = PsychImaging('OpenWindow', ptb.screenNumber, ptb.BackgroundColor);
    
    ptb.FontSize = Screen('TextSize', ptb.window, 30);
    ptb.widthMonitor = 400;
    ptb.heightMonitor = 350;
    %                                       The KeyList must be filled with doubles
    ptb.Keys.escape     = KbName ('ESCAPE');ptb.KeyList1(ptb.Keys.escape)= double(1);
    ptb.Keys.yes        = KbName ('/');     ptb.KeyList2(ptb.Keys.yes)   = double(1);
    ptb.Keys.no         = KbName ('*');     ptb.KeyList2(ptb.Keys.no)    = double(1);
    ptb.Keys.one        = KbName ('1');     ptb.KeyList2(ptb.Keys.one)   = double(1);
    ptb.Keys.two        = KbName ('2');     ptb.KeyList2(ptb.Keys.two)   = double(1);
    ptb.Keys.three      = KbName ('3');     ptb.KeyList2(ptb.Keys.three) = double(1);
    ptb.Keys.four       = KbName ('4');     ptb.KeyList2(ptb.Keys.four)  = double(1);
    ptb.Keys.five       = KbName ('5');     ptb.KeyList2(ptb.Keys.five)  = double(1);
    ptb.Keys.debug      = 1;
elseif SetUp == 2
    ptb.FontColor = [0 1 1];
    ptb.BackgroundColor = ptb.black;
    PsychImaging('PrepareConfiguration');
    PsychImaging('AddTask', 'General', 'UsePanelFitter', [1600 1080], 'Aspect');
    [ptb.window, ptb.windowRect] = PsychImaging('OpenWindow', ptb.screenNumber, ptb.BackgroundColor); 
    
    ptb.FontSize = Screen('TextSize', ptb.window, 40);
    ptb.widthMonitor = 399;
    ptb.heightMonitor = 224;
    % Because of the MR compatible keyboard we flip the order of button
    % presses. Additionally the buttons for the binary answers are index
    % (button 4$) and middle finger (button 3#)
    %                                       The KeyList must be filled with doubles
    ptb.Keys.escape     = KbName ('ESCAPE');ptb.KeyList1(ptb.Keys.escape)= double(1);
    ptb.Keys.yes        = KbName ('y');     ptb.KeyList1(ptb.Keys.yes)   = double(1);
    ptb.Keys.no         = KbName ('n');     ptb.KeyList1(ptb.Keys.no)    = double(1);
    ptb.Keys.yes_sub    = KbName ('0)');   %ptb.Keys.yes_sub    = KbName ('9(');
    ptb.Keys.no_sub     = KbName ('8*');
    ptb.Keys.one        = KbName ('6^');    ptb.KeyList2(ptb.Keys.one)   = double(1);
    ptb.Keys.two        = KbName ('7&');    ptb.KeyList2(ptb.Keys.two)   = double(1);
    ptb.Keys.three      = KbName ('8*');    ptb.KeyList2(ptb.Keys.three) = double(1);
    ptb.Keys.four       = KbName ('9(');    ptb.KeyList2(ptb.Keys.four)  = double(1);
    ptb.Keys.five       = KbName ('0)');    ptb.KeyList2(ptb.Keys.five)  = double(1);
    ptb.Keys.trg        = KbName ('w');     ptb.KeyList2(ptb.Keys.trg)   = double(1); % The scanner sends 'w' as USB keyboard input (from keyboard 2)
    ptb.Keys.debug      = 1;
end


%............................. KEYS ......................................%
% Get Keyboard indices
[keyboardIndices, productNames, ~] = GetKeyboardIndices();

if SetUp == 0 % Own computer
    ptb.Keys.kbrd2  = max(GetKeyboardIndices);
    ptb.Keys.kbrd1  = min(GetKeyboardIndices);
    fprintf('\n=> Experimentators keyboard Nr.: %u  %s \n',ptb.Keys.kbrd1, productNames{1});
    fprintf('\n=> Subjects keyboard Nr.: %u  %s \n',ptb.Keys.kbrd2, productNames{end});
elseif SetUp == 1 % CIN computer % PRG please find out which one is the right one!
    ptb.Keys.kbrd1 = -1;
    ptb.Keys.kbrd2 = -1;
%     for devi = 1:length(keyboardIndices)
%         if strcmp(productNames(devi), 'HID 04d9:1203') && keyboardIndices(devi) == 11 %Because there are 2 with the same name.
%             ptb.Keys.kbrd2  = keyboardIndices(devi);
%             fprintf('\n=> Subjects keyboard Nr.: %u  %s \n',ptb.Keys.kbrd2, productNames{devi});
%         elseif strcmp(productNames(devi), 'HID 046a:0023') && keyboardIndices(devi) == 10 %Because there are 2 with the same name. 
%             ptb.Keys.kbrd1  = keyboardIndices(devi);
%             fprintf('\n=> Experimentators keyboard Nr.: %u  %s \n',ptb.Keys.kbrd1, productNames{devi});
%         end
%     end
elseif SetUp == 2 % MPI computer!
    for devi = 1:length(keyboardIndices)
        if strcmp(productNames(devi), 'P.I. Engineering Xkeys')
            ptb.Keys.kbrd2  = keyboardIndices(devi);
            fprintf('\n=> Subjects keyboard Nr.: %u  %s \n',ptb.Keys.kbrd2, productNames{devi});
        elseif strcmp(productNames(devi), 'Dell Dell USB Entry Keyboard')
            ptb.Keys.kbrd1  = keyboardIndices(devi);
            fprintf('\n=> Experimentators keyboard Nr.: %u  %s \n',ptb.Keys.kbrd1, productNames{devi});
        end
    end
end

if ptb.Keys.kbrd1 == ptb.Keys.kbrd2
   warning('You are using only one keyboard!'); 
end

% create & start KbQueue for the first Keyboard.
KbQueueCreate(ptb.Keys.kbrd1, ptb.KeyList1);
KbQueueStart(ptb.Keys.kbrd1);

% create & start KbQueue for the second Keyboard.
KbQueueCreate(ptb.Keys.kbrd2, ptb.KeyList2);
KbQueueStart(ptb.Keys.kbrd2);


%.........................................................................%
% Query the inter-frame-interval. This refers to the minimum possible time
% between drawing to the screen
ptb.ifi = Screen('GetFlipInterval', ptb.window);


ptb.ResponseTime = 2 - ptb.ifi/2;    % After video presentation the sub has 2 seconds to press a botton to answer the question
ptb.VideoLength  = 14 - ptb.ifi/2;   % Every video should have a length of 14 seconds / after 14 seconds the question should be presented

% This function call will give use the same information as contained in
% "windowRect"
ptb.rect = Screen('Rect', ptb.window);

% Get the size of the on screen window in pixels, these are the last two
% numbers in "windowRect" and "rect"
[ptb.screenXpixels, ptb.screenYpixels] = Screen('WindowSize', ptb.window);

% Get the centre coordinate of the window in pixels.
% xCenter = screenXpixels / 2
% yCenter = screenYpixels / 2
[ptb.xCenter, ptb.yCenter] = RectCenter(ptb.windowRect);

% Retreive the maximum priority number
ptb.topPriorityLevel = MaxPriority(ptb.window);

% Length of time and number of frames we will use for each drawing test
ptb.numSecs = 1;
ptb.numFrames = round(ptb.numSecs / ptb.ifi);

% Numer of frames to wait when specifying good timing. Note: the use of
% wait frames is to show a generalisable coding. For example, by using
% waitframes = 2 one would flip on every other frame. See the PTB
% documentation for details. In what follows we flip every frame.
ptb.waitframes = 1;

% We can also determine the refresh rate of our screen. The
% relationship between the two is: ifi = 1 / hertz
ptb.hertz = FrameRate(ptb.window);

% We can also query the "nominal" refresh rate of our screen. This is
% the refresh rate as reported by the video card. This is rounded to the
% nearest integer. In reality there can be small differences between
% "hertz" and "nominalHertz"
% This is nothing to worry about. See Screen FrameRate? and Screen
% GetFlipInterval? for more information
ptb.nominalHertz = Screen('NominalFrameRate', ptb.window);

% Here we get the pixel size. This is not the physical size of the pixels
% but the color depth of the pixel in bits
ptb.pixelSize = Screen('PixelSize', ptb.window);

% Queries the display size in mm as reported by the operating system. Note
% that there are some complexities here. See Screen DisplaySize? for
% information. So always measure your screen size directly.
[ptb.width, ptb.height] = Screen('DisplaySize', ptb.screenNumber);

% Get the maximum coded luminance level (this should be 1)
ptb.maxLum = Screen('ColorRange', ptb.window);

% Set up alpha-blending for smooth (anti-aliased) lines
Screen('BlendFunction', ptb.window, 'GL_SRC_ALPHA', 'GL_ONE_MINUS_SRC_ALPHA');          

% Here we set the size of the arms of our fixation cross
% PRG somethings missing.

% Real world variable
ptb.DistToMonitor = 1050;   %Distance to monitor in mm
% Determine the number of horizontal and vertical pixels of your screen window 
[ptb.widthWindow, ptb.heightWindow]=Screen('WindowSize', ptb.window);

% Determine the width and hight of your monitor (Use caution as these values are not always correct!)
% [widthMonitor, heightMonitor]=Screen('DisplaySize', window); incorrect

ptb.pixWidth =  ptb.widthMonitor/ptb.widthWindow; % width of single pixel in mm 
ptb.pixHeight = ptb.heightMonitor/ptb.heightWindow; % height of single pixel in mm

ptb.DegPerPixWidth = 2*atand((0.5*ptb.pixWidth)/ptb.DistToMonitor);
ptb.PixPerDegWidth = 1/ptb.DegPerPixWidth;
ptb.DegPerPixHeight = 2*atand((0.5*ptb.pixHeight)/ptb.DistToMonitor);
ptb.PixPerDegHeight = 1/ptb.DegPerPixHeight;
end


% Alternative way of restricting which keys you actually check (without
% giving a KeyList
% KbCheckList = [ptb.Keys.escape, ptb.Keys.yes,ptb.Keys.no,ptb.Keys.one, ...
%                ptb.Keys.two, ptb.Keys.three, ptb.Keys.four, ptb.Keys.five, ptb.Keys.trg];
% DisableKeysForKbCheck(ptb.KeyList);




