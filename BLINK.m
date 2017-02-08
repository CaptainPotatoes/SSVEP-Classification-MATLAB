% SSVEP PTB
% James Norton
% University of Illinois
% This script depends on PsychToolBox, which must be obtained separately.
% v1.1 - 01/16/16
clear all
close all
clc

%% ---------------------------
% SETTINGS 
freqsToUse = [5 10 15 20 30] %& harmonics;
freqs = [7 10]; % what frequecies, should be frequencies that the monitor can achieve
iSi = 2; %inter stimulus interval in  seconds
fullscreen = false; % do you want this to be full screen
imageSizeW = 800; % ssvep stim size, must be an even number
imageSizeH = 800; % ssvep stim size, must be an even number   
startTrial = 1; % NOT IMPLEMENTED - SHOULD ALWAYS BE 1...future - you can load a trial order, and start from the middle, in case the program crashes 
numTrials = 1; % how many trials do you want?   
trial_length = 10; % how long do you want each trial to be? recommend at least 5 seconds
triggerSize = [200 0; 0 0; 0 200]; % we trigger using a photodiode...you can just set this to all 0s 

%% ---------------------------
% START PSYCHTOOLBOX
sca;
Screen('Preference', 'SkipSyncTests', 0);



%% ------------------------
% Open the screen, and set up some constants
%  ------------------------
% HideCursor; 
AssertOpenGL; 
screenNumber = max(Screen('Screens')); % Get the maximum screen number i.e. get an external screen if avaliable
if fullscreen
    [wPtr rect] = Screen('OpenWindow',screenNumber,200,[],32,2);  %opens a window,  monitor, with black background, with dimensions
else
    [wPtr rect] = Screen('OpenWindow',screenNumber,200,[20 20 1920 950],32,2);  %opens a window,  monitor, with black background, with dimensions
end
Priority(MaxPriority(wPtr));
[center(1) center(2)] = WindowCenter(wPtr); %find the fixation
font_size = 24;

Distance = -900;
Screen('TextSize',wPtr,font_size); %Set the font size

%% ------------------------
% Set the frequency of stimulation
%--------- ---------------
refresh = Screen('GetFlipInterval',wPtr);  %refresh period in seconds
slack = refresh /2; %compute slack to keep timing

%% ------------------------
% Create the stimuli 
%------------------------
black = [0; 0; 0];
white = [255; 255; 255];

%% ------------------------
% Display the flashing stimuli
%% Todo - change this!
%------------------------
stim_t = cell(numTrials,1);
for i_trial = startTrial:numTrials
      
    Screen('FillPoly',wPtr,black,triggerSize);  %for photodiode recording trigger
    Screen('Textsize', wPtr, 36); 
    Screen('DrawText', wPtr, '+', center(1), center(2), [0, 0, 0, 255]);
    Screen('Flip',wPtr,[],0);
    
    tempRand = randi(50)*.01;
    WaitSecs(iSi+tempRand);
    
   
    
    % Set Frequency of Left and Right for Trial
    %------------------------------------------
    f_left = freqs(1); %Frequency in Hz to approximate
    f_left_actual = (1/refresh)/round((1/refresh)/(f_left));  %closest refresh cycle approximation
    frameSwitchLeft =  ceil(round((1/(refresh*f_left_actual)))/2);% Frames on which the image switches i.e. 60 = the image switches on every 60th frame
    left_isi = 1/f_left_actual;

    f_right = freqs(2); %Frequency in Hz to approximate
    f_right_actual = (1/refresh)/round((1/refresh)/(f_right));  %closest refresh cycle approximation
    frameSwitchRight =  ceil(round((1/(refresh*f_right_actual)))/2);% Frames on which the image switches i.e. 60 = the image switches on every 60th frame
    right_isi = 1/f_right_actual;    
    
    % Determine fixation and more instructional text
    %------------------------------------------
    fixation = [center(1)-floor((imageSizeW)+100 ) center(2)-floor((imageSizeH/2))] ; %Remember to set original fixations (lines 23 & 24)
    fixation2 = [center(1)+ 100 center(2)-floor((imageSizeH/2))] ; %Remember to set original fixations (lines 23 & 24)   
    trial_start = GetSecs(); 
        
    % Start trial timing and cue the photodiode
    %------------------------------------------
    stim_t{i_trial} = zeros(2,2000);
    i_stimLeft = 0;
    i_stimRight = 0;     
    time_flagLeft = 0; 
    time_flagRight = 0; 
  
    % Display the Stimuli 
    %---------------- --------------------------    
    spriteRect = [0 0 imageSizeW imageSizeH];
    loc1 = [fixation(1) fixation(2) (fixation(1)+imageSizeW) (fixation(2)+imageSizeH) ];
    loc2 = [fixation2(1) fixation2(2) (fixation2(1)+imageSizeW) (fixation2(2)+imageSizeH) ];
    Screen('FillRect', wPtr ,black,loc1);
    Screen('FillRect', wPtr ,black,loc2);
    Screen('Textsize', wPtr, 36); 
    Screen('Flip', wPtr,[],    0); %flip it to the screen, leave the screen as it was besides 
    
    % Parameters to keep which image 
    %------------------------------------------
    tickerLeft = 0;
    tickerRight = 0;
    iter = 1; % "iter" just counts the frame number  

    while (GetSecs() - trial_start < trial_length)
        % Decide which image to show on the left hand side
        if ~mod(iter / frameSwitchLeft, 1) == 1
            if tickerLeft == 0
                tickerLeft = 1;
            elseif tickerLeft == 1
                tickerLeft = 0;
                time_flagLeft = 1;
            end
        end
        if ~mod(iter / frameSwitchRight, 1) == 1
            if tickerRight == 0
                tickerRight = 1;
            elseif tickerRight == 1
                tickerRight = 0;
                time_flagRight = 1;
            end
        end
       % Decide which of our two textures to show
        if tickerLeft == 0
            Screen('FillRect', wPtr, white, loc1);
            Screen('FillPoly',wPtr,white,triggerSize);  %for photodiode recording trigger
        elseif tickerLeft == 1
            Screen('FillRect', wPtr, black, loc1);
            Screen('FillPoly',wPtr,white,triggerSize);  %for photodiode recording trigger
        end

        if tickerRight == 0
            Screen('FillRect', wPtr, black, loc2);
            Screen('FillPoly',wPtr,white,triggerSize);  %for photodiode recording trigger
        elseif tickerRight == 1
            Screen('FillRect', wPtr, white, loc2);
            Screen('FillPoly',wPtr,white,triggerSize);  %for photodiode recording trigger
        end
        [VBL screen_on] = Screen('Flip', wPtr);  % Lazy flip to the screen
        iter = iter + 1 ;     % Increment the frame counter
        if time_flagLeft == 1   %to test the timing 
            i_stimLeft = i_stimLeft + 1;
            stim_t{i_trial}(1,i_stimLeft) = screen_on;
            time_flagLeft = 0;
        end                   
        if time_flagRight == 1   %to test the timing
            i_stimRight = i_stimRight + 1;
            stim_t{i_trial}(1,i_stimRight) = screen_on;
            time_flagRight = 0; 
        end
    end   
end

%% CLOSE EVERYTHING
KbWait; 
pause(1)
sca
