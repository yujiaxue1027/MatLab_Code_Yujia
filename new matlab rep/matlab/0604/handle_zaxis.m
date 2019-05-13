% clear; close all; clc;
global handle_zaxis_UPDN; % make h a global variable so it can be used outside the main
          % function. Useful when you do event handling and sequential           move
%% Create Matlab Figure Container
fpos    = get(0,'DefaultFigurePosition'); % figure default position
fpos(3) = 650; % figure window size;Width
fpos(4) = 450; % Height
 
f = figure('Position', fpos,...
           'Menu','None',...
           'Name','APT GUI');
%% Create ActiveX Controller
handle_zaxis_UPDN = actxcontrol('MGMOTOR.MGMotorCtrl.1',[20 20 600 400 ], f);
 
%% Initialize
% Start Control
handle_zaxis_UPDN.StartCtrl;
 
% Set the Serial Number
SN = 27252217; % put in the serial number of the hardware
set(handle_zaxis_UPDN,'HWSerialNum', SN);
 
% Indentify the device
handle_zaxis_UPDN.Identify;
 
pause(5); % waiting for the GUI to load up;
%% Controlling the Hardware
%h.MoveHome(0,0); % Home the stage. First 0 is the channel ID (channel 1)
                 % second 0 is to move immediately
%% Event Handling
handle_zaxis_UPDN.registerevent({'MoveComplete' 'MoveCompleteHandler'});
 
%% Sending Moving Commands
timeout = 10; % timeout for waiting the move to be completed
%h.MoveJog(0,1); % Jog
 
t1 = clock; % current time
while(etime(clock,t1)<timeout) 
% wait while the motor is active; timeout to avoid dead loop
    s = handle_zaxis_UPDN.GetStatusBits_Bits(0);
    if (IsMoving(s) == 0)
      pause(2); % pause 2 seconds;
      handle_zaxis_UPDN.MoveHome(0,0);
      disp('Home Started!');
      break;
    end
end

