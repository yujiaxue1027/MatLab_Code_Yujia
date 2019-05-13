%% example code to control thorlabs motorized translation stage

%% initialize handles for horizontal stage and vertical stage
handle_horizontal;
handle_vertical;

%% move to a certain location
horizontal_pos = -4.0292;
vertical_pos = -3.9771;

%% handle_horizontal_LPRN and handle_vertical_UPDN are 
%% global variables defined in handle_horizontal.m and handle_vertical.m
handle_horizontal_LPRN.SetAbsMovePos(0,horizontal_pos);% set absolute horizontal location
handle_horizontal_LPRN.MoveAbsolute(0,1==0);% move to the set location

handle_vertical_UPDN.SetAbsMovePos(0,vertical_pos);
handle_vertical_UPDN.MoveAbsolute(0,1==0);