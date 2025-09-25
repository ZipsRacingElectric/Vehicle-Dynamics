% calculates lateral load distribution stuff based on current parameters
% LLD is chosen to give peak lateral grip at a slightly negative stability
% index - 
% reference:


clear

%% load in parameters file
%filepath = 'C:\Users\benmo\OneDrive - The University of Akron\Documents - Zips Racing FSAE\ZR25\Vehicle Dynamics\System\Analysis\'
filepath = '\vehicle_data\'
filename = 'zr25_data.xlsx'

fileJawn = strcat(filepath,filename);

% Car parameters file, pull all data from here, don't input manually!
% this helps keep parameters consistent between iterations, update
% everything in the parameters excel file. Save after editing!
zr = vehicle(fileJawn);

g = zr.g;

%% prototype tire model, very basic for now
% Basic Tire Parameters, from chassissim guy
    a_x = [12.476 82.761 -1592.7 4985.8];

    fTire = tire(2.2,0.000443656,a_x);
    rTire = tire(2.2,0.000443656,a_x);



%% Values I wanna play with
mass = zr.mass_total;
cgHeight = zr.cg_height;
frontBias = zr.front_mass_distribution;
ClA = zr.cla_at_skidpad;
aeroBias = zr.cop_at_skidpad;
rwd_pct = 1;
rSA = 5;
targetLatG = 1.5;
frontSpringRate = zr.spring_rate_front;
reaspring_rate_rear = zr.spring_rate_rear;

roll_center_front = zr.roll_center_front;
roll_center_rear = zr.roll_center_rear;

frontBarRate = zr.bar_spring_rate_front;
reabar_spring_rate_rear = zr.bar_spring_rate_rear;



% Control Run- Current car parameters

    % roll axis to CG height at CG
    rcm = zr.roll_center_front + (1-zr.front_mass_distribution)*(zr.roll_center_rear-zr.roll_center_front);
    hsm = zr.cg_height-rcm;
    
    % Lateral Load Transfer for current trial - "Magic Number"
    prm = zr.track_width_front*zr.total_roll_rate_front/(zr.track_width_front*zr.total_roll_rate_front+zr.track_width_rear*zr.total_roll_rate_rear)
    prr = (zr.front_mass_distribution*zr.roll_center_front+prm*hsm)/zr.cg_height
    
    
    % Roll Distribution
    rollDistr = (zr.front_mass_distribution * zr.roll_center_front*zr.average_track_width/zr.track_width_front + prm*hsm)/zr.cg_height



%% Trial setup
% velocity reference value- keep this the same each trial cause I don't 
% want to deal with this problem as a diff eq and have to converge, etc.
velocity = 11;
zr.velocity_skidpad = velocity;

lltArray = linspace(.1,1,10);

%Trial Value - Start with target value, and run trials from -20% to +20%
% of target value. 
nTrials = 21;
trialPercent = 1-linspace(.2,-.2,nTrials);
trialRsr = zeros(2,nTrials);    % spring rates
trialFsr = zeros(2,nTrials);
trialRrc = zeros(2,nTrials);    % roll centers
rrialFrc = zeros(2,nTrials);
trialRbr = zeros(2,nTrials);    % bar rates
trialFbr = zeros(2,nTrials);

%% Rear Sping Rate
for i = 1:1:nTrials
    % Set trial Value and update
    zr.spring_rate_rear = trialPercent(i)*reaspring_rate_rear;
    zr.update;

    % roll axis to CG height at CG
    rcm = zr.roll_center_front + (1-zr.front_mass_distribution)*(zr.roll_center_rear-zr.roll_center_front);
    hsm = zr.cg_height-rcm;
    
    % Lateral Load Transfer for current trial - "Magic Number"
    prm = zr.track_width_front*zr.total_roll_rate_front/(zr.track_width_front*zr.total_roll_rate_front+zr.track_width_rear*zr.total_roll_rate_rear);
    
    % Roll Distribution
    rollDistr = (zr.front_mass_distribution * zr.roll_center_front*zr.average_track_width/zr.track_width_front + prm*hsm)/zr.cg_height;

    trialRsr(1,i) = prm;
    trialRsr(2,i) = rollDistr;
end
zr.reset;

%% Front Spring Rate
for i = 1:1:nTrials
    % Set trial Value and update
    zr.spring_rate_front = trialPercent(i)*frontSpringRate;
    zr.update;

    % roll axis to CG height at CG
    rcm = zr.roll_center_front + (1-zr.front_mass_distribution)*(zr.roll_center_rear-zr.roll_center_front);
    hsm = zr.cg_height-rcm;
    
    % Lateral Load Transfer for current trial - "Magic Number"
    prm = zr.track_width_front*zr.total_roll_rate_front/(zr.track_width_front*zr.total_roll_rate_front+zr.track_width_rear*zr.total_roll_rate_rear);
    
    % Roll Distribution
    rollDistr = (zr.front_mass_distribution * zr.roll_center_front*zr.average_track_width/zr.track_width_front + prm*hsm)/zr.cg_height;
 
    trialFsr(1,i) = prm;
    trialFsr(2,i) = rollDistr;
end
zr.reset;

%% Rear Roll Center
for i = 1:1:nTrials
    % Set trial Value and update
    zr.roll_center_rear = trialPercent(i)*roll_center_rear;
    zr.update;

    % roll axis to CG height at CG
    rcm = zr.roll_center_front + (1-zr.front_mass_distribution)*(zr.roll_center_rear-zr.roll_center_front);
    hsm = zr.cg_height-rcm;
    
    % Lateral Load Transfer for current trial - "Magic Number"
    prm = zr.track_width_front*zr.total_roll_rate_front/(zr.track_width_front*zr.total_roll_rate_front+zr.track_width_rear*zr.total_roll_rate_rear);
    
    % Roll Distribution

    rollDistr = (zr.front_mass_distribution * zr.roll_center_front*zr.average_track_width/zr.track_width_front + prm*hsm)/zr.cg_height;
 
    trialRrc(1,i) = prm;
    trialRrc(2,i) = rollDistr;
end
zr.reset;

%% Front Roll Center
for i = 1:1:nTrials
    % Set trial Value and update
    zr.roll_center_front = trialPercent(i)*roll_center_front;
    zr.update;

    % roll axis to CG height at CG
    rcm = zr.roll_center_front + (1-zr.front_mass_distribution)*(zr.roll_center_rear-zr.roll_center_front);
    hsm = zr.cg_height-rcm;
    
    % Lateral Load Transfer for current trial - "Magic Number"
    prm = zr.track_width_front*zr.total_roll_rate_front/(zr.track_width_front*zr.total_roll_rate_front+zr.track_width_rear*zr.total_roll_rate_rear);
    
    % Roll Distribution
    rollDistr = (zr.front_mass_distribution * zr.roll_center_front*zr.average_track_width/zr.track_width_front + prm*hsm)/zr.cg_height;
    
    trialFrc(1,i) = prm;
    trialFrc(2,i) = rollDistr;
end
zr.reset;

%% Rear Bar Rate
for i = 1:1:nTrials
    % Set trial Value and update
    zr.bar_spring_rate_rear = trialPercent(i)*reabar_spring_rate_rear;
    zr.update;

    % roll axis to CG height at CG
    rcm = zr.roll_center_front + (1-zr.front_mass_distribution)*(zr.roll_center_rear-zr.roll_center_front);
    hsm = zr.cg_height-rcm;
    
    % Lateral Load Transfer for current trial - "Magic Number"
    prm = zr.track_width_front*zr.total_roll_rate_front/(zr.track_width_front*zr.total_roll_rate_front+zr.track_width_rear*zr.total_roll_rate_rear);
    
    % Roll Distribution
    rollDistr = (zr.front_mass_distribution * zr.roll_center_front*zr.average_track_width/zr.track_width_front + prm*hsm)/zr.cg_height;

    trialRbr(1,i) = prm;
    trialRbr(2,i) = rollDistr;
end
zr.reset;

%% Front Bar Rate
for i = 1:1:nTrials
    % Set trial Value and update
    zr.bar_spring_rate_front = trialPercent(i)*frontBarRate;
    zr.update;

    % roll axis to CG height at CG
    rcm = zr.roll_center_front + (1-zr.front_mass_distribution)*(zr.roll_center_rear-zr.roll_center_front);
    hsm = zr.cg_height-rcm;
    
    % Lateral Load Transfer for current trial - "Magic Number"
    prm = zr.track_width_front*zr.total_roll_rate_front/(zr.track_width_front*zr.total_roll_rate_front+zr.track_width_rear*zr.total_roll_rate_rear);
    
    % Roll Distribution
    rollDistr = (zr.front_mass_distribution * zr.roll_center_front*zr.average_track_width/zr.track_width_front + prm*hsm)/zr.cg_height;
 
    trialFbr(1,i) = prm;
    trialFbr(2,i) = rollDistr;
end
zr.reset;



%% Plots
% Stable Lat G
figure(1)
subplot(2,1,1)
hold on
plot(trialPercent,trialRsr(1,:))
plot(trialPercent,trialFsr(1,:))
plot(trialPercent,trialRrc(1,:))
plot(trialPercent,trialFrc(1,:))
plot(trialPercent,trialRbr(1,:))
plot(trialPercent,trialFbr(1,:))
grid on

legend('Rear Spring Rate','Front Spring Rate','Rear Roll Center','FrontRollCenter','Rear Bar Rate','Front Bar Rate')
xlabel('% change from starting point')
ylabel('PRM')
title('Which parameters are a bigger wrench to affect car balance')
hold off


subplot(2,1,2)
hold on
plot(trialPercent,trialRsr(2,:))
plot(trialPercent,trialFsr(2,:))
plot(trialPercent,trialRrc(2,:))
plot(trialPercent,trialFrc(2,:))
plot(trialPercent,trialRbr(2,:))
plot(trialPercent,trialFbr(2,:))

legend('Rear Spring Rate','Front Spring Rate','Rear Roll Center','FrontRollCenter','Rear Bar Rate','Front Bar Rate')
xlabel('% change from starting point')
ylabel('Roll Distribution')
grid on
hold off


%% Functions to hide mess

function [frontSlipSlope,rearSlipSlope] = calcSlipSlopes(zr,fTire,rTire,lltArray,rSA,rwd_pct,totalLatForce)
    % Calculate front slip angle
    fSA = rSA*(zr.b*(rwd_pct*rTire.staticForce  - 2*((1-lltArray).^2).*rwd_pct*rTire.ka*rTire.kb .* (totalLatForce*9.81*zr.cgTrackRatio).^2));
    fSA = fSA./ (zr.a*(fTire.staticForce - 2*(lltArray.^2).*fTire.ka*fTire.kb.*(totalLatForce*9.81*zr.cgTrackRatio).^2));
   
    % Slope of slip angles
    frontSlipSlope = fTire.a_1 + 2*fTire.a_2.*(fSA./57.3) + 3*fTire.a_3.*(fSA./57.3).^2 + 4*fTire.a_4.*(fSA./57.3).^3;
    % not sure what these constants do
    rearSlipSlope = 14.323 -1306.2.*(rSA./57.3).^2 * ones(1,length(lltArray));  
end

