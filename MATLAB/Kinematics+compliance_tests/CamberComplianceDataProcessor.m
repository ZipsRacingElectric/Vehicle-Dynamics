% This Camber Compliance calculator uses deflection and moment measurements to
% calculate the camber compliance in degrees/Nm
% tutorial code Written by Abigail Tucker - 10/10/25 


% Process:
% Pull in Excel Data -> 
% Moment calculations -> 
% displacement calculations ->
% Find camber angle -> 
% Find compliance -> 
% plot results🐸
%% What you need to do if using this for compliance testing!!

% 1) Fill out camber compliance excel sheet with correct uints. dont change
%names or placement
% 2) Copy and paste your own filepath to that sheet below

%%
clear, clc, close all

%% 1. Pull in data from excel
% Here I use the "readtable" function. This pulls the data in as a table
% and automatically detects headers making it easy to call values or
% specific index values

filePath = "C:\Users\ATuck\OneDrive - The University of Akron\Zips Racing FSAE - ZR26\Vehicle Dynamics\200 Controls\Kinematics+Compliance Testing\Camber Compliance\Camber Compliance ZR25 (example).xlsx";

Excel = readtable(filePath);
FrontPositiveCamber = readtable(filePath, "Range", "FrontPositive");
FrontNegativeCamber = readtable(filePath, "Range", "FrontNegative");
RearPositiveCamber = readtable(filePath, "Range", "RearPositive");
RearNegativeCamber = readtable(filePath, "Range", "RearNegative");

%% 2. Calculate moments
% To calculate our applied moments we need to know two variables. our lever
% arm length (L) and our applied force (F). 

% For this test our lever arm
% length never changed and the values is defined in the excel table. We can
% start by calling this value from the taable

ArmLength = Excel.LengthOfMomentArm_m_All(1);
LeverArm =  ones(4,1).*ArmLength; % Lever arm length in meters.

% This is when we need to start being aware of what units we are working
% with. We know our final output should be in degrees/Nm. Because of this
% our distance values should all be in meters, our angle values should
% be in degrees, and our force values in newtons.

% Next we need to define our forces. Our force values are stored in a column titled Kg
% We can load the four differnet load cases into a matrix under one
% variable name. I used named data sets in excel and read in as a table to
% keep headers

ForcesFP = FrontPositiveCamber.AppliedForces_Kgf_FP;
ForcesFN = FrontNegativeCamber.AppliedForces_Kgf_FN;
ForcesRP = RearPositiveCamber.AppliedForces_Kgf_RP;
ForcesRN = RearNegativeCamber.AppliedForces_Kgf_RN;

% Now we are ready to find our moments. Lets use a function for this. 

CamberMomentsFP = GetMoment(ForcesFP,LeverArm);
CamberMomentsFN = GetMoment(ForcesFN,LeverArm);
CamberMomentsRP = GetMoment(ForcesRP,LeverArm);
CamberMomentsRN = GetMoment(ForcesRN,LeverArm);

%% 3. Find displacements of dial indicators
% to find the diplacements we need to calculate the difference between the
% movement of the two indicators for each of the four applied moment
% values. We do this for both sides of the car. This leaves us with
% 8 displacements to find from 16 numbers. Lets use another function to speed it up and
% then package all of the values in 2 easy to work with matrices, one for
% each wheel.

% First, assign variable names to the values of each dial indicator. 

% Front Positive
Dial1FP = FrontPositiveCamber.DialIndicator1__001In_FP;
Dial2FP = FrontPositiveCamber.DialIndicator2__001In_FP;
Dial3FP = FrontPositiveCamber.DialIndicator3__001In_FP;
Dial4FP = FrontPositiveCamber.DialIndicator4__001In_FP;

% Front Negative
Dial1FN = FrontNegativeCamber.DialIndicator1__001In_FN;
Dial2FN = FrontNegativeCamber.DialIndicator2__001In_FN;
Dial3FN = FrontNegativeCamber.DialIndicator3__001In_FN;
Dial4FN = FrontNegativeCamber.DialIndicator4__001In_FN;

% Rear Positive
Dial1RP = RearPositiveCamber.DialIndicator1__001In_RP;
Dial2RP = RearPositiveCamber.DialIndicator2__001In_RP;

% Rear Negative
Dial1RN = RearNegativeCamber.DialIndicator1__001In_RN;
Dial2RN = RearNegativeCamber.DialIndicator2__001In_RN;


% In the excel file these values were recorded in thousandths of an inch. We
% need to convert this to meters within our displacements function
% lets call the function and build our two, [4x1] matrices

DisplacementFRP = GetDisplacements(Dial1FP, Dial2FP) ;% displacement at corner
DisplacementFLP = GetDisplacements(Dial3FP, Dial4FP) ;% displacement at adjacent corner

DisplacementFRN = GetDisplacements(Dial1FN, Dial2FN) ;% displacement at corner
DisplacementFLN = GetDisplacements(Dial3FN, Dial4FN) ;% displacement at adjacent corner

DisplacementRP = GetDisplacements(Dial1RP, Dial2RP) ;% rear positive

DisplacementRN = GetDisplacements(Dial1RN, Dial2RN) ;% rear negative

% Using a function here has saved us a lot lines of calculations

%% 4. Find Camber Angles
% we almost have all the pieces we need to find our compliance. We have the
% values for the bottom of our triangles (Displacements), the second side length of our
% triangle is defined in the excel, We have our corresponding forces [CamberMoments], 
% the last piece is finding the angle change using trig and our two triangle 
% sides. 

% Pull distance between the two indicators as measured in the setup this is
% the second tringle side.
             

VerticalDistance = Excel.VerticalDistanceBetweenDialIndicators_m_All(1);
CamberBase =  ones(4,1).*VerticalDistance;

% now we can solve for the angles at each force interval. Naming convention
% Camber angle on front right side, positive = CAFRP

% Front Positive
CAFRP = GetComplianceAng(CamberBase, DisplacementFRP); % corner side theta calc
CAFLP = GetComplianceAng(CamberBase, DisplacementFLP); % adjacent side theta calc

% Front Negative
CAFRN = GetComplianceAng(CamberBase, DisplacementFRN);
CAFLN = GetComplianceAng(CamberBase, DisplacementFLN);

% Rear Positive
CARP = GetComplianceAng(CamberBase, DisplacementRP);

% Rear Negative
CARN = GetComplianceAng(CamberBase, DisplacementRN);


%% Find compliance for each side (degrees per Nm) 
% we can divide the angle change by the corresponding force to find the
% angle change per applied Nm of force.

% Front Positive
CFRP = GetCompliance(CAFRP, CamberMomentsFP); % compliance at test corner
CFLP = GetCompliance(CAFLP, CamberMomentsFP); % compliance at adjacent corner

% Front Negative
CFRN = GetCompliance(CAFRN, CamberMomentsFN);
CFLN = GetCompliance(CAFLN, CamberMomentsFN);

% Rear Postive
CRP = GetCompliance(CARP, CamberMomentsRP);

% Rear Negative
CRN = GetCompliance(CARN, CamberMomentsRN);


%% 5. Plot the reults (should be linearish)

figure(1)

% Define colors and markers for each curve
colors = {'b','r','g','m','c','k'};
markers = {'*','o','s','^','d','x'};

% Front Right Positive
subplot(2,3,1)
plot(CamberMomentsFP, CAFRP, 'Color', colors{1}, 'Marker', markers{1}, 'LineWidth',1.5)
title("Front Right Positive")
xlabel('Nm'); ylabel("Degrees"); grid on

% Front Left Positive
subplot(2,3,4)
plot(CamberMomentsFP, CAFLP, 'Color', colors{2}, 'Marker', markers{2}, 'LineWidth',1.5)
title("Front Left Positive")
xlabel('Nm'); ylabel("Degrees"); grid on

% Front Right Negative
subplot(2,3,2)
plot(CamberMomentsFN, CAFRN, 'Color', colors{3}, 'Marker', markers{3}, 'LineWidth',1.5)
title("Front Right Negative")
xlabel('Nm'); ylabel("Degrees"); grid on

% Front Left Negative
subplot(2,3,5)
plot(CamberMomentsFN, CAFLN, 'Color', colors{4}, 'Marker', markers{4}, 'LineWidth',1.5)
title("Front Left Negative")
xlabel('Nm'); ylabel("Degrees"); grid on

% Rear Positive
subplot(2,3,3)
plot(CamberMomentsRP, CARP, 'Color', colors{5}, 'Marker', markers{5}, 'LineWidth',1.5)
title("Rear Positive")
xlabel('Nm'); ylabel("Degrees"); grid on

% Rear Negative
subplot(2,3,6)
plot(CamberMomentsRN, CARN, 'Color', colors{6}, 'Marker', markers{6}, 'LineWidth',1.5)
title("Rear Negative")
xlabel('Nm'); ylabel("Degrees"); grid on

% Add a single title above all subplots
sgtitle('Camber Compliance by Applied Moment (abs)')

figure(2)
hold on

% Colors and markers
colors = {'b','r','g','m','c','k'};
markers = {'*','o','s','^','d','x'};

% Plot each curve
plot(CamberMomentsFP, CAFRP, 'Color', colors{1}, 'Marker', markers{1}, 'LineWidth',2)
plot(CamberMomentsFP, CAFLP, 'Color', colors{2}, 'Marker', markers{2}, 'LineWidth',2)
plot(CamberMomentsFN, CAFRN, 'Color', colors{3}, 'Marker', markers{3}, 'LineWidth',2)
plot(CamberMomentsFN, CAFLN, 'Color', colors{4}, 'Marker', markers{4}, 'LineWidth',2)
plot(CamberMomentsRP, CARP, 'Color', colors{5}, 'Marker', markers{5}, 'LineWidth',2)
plot(CamberMomentsRN, CARN, 'Color', colors{6}, 'Marker', markers{6}, 'LineWidth',2)

% Labels, title, and grid
xlabel('Nm')
ylabel('Camber Angle [deg]')
title('Camber Compliance Comparison All Corners')
grid on

% Legend
legend({'Front Right Positive', 'Front Left Positive', ...
        'Front Right Negative', 'Front Left Negative', ...
        'Rear Positive', 'Rear Negative'}, 'Location','best')

hold off


%% displacement function
function distance = GetDisplacements(D1in, D2in)

    distance = abs(D1in)+abs(D2in);
    % convert to meters
    distance = distance*.0000254; % from thousandths to meters

end

%% compliance function

function C = GetCompliance(ang,Moment) % outputs degrees per Kg

    C = ang./Moment;

end

%% Compliance angle function
function theta = GetComplianceAng(Dispx,Dispy)

    theta = rad2deg(atan(Dispy./Dispx)); 
    % here we use ./ to divide each element of the matrix by the scalar. 
end

%% moment function
% This function needs fed mass (Kg) and lever arm length (LAD)
% We start by writing what the output will be, in this case "moments".
% Next we define the name of our function "GetMoment" and put the
% arguements it will need in  parenthesis (FKg, LAD)

function moments = GetMoment(Kg, LAD) % get moments back in Nm

% We then write the driving equation using the variables defined in the title. 
    moments = Kg.*LAD*9.81; % {kgf*m*m/s^2) = Nm
end