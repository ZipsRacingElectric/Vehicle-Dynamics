
% This Toe Compliance calculator uses deflection and moment measurements to
% calculate the toe compliance in degrees/Nm
% Written by Abigail Tucker - 10/10/25


% Process:
% Pull in Excel Data -> 
% Moment calculations -> 
% displacement calculations ->
% Find toe angle -> 
% Find compliance -> 
% plot results🐸

%% This starts with the front left [FL] wheel and a moment causing 
% positive toe angle

clear, clc
close all

%% 1. Pull in data from excel

filePath = "C:\Users\ATuck\OneDrive - The University of Akron\Zips Racing FSAE - ZR26\Vehicle Dynamics\200 Controls\Kinematics+Compliance Testing\Toe Compliance\Toe Compliance ZR25 (example).xlsx";

Excel = readtable(filePath);
FrontPositiveToe = readtable(filePath, "Range", "FrontPositive");
FrontNegativeToe = readtable(filePath, "Range", "FrontNegative");
RearPositiveToe = readtable(filePath, "Range", "RearPositive");
RearNegativeToe = readtable(filePath, "Range", "RearNegative");

%% 2. Calculate moments

ArmLength = Excel.LengthOfMomentArm_m_All(1);
LeverArm =  ones(4,1).*ArmLength; % Lever arm length in meters.

% Define Forces
ForcesFP = FrontPositiveToe.AppliedForces_Kgf_FP;
ForcesFN = FrontNegativeToe.AppliedForces_Kgf_FN;
ForcesRP = RearPositiveToe.AppliedForces_Kgf_RP;
ForcesRN = RearNegativeToe.AppliedForces_Kgf_RN;

% Now we are ready to find our moments. Lets use a function for this. 

ToeMomentsFP = GetToeMoment(ForcesFP,LeverArm);
ToeMomentsFN = GetToeMoment(ForcesFN,LeverArm);
ToeMomentsRP = GetToeMoment(ForcesRP,LeverArm);
ToeMomentsRN = GetToeMoment(ForcesRN,LeverArm);


%% 3. Find displacements of dial indicators

% Front Positive
Dial1FP = FrontPositiveToe.DialIndicator1__001In_FP;
Dial2FP = FrontPositiveToe.DialIndicator2__001In_FP;
Dial3FP = FrontPositiveToe.DialIndicator3__001In_FP;
Dial4FP = FrontPositiveToe.DialIndicator4__001In_FP;

% Front Negative
Dial1FN = FrontNegativeToe.DialIndicator1__001In_FN;
Dial2FN = FrontNegativeToe.DialIndicator2__001In_FN;
Dial3FN = FrontNegativeToe.DialIndicator3__001In_FN;
Dial4FN = FrontNegativeToe.DialIndicator4__001In_FN;

% Rear Positive
Dial1RP = RearPositiveToe.DialIndicator1__001In_RP;
Dial2RP = RearPositiveToe.DialIndicator2__001In_RP;

% Rear Negative
Dial1RN = RearNegativeToe.DialIndicator1__001In_RN;
Dial2RN = RearNegativeToe.DialIndicator2__001In_RN;


% In the excel file these values were recorded in thousandths of an inch. We
% need to convert this to meters within our displacements function
% lets call the function and build our two, [4x1] matrices

DisplacementFRP = GetToeDisplacements(Dial1FP, Dial2FP) ;% displacement at corner
DisplacementFLP = GetToeDisplacements(Dial3FP, Dial4FP) ;% displacement at adjacent corner

DisplacementFRN = GetToeDisplacements(Dial1FN, Dial2FN) ;% displacement at corner
DisplacementFLN = GetToeDisplacements(Dial3FN, Dial4FN) ;% displacement at adjacent corner

DisplacementRP = GetToeDisplacements(Dial1RP, Dial2RP) ;% rear positive

DisplacementRN = GetToeDisplacements(Dial1RN, Dial2RN) ;% rear negative




%% 4. Find Toe Angles
HoriontalDistance = Excel.HorizontalDistanceBetweenDialIndicators_m_All(1);
ToeBase =  ones(4,1).*HoriontalDistance;

% now we can solve for the angles at each force interval. Naming convention
% Toe angle on front right side, positive = TAFRP

% Front Positive
TAFRP = GetToeComplianceAng(DisplacementFRP, ToeBase); % corner side theta calc
TAFLP = GetToeComplianceAng(DisplacementFLP, ToeBase); % adjacent side theta calc

% Front Negative
TAFRN = GetToeComplianceAng(DisplacementFRN, ToeBase);
TAFLN = GetToeComplianceAng(DisplacementFLN, ToeBase);

% Rear Positive
TARP = GetToeComplianceAng(DisplacementRP, ToeBase);

% Rear Negative
TARN = GetToeComplianceAng(DisplacementRN, ToeBase);


%% Find compliance for each side (degrees per Nm) 

% Front Positive
CFRP = GetToeCompliance(TAFRP, ToeMomentsFP); % compliance at test corner
CFLP = GetToeCompliance(TAFLP, ToeMomentsFP); % compliance at adjacent corner

% Front Negative
CFRN = GetToeCompliance(TAFRN, ToeMomentsFN);
CFLN = GetToeCompliance(TAFLN, ToeMomentsFN);

% Rear Postive
CRP = GetToeCompliance(TARP, ToeMomentsRP);

% Rear Negative
CRN = GetToeCompliance(TARN, ToeMomentsRN);


%% 5. Plot the reults (should be linearish)

figure(1)

% Define colors and markers for each curve
colors = {'b','r','g','m','c','k'};
markers = {'*','o','s','^','d','x'};

% Front Right Positive
subplot(2,3,1)
plot(ToeMomentsFP, TAFRP, 'Color', colors{1}, 'Marker', markers{1}, 'LineWidth',1.5)
title("Front Right Positive")
xlabel('Nm'); ylabel("Degrees"); grid on

% Front Left Positive
subplot(2,3,4)
plot(ToeMomentsFP, TAFLP, 'Color', colors{2}, 'Marker', markers{2}, 'LineWidth',1.5)
title("Front Left Positive")
xlabel('Nm'); ylabel("Degrees"); grid on

% Front Right Negative
subplot(2,3,2)
plot(ToeMomentsFN, TAFRN, 'Color', colors{3}, 'Marker', markers{3}, 'LineWidth',1.5)
title("Front Right Negative")
xlabel('Nm'); ylabel("Degrees"); grid on

% Front Left Negative
subplot(2,3,5)
plot(ToeMomentsFN, TAFLN, 'Color', colors{4}, 'Marker', markers{4}, 'LineWidth',1.5)
title("Front Left Negative")
xlabel('Nm'); ylabel("Degrees"); grid on

% Rear Positive
subplot(2,3,3)
plot(ToeMomentsRP, TARP, 'Color', colors{5}, 'Marker', markers{5}, 'LineWidth',1.5)
title("Rear Positive")
xlabel('Nm'); ylabel("Degrees"); grid on

% Rear Negative
subplot(2,3,6)
plot(ToeMomentsRN, TARN, 'Color', colors{6}, 'Marker', markers{6}, 'LineWidth',1.5)
title("Rear Negative")
xlabel('Nm'); ylabel("Degrees"); grid on

% Add a single title above all subplots
sgtitle('Toe Compliance by Applied Moment (abs)')

figure(2)
hold on

% Colors and markers
colors = {'b','r','g','m','c','k'};
markers = {'*','o','s','^','d','x'};

% Plot each curve
plot(ToeMomentsFP, TAFRP, 'Color', colors{1}, 'Marker', markers{1}, 'LineWidth',2)
plot(ToeMomentsFP, TAFLP, 'Color', colors{2}, 'Marker', markers{2}, 'LineWidth',2)
plot(ToeMomentsFN, TAFRN, 'Color', colors{3}, 'Marker', markers{3}, 'LineWidth',2)
plot(ToeMomentsFN, TAFLN, 'Color', colors{4}, 'Marker', markers{4}, 'LineWidth',2)
plot(ToeMomentsRP, TARP, 'Color', colors{5}, 'Marker', markers{5}, 'LineWidth',2)
plot(ToeMomentsRN, TARN, 'Color', colors{6}, 'Marker', markers{6}, 'LineWidth',2)

% Labels, title, and grid
xlabel('Nm')
ylabel('Toe Angle [deg]')
title('Toe compliance Comparison All Corners')
grid on

% Legend
legend({'Front Right Positive', 'Front Left Positive', ...
        'Front Right Negative', 'Front Left Negative', ...
        'Rear Positive', 'Rear Negative'}, 'Location','best')

hold off



%% displacement function Positive
function distance = GetToeDisplacements(D1in, D2in)

    distance = abs(D1in)+abs(D2in);
    % convert to meters
    distance = distance*.0000254; % from thousandths to meters

end

%% compliance function

function C = GetToeCompliance(ang,Moment) % outputs degrees per Kg

    C = ang./Moment;

end
%% Compliance angle function
function theta = GetToeComplianceAng(Dispx,Dispy)

    theta = rad2deg(atan(Dispx./Dispy)); 
    % here we use ./ to divide each element of the matrix by the scalar. 
end

%% moment function
% This function needs fed mass (Kg) and lever arm length (LAD)
% We start by writing what the output will be, in this case "moments".
% Next we define the name of our function "GetMoment" and put the
% arguements it will need in  parenthesis (FKg, LAD)

function moments = GetToeMoment(Kg, LAD) % get moments back in Nm

% We then write the driving equation using the variables defined in the title. 
    moments = Kg.*LAD*9.81; % {kgf*m*m/s^2) = Nm
end
% This function needs fed mass (Kg) and lever arm length (LAD)
% We start by writing what the output will be, in this case "moments".
% Next we define the name of our function "GetMoment" and put the
% arguements it will need in  parenthesis (FKg, LAD)

function moments = GetMoment(Kg, LAD) % get moments back in Nm

% We then write the driving equation using the variables defined in the title. 
    moments = Kg*LAD*9.81; % {kg*m*m/s^2) = Nm
end