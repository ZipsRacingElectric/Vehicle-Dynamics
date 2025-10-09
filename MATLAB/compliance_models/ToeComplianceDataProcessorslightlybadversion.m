

% This Toe Compliance calculator uses deflection and moment measures to
% calculate the toe compliance in degrees/Nm

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

%% 1. Pull in data from excel
% Here I use the "readtable" function. This pulls the data in as a table
% and automatically detecs headers making it easy to call values or
% specific index values

Excel = readtable("C:\Users\ATuck\Documents\GitHub\Vehicle-Dynamics\MATLAB\compliance_models\Compliance Testing ZR25.xlsx")

%% 2. Calculate moments
% To calculate our applied moments we need to know two variables. our lever
% arm length (L) and our applied force (F). 

% For this test our lever arm
% length never changed and the values is defined in the excel table. We can
% start by assinging this value a variable name.

dLeverArm = 1 % Lever arm length in meters.

% This is when we need to start being aware of what units we are working
% with. We know our final output should be in degrees/Nm. Because of this
% our distance values should all be in meters, our angle values should
% be in degrees, and our force values in newtons.

% Next we need to define our forces. These are stored in a column titled Kg
% We can load the four differnet load cases into a matrix under one
% variable name.

Forces = Excel.kg(1:4)

% This is an alternative way to call the same values by indexing their
% cells. Forces = Excel{1:4,1} 1:4 means rows 1 through 4 the , 1 pulls
% those rows from column 1

% Now we are ready to find our moments. Lets use a function for this. 

ToeMoments = GetMoment(Forces,dLeverArm)

%% moment function
% This function needs fed forces (FKg) and lever arm length (LAD)
% We start by writing what the output will be, in this case "moments".
% Next we define the name of our function "GetMoment" and put the
% arguements it will need in  parenthasis (FKg, LAD)

function moments = GetMoment(FKg, LAD) % get moments back in Nm

% We then write the driving equation using the variables defined in the title. 
    moments = FKg*LAD*9.81 % {kg*m*m/s^2) = Nm

end

%% 3. Find displacements of dial indicators
% to find the diplacements we need to calulate the difference between the
% movement of the two indicators for each of the four applied moment
% values. We do this for both sides of the car. This leaves us with
% 8 displacements to find from 16 numbers. Lets use another function to speed it up and
% then package all of the values in 2 easy to work with matrices, one for
% each wheel.

% First, assign variable names to the values of each dial indicator. The
% naming convention will be DI1 (displacement indicator one) etc...

DI1 = Excel{1:4,2}
DI2 = Excel{1:4,3}
DI3 = Excel{1:4,4}
DI4 = Excel{1:4,5}

% now lets build our two matrices both 

% In the excel file these values were recorded in hundredths of an inch. We
% need to convert this to meters within our displacements function
% lets call the function and build our two, [4x1] matrices

DR = GetDisplacements(DI1, DI2) % displacement at corner
DL = GetDisplacements(DI3, DI4) % displacement at adjacent corner

% Using a function here has saved us 16 calculations

%% diaplacement function
function distance = GetDisplacements(D1in, D2in)

distance = D1in-D2in
% convert to meters
distance = distance*.000254 % from hundredths to meters

end

%% 4. Find Toe Angles
% we almost have all the pieces we need to find our compliance. We have the
% values for the bottom of our triangles [DR, DL], the second side length of our
% triangle is defined in the excel, We have our corresponding forces [ToeMoments], 
% the last piece is finding the angle change using trig and our two triangle 
% sides. 

% Pull distance between the two indicators as measured in the setup this is
% the second tringle side.This was measured in inches so must also be
% converted to meters.

d = Excel{3,6} 

% this pulls a single cell which cant easily be divided with the 4x1
% displacement matrices. Lets fix this by making into a repeating 4x1
% matric so dimensions match
             
val = str2double(d{1});       % convert '15.5' → 15.5 numeric
DY = repmat(val, 4, 1);       % make 4x1 numeric column



% now we can solve for the angles at each force interval. Naming convention
% toe angle on right side = TARS
TARS = GetCompAng(DY, DR); % corner side theta calc
TALS = GetCompAng(DY, DL); % adjacent side theta calc

%% Compliance angle function
function theta = GetCompAng(Dispx,Dispy)

theta = rad2deg(atan(Dispy./Dispx));


end
%% Find compliance for each side (degrees per Nm) 
% we can divide the angle change by the corresponding force to find the
% angle change per applied Nm of force.

CRS = GetToeCompliance(TARS, ToeMoments); % compliance at test corner
CLS = GetToeCompliance(TALS, ToeMoments) ;% compliance at adjacent corner
%% compliance function

function TC = GetToeCompliance(ang,Moment) % outputs degrees per Kg

    TC = ang./Moment

end

%% 5. Plot the reults (should be linearish)

figure

subplot(2,1,1)
plot(ToeMoments, TARS)
title("Toe Compliance RS")
xlabel('Nm')
ylabel("Degrees")

subplot(2,1,2)
plot(ToeMoments, TALS)
title("Toe Compliance LS")
xlabel('Nm')
ylabel("Degrees")
 









