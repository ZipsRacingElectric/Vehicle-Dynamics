clear
%% Suspension force calc, even more slightly polished version
    % Written by Ben Model, 2025 VD Lead

    % Calculates the forces in each suspension member (CA's, push/tie rods)
    % for a number of defined load cases. These load cases represent
    % scenarios that we do not expect to frequently see on the car, but
    % ones for which that the suspension should be strong enough to take!

    % This matlab script takes values from a suspension design file, points
    % from a points file, and optionally writes the output forces to the
    % points file (zr25_SuspensionForces_heavy)

writeForcesToSpreadsheet = true;

addpath vehicle_data ;

folderWithDesignSpreadsheets = 'C:\Users\ATuck\OneDrive - The University of Akron\Zips Racing FSAE - ZR26\Vehicle Dynamics\100 Suspension\Suspension Points\';

githubFolder = '\vehicle_data\';
parameterSpreadsheet = strcat(githubFolder,'zr26_data.xlsx');
pointsSpreadsheet = strcat(folderWithDesignSpreadsheets,'ZR26_SuspensionForces.xlsx');



%% THERE IS A HUMAN IN THE CAR YOU ARE RESPONSIBLE FOR!!!!!!!!
% These can be tweaked once the team is confident that the current
% cases (from 2025, the first year this was used) provide sufficient
% margins of safety. Don't put a control arm through your friend's leg
% :)

%% Define Load Cases
    % Load cases, applied at contact patch:
        % 2G brake, FWD and RWD, aero
        % 2G accel, FWD and RWD. No Aero (can't do 2g accel at 90mph)
        % 1.9G Corner, aero
        % 3G Bump, aero
        % 2G Brake, backwards (with hub motors, this is the same as accel)
        % 1.75G corner, 3G bump, 2G brake
        % 3G Bump, 2G Brake backwards (50/50 bias)



LoadCases = [-1.6 0 1;        % 2g break, FWD
            1.6 0 1;          % 2g accel, FWD
            -1.6 0 1;         % 2g break, RWD
            1.6 0 1;          % 2g accel, RWD
            0 1.9 1;        % 1.9g Corner
            0 0 3;          % 3g bump
            -1.8 1.8 3;       % 1.8g brake, 1.8g corner, 3g bump, 50/50 bias
            2 1 3;          % backwards 2g brake, 1g corner, 3g bump, 50/50 bias
            -1.75 1.75 1;  % normal braking and turn, for steering forces
            0 1 1;          % 1g turn
            -1 0 1;         % 1g brake
            0 0 3];          % Pure 3g bump

bias = [1 1 0 0 .5 .5 .5 .5 .65 0 .7 0];        % Brake/torque bias: all forward (1), all rearwards (0)
direction = [1 1 1 1 1 1 1 -1 1 1 1 1];         % car traveling forwards (1) or backwards (-1)
aero = [1 0 1 0 1 1 1 0 1 1 1 1];               % Aero loads on (1), off (0), bias automatically calculated from CoP



%% Load Car Data
% Car parameters file, pull all data from here, don't input manually!
% this helps keep parameters consistent between iterations
% update everything in the parameters excel file. Save after editing!
% File also must be closed to run!!!!
zr = vehicle(parameterSpreadsheet);


%% Load Points Locations
FL = SuspensionPoints_better(pointsSpreadsheet,true);
FR = SuspensionPoints_better(pointsSpreadsheet,true);
FR.mirror;

RL = SuspensionPoints_better(pointsSpreadsheet,false);
RR = SuspensionPoints_better(pointsSpreadsheet,false);
RR.mirror;


%% Run Each Load Case
        % quick validation to yell at you if your load case vectors aren't the same
        % length
        if (length(LoadCases)-length(bias))~= 0 || (length(LoadCases)-length(direction))~= 0 || (length(LoadCases)-length(aero))~= 0
            error('Load Case Vector Lengths not consistent, dummy!!')
        end

for i=1:1:length(LoadCases)
    
    loadCase = LoadCases(i,:);

    % Loads at CP
    tireFz = getTireFz(zr,loadCase,aero(i));
    tireFy = getTireFy(zr,loadCase,tireFz); 
    tireFx = getTireFX(zr,loadCase,tireFz,bias(i),direction(i));
    % rewrite to array for each tire
    forcesFL = [tireFx.FL,tireFy.FL,tireFz.FL];
    forcesFR = [tireFx.FR,tireFy.FR,tireFz.FR];
    forcesRL = [tireFx.RL,tireFy.RL,tireFz.RL];
    forcesRR = [tireFx.RR,tireFy.RR,tireFz.RR];

        % Tuck into a struct for tracking
        forcesStruct.FL(i,:) = forcesFL;
        forcesStruct.FR(i,:) = forcesFR;
        forcesStruct.RL(i,:) = forcesRL;
        forcesStruct.RR(i,:) = forcesRR;

    % Solve for forces in suspension members
        % Left Front
        linkForcesFL(i) = solvySolve(forcesFL,FL);
    
        % Right Front
        linkForcesFR(i) = solvySolve(forcesFR,FR);
    
        % Left Rear
        linkForcesRL(i) = solvySolve(forcesRL,RL);
    
        % Right Rear
        linkForcesRR(i) = solvySolve(forcesRR,RR);
end


%% Paste results into excel spreadsheet
if writeForcesToSpreadsheet
    nLoadCases = length(LoadCases);
    startCell = 4;
     endCell = startCell + nLoadCases;

    cell =@(letter,number) strcat(letter,string(number));

    % also multiply by -1, for conventional tension/compression sign
    s = -1.*struct2table(linkForcesFL);
    writetable(s,pointsSpreadsheet,'Sheet','Suspension Forces' ,'Range',cell('J',startCell), 'WriteVariableNames', false)
    
    s = -1.*struct2table(linkForcesFR);
    writetable(s,pointsSpreadsheet,'Sheet','Suspension Forces','Range',cell('J',endCell),'WriteVariableNames',false)
    
    s = -1.*struct2table(linkForcesRL);
    writetable(s,pointsSpreadsheet,'Sheet','Suspension Forces','Range',cell('Q',startCell),'WriteVariableNames',false)
    
    s = -1.*struct2table(linkForcesRR);
    writetable(s,pointsSpreadsheet,'Sheet','Suspension Forces','Range',cell('Q',endCell),'WriteVariableNames',false)
    
    % wheel loads
    s = struct2table(forcesStruct)
    writetable(s,pointsSpreadsheet,'Sheet','Wheel Loads','Range','E2','WriteVariableNames',true)
    
    % Load Cases
    sp = char(ones([length(LoadCases),1]).*32);       % array of spaces
    c = char(ones([length(LoadCases),1]).*44);
    a = string(LoadCases(:,1)) +c+sp+ string(LoadCases(:,2)) +c+sp+ string(LoadCases(:,3));      % this is so dumb lmao
  
    s = array2table(a);
    writetable(s, pointsSpreadsheet, 'Sheet', 'Wheel Loads', 'Range', 'B3', 'WriteVariableNames', false);

    % --- Load Case Struct ---
    loadCaseStruct.Loads = a;
    loadCaseStruct.brakeBias = bias';
    loadCaseStruct.travelDirection = direction';
    loadCaseStruct.AeroApplied = aero';

   s = struct2table(loadCaseStruct);
    writetable(s, pointsSpreadsheet, 'Sheet', 'Suspension Forces', 'Range', 'C4', 'WriteVariableNames', false);

  
end



%% Functions
function susLinkForces = solvySolve(tireForces,corner)
    % Solve Linear System consisting of matrix and forces/moments on CP

    % the corner matrix is composed of the unit vectors of each suspension
    % member, and the cross product of their positions and unit vectors (this
    % comes straight from the static force/moment balance
    cornerMatrix = getCornerMatrix(corner);
    knowns = [tireForces,cross(corner.contactPatch',tireForces)]';

    solved = linsolve(cornerMatrix,knowns);

    susLinkForces.lcaFore = solved(1);
    susLinkForces.lcaAft = solved(2);
    susLinkForces.ucaFore = solved(3);
    susLinkForces.ucaAft = solved(4);
    susLinkForces.toe = solved(5);
    susLinkForces.push = solved(6);


end

function A = getCornerMatrix(corner)
    % takes a suspension corner object and returns matrix to calculate
    % forces for arbitrary applied load
    % transpose just for compactness in display
    
    % on chassis
    reactionPoints = [corner.lcaFore,corner.lcaAft,corner.ucaFore,corner.ucaAft,corner.toeInboard,corner.pushInboard];
    
    % make unit vectors for directions at reaction points
    lcaVectors = [corner.lcaFore,corner.lcaAft] - corner.lcaOutboard;
    ucaVectors = [corner.ucaFore,corner.ucaAft] - corner.ucaOutboard;
    toeVector = corner.toeInboard - corner.toeOutboard;
    pushVector = corner.pushInboard-corner.pushOutboard;
    
    linkUnitVectors = getUnitVectors([lcaVectors,ucaVectors,toeVector,pushVector]);
    
    momentVectors = cross(reactionPoints,linkUnitVectors);

    % build matrix
    A = [linkUnitVectors;momentVectors];

end


function unitVectors = getUnitVectors(arrayOfVectors)
% takes an 3 by n matrix, returns the unit vectors of each column
unitVectors = zeros(size(arrayOfVectors));
    for i =1:length(arrayOfVectors)
        unitVectors(:,i) = arrayOfVectors(:,i)./norm(arrayOfVectors(:,i));
    end
end


%% Vertical Load on Tire
function tireFz = getTireFz(car,loadCase,aeroBoolean)
    % this applies additional bump over 1g only to sprung mass!!
    % kind of still a rough simplification, tire loads are used for
    % calculating long. and lat. loads, but also vertical loads reacted by
    % control arms and push rods. This is a halfway point of using the
    % complete load at the contact patch, and recognizing that not all of
    % that load is reacted thru the CA's
    % more accuracy would require somewhat of a multi-body approach-
    % that sounds really hard.
    
    longG = loadCase(1);
    latG = loadCase(2);
    bump = loadCase(3);
    %% Static + Aero Loads + bump
    % 9.81*(mass over front corner * (addl. bump)*sprung mass) * aero load
    loadFL = car.g*(car.front_corner_mass + (bump-1)*car.sprung_mass_front) + aeroBoolean*car.downforce_at_max_velocity_front/2; 
    loadFR = car.g*(car.front_corner_mass + (bump-1)*car.sprung_mass_front) + aeroBoolean*car.downforce_at_max_velocity_front/2;
    loadRL = car.g*(car.rear_corner_mass +  (bump-1)*car.sprung_mass_rear)  + aeroBoolean*car.downforce_at_max_velocity_rear/2;
    loadRR = car.g*(car.rear_corner_mass +  (bump-1)*car.sprung_mass_rear)  + aeroBoolean*car.downforce_at_max_velocity_rear/2;

    %% Lateral Weight Transfer
    loadFL = loadFL + car.lltd_front *latG;
    loadFR = loadFR - car.lltd_front *latG;
    loadRL = loadRL + car.lltd_rear*latG;
    loadRR = loadRR - car.lltd_rear*latG;

    %% Longditudinal Weight Transfer
    % also tuck into a nice struct
    tireFz.FL = loadFL - car.long_load_transfer*longG/2;
    tireFz.FR = loadFR - car.long_load_transfer*longG/2;
    tireFz.RL = loadRL + car.long_load_transfer*longG/2;
    tireFz.RR = loadRR + car.long_load_transfer*longG/2;

    %% no fake downforce- make zero if negative
    if tireFz.FL <=0;tireFz.FL=0;end
    if tireFz.FR <=0;tireFz.FR=0;end
    if tireFz.RL <=0;tireFz.RL=0;end
    if tireFz.RR <=0;tireFz.RR=0;end


    
end
 
%% Lateral Load on Tire
function tireFy = getTireFy(car,loadCase,loads) 
    % calculates lateral forces generated by each tires by:
    % > find total lateral force, front and rear axle forces
    % Assume outside to inside Fy ratio is proportional to ratio of Fz
    % this should be a conservative assumption, as the derivative of the
    % max possible force decreases 

    totalFy = -loadCase(2)*car.mass_total*car.g;


    % solve simple force/moment balance for car at SS lat g
    eq = rref([car.a,-car.b,0;1,1,totalFy]);
    
    frontAxleFy = eq(1,3);
    rearAxleFy = eq(2,3);

    tireFy.FL = frontAxleFy * loads.FL / (loads.FL+loads.FR);
    tireFy.FR = frontAxleFy * loads.FR / (loads.FL+loads.FR);
    tireFy.RL = rearAxleFy * loads.RL / (loads.RL+loads.RR);
    tireFy.RR = rearAxleFy * loads.RR / (loads.RL+loads.RR);
end

%% Longditudinal Load on Tire
function tireFX = getTireFX(car,loadCase,loads,biasFront,direction) 
    % calculates longd. forces generated by each tires by:
    % > find total longd. force, 
    % Assume outside to inside Fy ratio is proportional to ratio of Fz
    % this should be a conservative assumption, as the derivative of the
    % max possible force decreases 
    
    Ax = loadCase(1);
    totalFx = Ax*car.mass_total*car.g;
    
    % longd. weight transfer on axles
    frontAxleFx = direction*biasFront*totalFx;
    rearAxleFx = direction*(1-biasFront)*totalFx;

    % scale by proportion of normal force on each corner to total axle
    % norm. force.
    tireFX.FL = frontAxleFx * loads.FL / (loads.FL+loads.FR);
    tireFX.FR = frontAxleFx * loads.FR / (loads.FL+loads.FR);
    tireFX.RL = rearAxleFx * loads.RL / (loads.RL+loads.RR);
    tireFX.RR = rearAxleFx * loads.RR / (loads.RL+loads.RR);
    
end








