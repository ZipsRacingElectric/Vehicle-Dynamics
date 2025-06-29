%% BUMP STEER ELIMINATOR
% written by Ben Model and Kenneth Dubos, Oct. 2024
% gib points spreadsheet as per example. Choose front or rear, choose
% inboard or outboard point as optimization point. Can play around with
% values manually. Run. Remove Claude from Design Process, one step at a
% time. Eventually.

% UNIFIED VERSION! 
% now handles front and rear.
% -BMM, 6/25


%% Load Points
filepath = 'C:\Users\benmo\OneDrive - The University of Akron\Documents - Zips Racing FSAE\ZR25\Vehicle Dynamics\System\Analysis\';
pointsSpreadsheet = strcat(filepath,'ZR25_SuspensionForcesHeavy.xlsx');


%% IS FRONT OR REAR????
% rotations weird, make true if front left, false if rear left
% pls no try be clever and make right. I'm not testing that.
isFront = false;
    
    if isFront 
        zTravel = -15:2:45; 
        optPoint = "inboard";       % what I changed for zr25, can also try other way around! tho outboard will make more changes to ackerman
    end
    if not(isFront) 
        zTravel = -20:2:40; 
        optPoint = "outboard";
    end
    % ....su, this works

corner = SuspensionPoints_better(pointsSpreadsheet,isFront);
    

% make changes to inboard or outboard point???
    if optPoint == "inboard"
        initialValue = corner.toeInboard;
    elseif optPoint == "outboard"
        initialValue = corner.toeOutboard;
    else
        warning("are tweaking inboard point or outboard point? no caps!")
    end


% play around with different positions: zero*(tweaks)
corner.toeOutboard = corner.toeOutboard + 0 * [0;-10;48]
corner.toeInboard = corner.toeInboard + 0 * [21;0;0]

%% Simulate each control arms swinging independently
    n = 50000;
    halfAngleDisplacement = 10;         % 20 deg swing, much larger than CAs will swing in actuality
    % vectors to hold trace of outboard points
    lbj_path = getBallJointPath(n,halfAngleDisplacement,corner,true);

    %% want positions of lower ball joint at every half mm from -30 to +30 from start loc.
    % reduces the size of the problem to simulate, but still allows for
    % high precision by searching for matching point on UBJ path from
    % n=50,000 sized array

        indexArrayLBJ = getZDisplacementPointsOnPath(corner,zTravel,lbj_path);

        lbj_path = lbj_path(:,indexArrayLBJ(:));
        ubj_path = getUpperBallJointPath(corner,lbj_path);
        tierod_path = getTieRodPath(corner,lbj_path,ubj_path);
        

%% Calculate Toe angle as function of travel
tic
toeAngleInitial = getToeArray(corner,lbj_path,ubj_path,tierod_path);
toc
figure(1)
plot(-zTravel,toeAngleInitial)
xlabel('Z Travel (mm)')
ylabel('Toe Angle (degrees)')


%% Optimization Stuff!!!!

% This lambda function acts as the objective function to be optimized, and
% handles the fact that the real objective function (toeAngleOptFunc)
% requires extra parameters (corner struct and the ball joint paths), the
% matlab optimization stuff wants it purely as a function on the
% optimization domain (in this case, X defined in R3 space, literally the
% toe outboard point)
    optFunction =@(X) toeAngleOptFunc(corner,optPoint,lbj_path,ubj_path,X);

    if optPoint == "inboard"
        x0 = corner.toeInboard + 1.*[0;10;0];          % initial guess, move around some, idk
    elseif optPoint == "outboard"
        x0 = corner.toeOutboard + 1.*[0;-10;0];
    end

%% optimization problem definition, uses the bulkier struct rather than
% passing arguments directly. Playing with the tolerances determines if it
% converges or errors out, as the binary search thing converges to its own
% defined tolerance
    options = optimoptions("fminunc");
    options.StepTolerance = 1e-7;
    options.OptimalityTolerance = 1e-20;
    options.FunctionTolerance = 1e-20;
    
    problem.options = options;
    problem.x0 = x0;
    problem.objective = optFunction;
    problem.solver = 'fminunc';


%% start optimizing!
tic
[optimizedToePoint,fval,exitflag,output] = fminunc(problem);
toc

howMuchMoved = optimizedToePoint - initialValue        % relative vector to found point from initial pt.
%fval;                        % SSE of toe angle across travel for new point


%% update new point after optimization
if optPoint == "inboard"
    corner.toeInboard = optimizedToePoint;
elseif optPoint == "outboard"
    corner.toeOutboard = optimizedToePoint;
end

tierod_path2 = getTieRodPath(corner,lbj_path,ubj_path);
toeAngleOptimized = getToeArray(corner,lbj_path,ubj_path,tierod_path2);

hold on
plot(-zTravel,toeAngleOptimized)
scatter([10 0 -10 -20 -30 -40],[0 0 0 0 0 0])
title('Bump Steer, Front')
legend('Previous Point', 'New Optimized Point', 'Measured Data (laser line method)')
ylim([-.05 .05])
grid on
hold off

figure(2)
    plot3(lbj_path(1,:),lbj_path(2,:),lbj_path(3,:),ubj_path(1,:),ubj_path(2,:),ubj_path(3,:),tierod_path(1,:),tierod_path(2,:),tierod_path(3,:),tierod_path2(1,:),tierod_path2(2,:),tierod_path2(3,:))
    hold on
    xlim([corner.tireCenter(1) - 200,corner.tireCenter(1) + 200])
    ylim([corner.tireCenter(2) - 200,corner.tireCenter(2) + 200])
    zlim([corner.tireCenter(3) - 200,corner.tireCenter(3) + 200])
    grid on
    hold off


%% functions 
function sse = toeAngleOptFunc(corner,optPoint,lbj_path,ubj_path,X)
    % this is the actual objective function. getToeArray returns an array
    % of angles corresponding to the toe displacement across the range of
    % travel. To distill this to one number, I am taking the SSE of this.
    % To optimize for a specific toe curve, one could take the SSE against
    % a defined function or something.

    % corner struct is used to pass geometrical data around. Trial outboard
    % points are set here, and sent to the getToeArray function
    if optPoint == "inboard"
        corner.toeInboard = X;
    elseif optPoint == "outboard"
        corner.toeOutboard = X;
    end
    tierod_path = getTieRodPath(corner,lbj_path,ubj_path);
    toeAngle = getToeArray(corner,lbj_path,ubj_path,tierod_path);

    sse = sum(toeAngle.^2);
end

function toeAngle = getToeArray(corner,lbj_path,ubj_path,tierod_path)
    % this function calculates the toe curve as fn of sus travel at wheel.
    % It takes two paths, of the lower ball joint at half mm increments (or
    % otherwise defined), and the point on the upper ball joint path
    % corresponding to the lower ball joint's position. (These are
    % calculated in getReducedUBJPathIndex). 
    toeAngle = zeros(1,length(lbj_path));
   

    wheelAxisPoint = corner.tireCenter;
    wheelAxisDir = wheelAxisPoint+[0 20 0]';

    r123 = [norm(corner.lcaOutboard-wheelAxisPoint);norm(corner.ucaOutboard-wheelAxisPoint);norm(corner.toeOutboard-wheelAxisPoint)];
    r456 = [norm(corner.lcaOutboard-wheelAxisDir);norm(corner.ucaOutboard-wheelAxisDir);norm(corner.toeOutboard-wheelAxisDir)];
    trilaterationOptions = optimoptions("fsolve");
    trilaterationOptions.Display = "off";

        for i=1:1:length(lbj_path)
           
            ri = trilateration(lbj_path(:,i),ubj_path(:,i),tierod_path(:,i),r123,wheelAxisPoint,trilaterationOptions);
            rj = trilateration(lbj_path(:,i),ubj_path(:,i),tierod_path(:,i),r456,wheelAxisPoint,trilaterationOptions);
            
            % project onto ground plane, make unit length
            ri = [1 0 0;0 1 0;0 0 0]*ri;
            rj = [1 0 0;0 1 0;0 0 0]*rj;

            rij = (rj-ri)./norm(ri-rj);
            
            % get angle to y axis
            toeAngle(i) = acosd(dot(rij,[0;1;0]));
        end
end

function indexArray = getZDisplacementPointsOnPath(corner,zTravel,path)
% returns an array with indexes matching the location of the closes point
% in the path to the z travel vector

    % z height of LBJ point at static displacement (as defined in geometry)
    refHeight = corner.lcaOutboard(3,1);

    searchArray = path(3,:) - refHeight;
    indexArray = zTravel.*0;
    
    for i = 1:length(zTravel)        
        [m,indexArray(i)] = min(abs(searchArray-zTravel(i)));
    end

end

function bj_path = getBallJointPath(nPoints,thetaSwing,corner,isLCA)
% Returns path of upper or lower ball joint along a +/- thetaSwing swing

    bj_path = zeros(3,nPoints);
    theta = deg2rad(linspace(-thetaSwing,thetaSwing,nPoints));
    
    % lower or upper arm
    if isLCA
        p1 = corner.lcaAft;
        p2 = corner.lcaFore;
        p3 = corner.lcaOutboard;
    else 
        p1 = corner.ucaAft;
        p2 = corner.ucaFore;
        p3 = corner.ucaOutboard;
    end

    % define axis and matrices for rotation matrix
    u_axis = (p2-p1)./norm(p2-p1);
    u_op = getOuterProduct(u_axis,u_axis);  
    u_cpm = getCrossProductMatrix(u_axis);

    % Rotation matrix Lambda Function
    getR =@(theta,u_cpm,u_op) cos(theta)*eye(3) + sin(theta)*u_cpm + (1-cos(theta))*u_op;

    
    for i = 1:1:length(theta)
        % R is rotation matrix
        R = getR(theta(i),u_cpm,u_op);    
        bj_path(:,i) = R*(p3-p1) + p1;        
    end

end

function ubj_path = getUpperBallJointPath(corner,lbj_path)
% Returns path of upper or lower ball joint along a +/- thetaSwing swing
    tolerance = 1e-8;
    ubj_path = zeros(3,length(lbj_path));

    uprightLength = norm(corner.ucaOutboard-corner.lcaOutboard);
    
    % lower or upper arm
    p1 = corner.ucaAft;
    p2 = corner.ucaFore;
    p3 = corner.ucaOutboard;

    % define axis and matrices for rotation matrix
    u_axis = (p2-p1)./norm(p2-p1);
    u_op = getOuterProduct(u_axis,u_axis);  
    u_cpm = getCrossProductMatrix(u_axis);
    
    % find first angle of rotation of UCA corresponding to LCA position
    % with wider search tolerance, so next iterations can be faster
    angle = rotationBinarySearch(tolerance,-15,0,u_axis,p1,uprightLength,lbj_path(:,1),p3);
    R = getFasterRotationMatrix(angle,u_op,u_cpm);    
    ubj_path(:,1) = R*(p3-p1) + p1;


    for i = 1:1:length(lbj_path)
        % find angle of rotation of UCA corresponding to LCA position
        angle = rotationBinarySearch(tolerance,-1,rad2deg(angle),u_axis,p1,uprightLength,lbj_path(:,i),p3);
        
        % R is rotation matrix
        R = getFasterRotationMatrix(angle,u_op,u_cpm);    
        ubj_path(:,i) = R*(p3-p1) + p1;
        %ubj_path(:,i) = R*(p3);
    end

end

function tieRodPath = getTieRodPath(corner,lbjPath,ubjPath)
    tieRodPath = lbjPath.*0;
    referenceLengths = [norm(corner.lcaOutboard-corner.toeOutboard);norm(corner.ucaOutboard-corner.toeOutboard);norm(corner.toeInboard-corner.toeOutboard)];
    
    % passing this to the function speeds it up a lot
    trilaterationOptions = optimoptions("fsolve");
    trilaterationOptions.Display = "off";

    tieRodPath(:,1) = trilateration(lbjPath(:,1),ubjPath(:,1),corner.toeInboard,referenceLengths,corner.toeOutboard,trilaterationOptions);
    
        for i=2:1:length(lbjPath)
            pointA = lbjPath(:,i);
            pointB = ubjPath(:,i);
            pointC = corner.toeInboard;
    
            tieRodPath(:,i) = trilateration(pointA,pointB,pointC,referenceLengths,tieRodPath(:,i-1),trilaterationOptions);
    
        end
end

function p4 = trilateration(p1,p2,p3,r123,x0,trilaterationOptions)
    refLengthsSquared = r123.^2;

    f =@(p_4) p_4'*p_4 - 2*[p1';p2';p3']*p_4 + [p1'*p1;p2'*p2;p3'*p3] -refLengthsSquared;

    p4 = fsolve(f,x0,trilaterationOptions);
end

function angle = rotationBinarySearch(tolerance,initialAngle,angleOffset,axis,origin,refLength,stationaryPoint,pointVector)
%% Binary search to find rotation angle that matches reference length
        % apply a binary search by rotating the point about the given axis
        % changing the angle until it matches to tolerance
        L_angle = deg2rad(-initialAngle+angleOffset);
        R_angle = deg2rad(initialAngle+angleOffset);
        searchAngle = .5*(L_angle+R_angle);

        % Enforce axis to be a unit vector!
        axis = axis./norm(axis);

        % rotation matrix definition
        axis_op = getOuterProduct(axis,axis);
        axis_cpm = getCrossProductMatrix(axis);

        % get points relative to origin
        stationaryPoint = stationaryPoint - origin;
        pointVector = pointVector - origin;
        
        % Determine if positive rotation makes positive or negative length
        % increase
        %a = norm(pointVector-stationaryPoint);
        %rDirection = getRotationLengthDelta(axis_op,axis_cpm,a,stationaryPoint,pointVector,deg2rad(1));
        %rDirection = rDirection/abs(rDirection);
       
        % function to measure length. calculates rotation matrix,
        % applies rotation, measures length, and subtracts toe link length all
        % in one
        lengthDelta = getRotationLengthDelta(axis_op,axis_cpm,refLength,stationaryPoint,pointVector,searchAngle);
    
        while abs(lengthDelta) >= tolerance
            searchAngle = .5*(L_angle+R_angle);
            lengthDelta = getRotationLengthDelta(axis_op,axis_cpm,refLength,stationaryPoint,pointVector,searchAngle);
            
            if lengthDelta >= 0
                % too long
                L_angle = searchAngle; 
            else
                % too short
                R_angle = searchAngle;
            end
        end
        angle = searchAngle;
end

function delta = getRotationLengthDelta(u_op,u_cpm,refLength,stationaryPoint,vector,angle)
    % function to measure length delta. calculates rotation matrix,
    % applies rotation, measures length, and subtracts ref length all
    % in one
    R = getFasterRotationMatrix(angle,u_op,u_cpm);
    delta = norm(R*vector - stationaryPoint) - refLength;
end

function R = getRotationMatrix(theta,u)
    % returns the axis-angle rotation matrix for an angle theta in radians
    % and an axis u
    u = u./norm(u);
    u_cpm = [cross(u,[1 0 0]') cross(u,[0 1 0]') cross(u,[0 0 1]')];
    u_op = u*u';

    R = cos(theta)*eye(3) + sin(theta)*u_cpm + (1-cos(theta))*u_op;

end

function R = getFasterRotationMatrix(theta,u_op,u_cpm)
    R = cos(theta)*eye(3) + sin(theta)*u_cpm + (1-cos(theta))*u_op;
end

function u_cpm = getCrossProductMatrix(u)
    u_cpm = [cross(u,[1 0 0]') cross(u,[0 1 0]') cross(u,[0 0 1]')];
end

function u_op = getOuterProduct(u,v)
% Outer Product Matrix
    u_op = u*v';
end