classdef SuspensionPoints_better < handle
    % Suspension points object, contains front or rear sus points.
    %   Sure, objects are not the most efficient data structure,
    %   but whatever.

    % All this does is read points, and mirror for left/right.
    % All math is done in SuspensionForces 

    %
    %   create new one with:
    %   frontSuspension = SuspensionPoints('.zsf file path and name')
    %
    %   Access variables like:obj.
    %   frontUCAFore = frontSuspension.ucaFore


    properties (Access=public)
        %% 
        ucaFore
        ucaAft
        lcaFore
        lcaAft

        ucaOutboard
        lcaOutboard
        
        toeOutboard
        toeInboard

        pushOutboard
        pushInboard
       
        contactPatch
        tireCenter
    end

    methods
        function obj = SuspensionPoints_better(pointsFile,isFront)
            %UNTITLED Construct an instance of this class
            %   Constructor function - 

            pointsTable = readtable(pointsFile,'Sheet','points');
            pointsTable = containers.Map(pointsTable.Var1, pointsTable.Var2);

            if isFront
                frontRear = 'F_';
            else
                frontRear = 'R_';
            end
                
            % UCA
            link = 'UCA_';
                loc = ["FORE_";"AFT_";"OUT_"];
                obj.ucaFore = getPoint(obj,pointsTable,frontRear,link,loc(1));    
                obj.ucaAft = getPoint(obj,pointsTable,frontRear,link,loc(2));    
                obj.ucaOutboard = getPoint(obj,pointsTable,frontRear,link,loc(3));
            % LCA
            link = 'LCA_';
                obj.lcaFore = getPoint(obj,pointsTable,frontRear,link,loc(1));
                obj.lcaAft = getPoint(obj,pointsTable,frontRear,link,loc(2));    
                obj.lcaOutboard = getPoint(obj,pointsTable,frontRear,link,loc(3));            
            % TIE
            link = 'TIE_';
                loc = ["IN_";"OUT_"];
                obj.toeInboard = getPoint(obj,pointsTable,frontRear,link,loc(1));
                obj.toeOutboard = getPoint(obj,pointsTable,frontRear,link,loc(2)); 
            % PUSH/PULL
            link = 'PUSH_';
                obj.pushInboard = getPoint(obj,pointsTable,frontRear,link,loc(1));
                obj.pushOutboard = getPoint(obj,pointsTable,frontRear,link,loc(2)); 
            % CONTACT PATCH
            link = 'TIRE_';
                loc = 'CP_';
                obj.contactPatch = getPoint(obj,pointsTable,frontRear,link,loc);
                
            % TIRE CENTER
                loc = 'CENTER_';
                obj.tireCenter = getPoint(obj,pointsTable,frontRear,link,loc);
        end

        function point = getPoint(obj,pts,fr,link,loc)
            ptStr = strcat(fr,link,loc);
            point = [pts(strcat(ptStr,'X'));pts(strcat(ptStr,'Y'));pts(strcat(ptStr,'Z'))];
        end

        function r = mirror(obj)
            flipY = [1,-1,1]';

            obj.ucaFore = obj.ucaFore .* flipY;
            obj.ucaAft = obj.ucaAft .* flipY;
            obj.lcaFore = obj.lcaFore .* flipY;
            obj.lcaAft = obj.lcaAft .* flipY;
    
            obj.ucaOutboard = obj.ucaOutboard .* flipY;
            obj.lcaOutboard = obj.lcaOutboard .* flipY;
            
            obj.toeOutboard = obj.toeOutboard .* flipY;
            obj.toeInboard = obj.toeInboard .* flipY;
    
            obj.pushOutboard = obj.pushOutboard .* flipY;
            obj.pushInboard = obj.pushInboard .* flipY;
           
            obj.contactPatch = obj.contactPatch .* flipY;   
            obj.tireCenter = obj.tireCenter .* flipY; 
            % lmao.
            r =0;

        end    
    end
end


