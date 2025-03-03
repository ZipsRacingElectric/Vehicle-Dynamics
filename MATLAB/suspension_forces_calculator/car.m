classdef car < handle
    %CAR Object that contains all simulation parameters of the car
    %   Use this class to instanciate a CAR object that you will
    %   assign relevant parameters. Use the CAR object to recall
    %   these parameters easily in your simulation
    %
    %   create new one with:
    %   ZR25 = car('parameters file path and name')
    %
    %   Access variables like:
    %   mass = ZR25.mass


    properties (Access=private)
        % no touchy!
        parameters
    end

    properties (Access=public)
        %% Other Objects (not really implemented)
        tire
        motor

        %% non-calculated- pulled from excel when object is created

        % mass        
        carMass
        frontCornerMassUS
        rearCornerMassUS
        driverMass
        
        weightBias
        
        % size
        wheelbase
        trackFront
        trackRear
        wheelRadius
        cgZ
        
        damperTravel
        
        % Suspension/Kinematics
        fTireSpringRate
        rTireSpringRate
        
        fMotionRatio
        rMotionRatio
        fSpringRate
        rSpringRate
        
        fBarMotionRatio
        rBarMotionRatio
        fBarRate
        rBarRate
        
        fRollCenter
        rRollCenter
        
        % Aero
        vSkid
        claVskid
        copVskid
        
        vMax
        claVmax
        copVmax
       
        % these are sorta references
        fRideFreq
        rRideFreq

        %% Calculated - call update method to calculate
        % mass stuff
        mass
        frontMassUS
        rearMassUS
        unsprungMass
        
        sprungMass      
        frontMassS
        rearMassS
        cornerMassAvg
        
        % size/dimensions
        a
        b
        bSprung
        aSprung
        cgZSprung
        meanTrack
        cgTrackRatio
        
        % Axle/Corner Loads (in kg!!!)
        frontWeight
        rearWeight
        frontCornerWeight
        rearCornerWeight
        
       
        % Rates
        fWheelRate
        rWheelRate
        
        fBarWheelRate
        rBarWheelRate
        fTotalSpringRate
        rTotalSpringRate
        fRollWheelRate
        rRollWheelRate
        fTotalRollRate
        rTotalRollRate
        
        % Aero     
        downforceVskidFront
        downforceVskidRear
        downforceVmaxFront
        downforceVmaxRear

        % Lateral Load Transfer Distribution Stuff
        heightNRAtoSMCenter
        Kf
        Kr
        KfPrime
        KrPrime
        fLLTD
        rLLTD
        TLLTD

        % Longditudinal Load Transfer Stuff
        LongLoadTransfer



        %% Constants
        % Environment
        rhoAir = 1.225;
        rhoAirHot = 1.164;      % at 30 degrees C
        g = 9.81;
    end

    methods
        function obj = car(parametersFile)
            %UNTITLED Construct an instance of this class
            %   Constructor function - pass parameters file, tires object,
            %   
            parameterTable = readtable(parametersFile,'Sheet','parameters');
            
            % stores unchanged parameters from Excel file for reset
            % private variable, no touchy!
            obj.parameters = containers.Map(parameterTable.Var1, parameterTable.Var2);
            
            % set non-calculated variables in class to values from
            % parameters container
            obj.setParameters();            

            % evaluate calculated parameters
            obj.update();          

        end
        function r = setParameters(obj)
            % sets parameters of object to parameters from container read
            % from excel.
            % slightly confusing way but this is called in object
            % constructor and by obj.reset() method

            %% Set Non-Calculated Parameters
            % dont worry, I made the computer write this

            % Mass and Bias stuff
            obj.carMass = obj.parameters('carMass');
            obj.frontCornerMassUS = obj.parameters('frontCornerMassUS');
            obj.rearCornerMassUS = obj.parameters('rearCornerMassUS');
            obj.driverMass = obj.parameters('driverMass');
            
            obj.weightBias = obj.parameters('weightBias');
            
            % Size
            obj.wheelbase = obj.parameters('wheelbase');
            obj.trackFront = obj.parameters('trackFront');
            obj.trackRear = obj.parameters('trackRear');
            obj.wheelRadius = obj.parameters('wheelRadius');
            obj.cgZ = obj.parameters('cgZ');
            
            % Suspension stuff
            obj.damperTravel = obj.parameters('damperTravel');

            obj.fTireSpringRate = obj.parameters('fTireSpringRate');
            obj.rTireSpringRate = obj.parameters('rTireSpringRate');
            
            obj.fMotionRatio = obj.parameters('fMotionRatio');
            obj.rMotionRatio = obj.parameters('rMotionRatio');
            obj.fSpringRate = obj.parameters('fSpringRate');
            obj.rSpringRate = obj.parameters('rSpringRate');
            
            obj.fBarMotionRatio = obj.parameters('fBarMotionRatio');
            obj.rBarMotionRatio = obj.parameters('rBarMotionRatio');
            obj.fBarRate = obj.parameters('fBarRate');
            obj.rBarRate = obj.parameters('rBarRate');
            
            obj.fRollCenter = obj.parameters('fRollCenter');
            obj.rRollCenter = obj.parameters('rRollCenter');
            
            % Aero stuff
            obj.vSkid = obj.parameters('vSkid');
            obj.claVskid = obj.parameters('claVskid');
            obj.copVskid = obj.parameters('copVskid');
            
            obj.vMax = obj.parameters('vMax');
            obj.claVmax = obj.parameters('claVmax');
            obj.copVmax = obj.parameters('copVmax');
            
            % These really should be calculated
            obj.fRideFreq = obj.parameters('fRideFreq');
            obj.rRideFreq = obj.parameters('rRideFreq');

        end
        function r = update(obj)
            % Call this if you change a property that calculations rely on
            % otherwise, dependent properties will be stale (aka wrong)

            % yes, a properly implemented class would protect those
            % properties but oh well. Might do that later.

            %% Evaluates calculated parameters
            % I asked the computer nicely to write these as well

            % unsprung mass
            obj.mass = obj.carMass+obj.driverMass;
            obj.frontMassUS = 2*obj.frontCornerMassUS;
            obj.rearMassUS = 2*obj.rearCornerMassUS;
            obj.unsprungMass = obj.frontMassUS+obj.rearMassUS;
            obj.sprungMass = obj.carMass+obj.driverMass-obj.unsprungMass;
            
            % Distribution
            obj.a = obj.wheelbase*(1-obj.weightBias);
            obj.b = obj.wheelbase*(obj.weightBias);
            
            % Sprung Mass
            obj.frontMassS = (obj.mass*obj.b/obj.wheelbase)-obj.frontMassUS;
            obj.rearMassS = obj.sprungMass-obj.frontMassS;
            obj.cornerMassAvg = obj.unsprungMass/4;
            
            % Size            
            obj.bSprung = (obj.mass*obj.b-obj.frontMassUS*obj.wheelbase)/obj.sprungMass;
            obj.aSprung = obj.wheelbase-obj.bSprung;
            obj.cgZSprung = ((obj.mass*obj.cgZ)-(obj.frontMassUS*obj.wheelRadius)-(obj.rearMassUS*obj.wheelRadius))/obj.sprungMass;
            obj.meanTrack = mean([obj.trackFront, obj.trackRear]);
            obj.cgTrackRatio = obj.cgZ/obj.meanTrack;
            
            % wheel loads (in kg!!!!!!!!)
            obj.frontWeight = obj.mass*obj.weightBias;
            obj.rearWeight = obj.mass*(1-obj.weightBias);
            obj.frontCornerWeight = obj.frontWeight/2;
            obj.rearCornerWeight = obj.rearWeight/2;
            
            % Spring rates and stuff
            obj.fWheelRate = obj.fSpringRate*obj.fMotionRatio^2;
            obj.rWheelRate = obj.rSpringRate*obj.rMotionRatio^2;
            
            obj.fBarWheelRate = obj.fBarRate*obj.fBarMotionRatio^2;
            obj.rBarWheelRate = obj.rBarRate*obj.rBarMotionRatio^2;
            obj.fTotalSpringRate = (obj.fWheelRate*obj.fTireSpringRate)/(obj.fWheelRate+obj.fTireSpringRate);
            obj.rTotalSpringRate = (obj.rWheelRate*obj.rTireSpringRate)/(obj.rWheelRate+obj.rTireSpringRate);
            obj.fRollWheelRate = obj.fWheelRate+obj.fBarWheelRate;
            obj.rRollWheelRate = obj.rWheelRate+obj.rBarWheelRate;
            obj.fTotalRollRate = (obj.fRollWheelRate*obj.fTireSpringRate)/(obj.fRollWheelRate+obj.fTireSpringRate);
            obj.rTotalRollRate = (obj.rRollWheelRate*obj.rTireSpringRate)/(obj.rRollWheelRate+obj.rTireSpringRate);
            
            % Aero stuff 
            obj.downforceVskidFront = (obj.claVskid*obj.rhoAir*0.5*obj.vSkid^2)*obj.copVskid/obj.g;
            obj.downforceVskidRear = (obj.claVskid*obj.rhoAir*0.5*obj.vSkid^2)*(1-obj.copVskid)/obj.g;
            obj.downforceVmaxFront = (obj.claVmax*obj.rhoAir*0.5*obj.vMax^2)*obj.copVmax/obj.g;
            obj.downforceVmaxRear = (obj.claVmax*obj.rhoAir*0.5*obj.vMax^2)*(1-obj.copVmax)/obj.g;

            % Lateral Load Transfer Stuff
            obj.heightNRAtoSMCenter = abs((obj.rRollCenter-obj.fRollCenter)*(obj.a-obj.aSprung) - (-obj.b-obj.a)*obj.cgZSprung + (-obj.b)*obj.fRollCenter - obj.rRollCenter*obj.a) / sqrt((obj.fRollCenter-obj.rRollCenter)^2 + (-obj.b-obj.a)^2);
            obj.Kf = 0.5*obj.fTotalRollRate*obj.trackFront^2;
            obj.Kr = 0.5*obj.rTotalRollRate*obj.trackRear^2;
            obj.KfPrime = obj.Kf-(obj.wheelbase-obj.aSprung)*obj.sprungMass*obj.g*obj.heightNRAtoSMCenter/obj.wheelbase;
            obj.KrPrime = obj.Kr-obj.aSprung*obj.sprungMass*obj.heightNRAtoSMCenter/obj.wheelbase;
            obj.fLLTD = obj.getFLatLoadTransferDist;
            obj.rLLTD = obj.getRLatLoadTransferDist;
            obj.TLLTD = obj.fLLTD / (obj.fLLTD+obj.rLLTD);

            % Longditudinal Load Transfer Stuff (increase in force per G)
            obj.LongLoadTransfer = obj.cgZ/obj.wheelbase*obj.mass*obj.g;
        end

        function r = reset(obj)
            % sets non-calculated objects back to default from excel file
            obj.setParameters()
            obj.update()
        end
        
        function fLLTD = getFLatLoadTransferDist(obj)
            l = obj.wheelbase;
            frc = obj.fRollCenter;
  
            % assume same as cg height for now
            fSMHeight = obj.cgZSprung;
            fLLTD = (obj.sprungMass*obj.g/obj.trackFront)*...
                ((obj.heightNRAtoSMCenter*obj.KfPrime)/(obj.Kf+obj.Kr-obj.sprungMass*obj.g*obj.heightNRAtoSMCenter)+(l-obj.aSprung)/l*frc) +...
                obj.frontMassUS*obj.g/obj.trackFront*fSMHeight;
        end
        
        function rLLTD = getRLatLoadTransferDist(obj)
            l = obj.wheelbase;
            rrc = obj.rRollCenter;

  
            % assume same as cg height for now
            rSMHeight = obj.cgZSprung;
            rLLTD = (obj.sprungMass*obj.g/obj.trackRear)*...
                ((obj.heightNRAtoSMCenter*obj.KrPrime)/(obj.Kf+obj.Kr-obj.sprungMass*obj.g*obj.heightNRAtoSMCenter)+(obj.aSprung)/l*rrc) +...
                obj.rearMassUS*obj.g/obj.trackRear*rSMHeight;
        end


        %% Tire functions


        
        %% Aero Functions
        function F_d = getDrag(obj,Vx)
            %F_d = 0.5 * obj.rhoAir * Vx.^2 * obj.Cd * obj.frontalArea;
        end

        function F_L = getDownForce(obj,Vx)
            F_L = 0.5 * obj.rhoAir * Vx.^2 * obj.claVmax;
        end

    end
end