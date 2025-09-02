classdef tire < handle
    %TIRE Object that contains some stuff related to tires
    %   This is pretty old and I didn't end up using it a ton, but the
    %   structure is there for you guys to adapt. Super basic rn, kinda
    %   overkill
    %   Currently I think only the LLT stuff uses it
    %
    %   create new one with:
    %   (Basic Tire Parameters, from chassissim guy)
    %   a_x = [12.476 82.761 -1592.7 4985.8];
    %
    %   front = tire(2.2,0.000443656,a_x);
    %
    %   Access variables like:
    %   mu = front.ka

    properties
        % Initial Coeff. of Friction
        ka
       
        % Load Sensitivity
        kb

        % ChassisSim Guy Coefficients
        a_1
        a_2
        a_3
        a_4

        staticLoad      % normal load
        staticForce     % lat force

        g = 9.81;
    end

    methods
        function obj = tire(a,b,coeffs)
            %UNTITLED Construct an instance of this class
            %   Constructor function - pass parameters file, tires object,
    
            obj.ka = a;
            obj.kb = b;
            obj.a_1 = coeffs(1);
            obj.a_2 = coeffs(2);
            obj.a_3 = coeffs(3);
            obj.a_4 = coeffs(4);
            

            

        end

        %% Tire functions
        function obj = set.staticLoad(obj,staticLoad)
            % in kgf???
            % dumbass moment

            % 8/27/25 - nope, mandating newtons out of principle
            obj.staticLoad = staticLoad;
        end
        function staticLatForce = get.staticForce(obj)
            % returns in Newtons
            % ^ written by a dumbass
            staticLatForce = 2*obj.ka*(1 - obj.kb*obj.staticLoad)*obj.staticLoad;

        end

    end
end