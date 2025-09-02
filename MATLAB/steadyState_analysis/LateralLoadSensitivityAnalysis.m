% calculates lateral load distribution stuff based on current parameters
% LLD is chosen to give peak lateral grip at a slightly negative stability
% Varies some parameters to determine sensitivities to max lateral grip
% index - 
% reference:


clear

%% load in parameters file
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

    % creates tire object
    fTire = tire(2.2,0.000443656,a_x);
    rTire = tire(2.2,0.000443656,a_x);
    
    % front and rear coefficients. yeah this is a really stupid way to do
    % this. Built for future expandability
    kaf = fTire.ka;
    kbf = fTire.kb;
    
    kar = rTire.ka;
    kbr = rTire.kb;



%% Values I wanna play with
% records the unchanged values here
    mass = zr.mass_total;
    roll_center_rear = zr.roll_center_rear;
    roll_center_front = zr.roll_center_front;
    cgHeight = zr.cg_height;
    frontBias = zr.front_mass_distribution;
    ClA = zr.Cl*zr.frontal_area;
    aeroBias = zr.center_of_pressure_distribution;
    rwd_pct = 1;
    rSA = 5;
    targetLatG = 1.5;
    springRate = zr.spring_rate_rear;



%% Trial setup
% velocity reference value- keep this the same each trial cause I don't 
% want to deal with this problem as a diff eq and have to converge, etc.
velocity = 30;
zr.velocity_skidpad = velocity;

lltArray = linspace(.1,1,200);

   % Control?
    % Calculate total force and Stability Index vs LLT
    [totalForceAvail,stabilityIndex] = sensitivityJawn(zr,fTire,rTire,lltArray,targetLatG,rwd_pct,rSA,velocity);

    controlMaxG = maxPredictedG(zr.mass_total,max(totalForceAvail))
    controlMaxStableG = maxPredictedG(zr.mass_total,getMaxStableLatForce(lltArray,totalForceAvail,stabilityIndex))
    controlMaxStableLTT = getmaxStableLLT(lltArray,stabilityIndex)
    
    figure(1)
    plotLLTvsFy(lltArray,totalForceAvail,stabilityIndex)
    CurrentTLLTD = zr.tlltd

%Trial Value - Start with target value, and run trials from -20% to +20%
% of target value. 
% this turned out to only be a shitty way of taking a partial derivative :)

if true         % turns off all the trials lol. 

nTrials = 101;
trialPercent = 1-linspace(.2,-.2,nTrials);

%% Mass
    originalCarMass = zr.vehicle_mass;
    originalDriverMass = zr.driver_mass;
    trialMass = mass*trialPercent;
    trialMassMaxG = trialMass.*0;
    trialMassMaxGStable = trialPercent.*0;
    trialMassStable = trialPercent.*0;
    for i = 1:1:nTrials
    
        % Set trial mass and update
        zr.vehicle_mass = originalCarMass*trialPercent(i);
        %zr.driver_mass = originalDriverMass*trialPercent(i);

        zr.update;
    
        % Calculate total force and Stability Index vs LLT
        [totalForceAvail,stabilityIndex] = sensitivityJawn(zr,fTire,rTire,lltArray,targetLatG,rwd_pct,rSA,velocity);
    
        maxFy = max(totalForceAvail);
        trialMassMaxG(i) = maxPredictedG(zr.mass_total,maxFy);

        maxFyStable = getMaxStableLatForce(lltArray,totalForceAvail,stabilityIndex);

        trialMassMaxGStable(i) = maxPredictedG(zr.mass_total,maxFyStable);
        trialMassStable(i) = getmaxStableLLT(lltArray,stabilityIndex);
       
    
    end
    
    zr.reset

%% CG
    trialCG = cgHeight*trialPercent;
    trialCGMaxG = trialPercent.*0;
    trialCGMaxGStable = trialPercent.*0;
    trialCGStable = trialPercent.*0;
    for i = 1:1:nTrials
    
        % Set trial mass and update
        zr.cg_height = trialCG(i);
        zr.update;
    
        % Calculate total force and Stability Index vs LLT
        [totalForceAvail,stabilityIndex] = sensitivityJawn(zr,fTire,rTire,lltArray,targetLatG,rwd_pct,rSA,velocity);
    
        maxFy = max(totalForceAvail);
        trialCGMaxG(i) = maxPredictedG(zr.mass_total,maxFy);
        maxFyStable = getMaxStableLatForce(lltArray,totalForceAvail,stabilityIndex);
        trialCGMaxGStable(i) = maxPredictedG(zr.mass_total,maxFyStable);
        trialCGStable(i) = getmaxStableLLT(lltArray,stabilityIndex);
    
    end
    
    zr.reset

%% Weight Distribution
    trialWeightDist = frontBias*trialPercent;
    trialWDMaxG = trialPercent.*0;
    trialWDMaxGStable = trialPercent.*0;
    trialWDStable = trialPercent.*0;

    for i = 1:1:nTrials
    
        % Set trial mass and update
        zr.front_mass_distribution = trialWeightDist(i);
        zr.update;
    
        % Calculate total force and Stability Index vs LLT
        [totalForceAvail,stabilityIndex] = sensitivityJawn(zr,fTire,rTire,lltArray,targetLatG,rwd_pct,rSA,velocity);
    
        maxFy = max(totalForceAvail);
        trialWDMaxG(i) = maxPredictedG(zr.mass_total,maxFy);

        maxFyStable = getMaxStableLatForce(lltArray,totalForceAvail,stabilityIndex);
        trialWDMaxGStable(i) = maxPredictedG(zr.mass_total,maxFyStable);
        trialWDStable(i) = getmaxStableLLT(lltArray,stabilityIndex);
    
    end
    
    zr.reset

%% Cl*A
    trialClA = ClA*trialPercent;
    trialClAMaxG = trialPercent.*0;
    trialClAMaxGStable = trialPercent.*0;
    trialClAStable = trialPercent.*0;

    for i = 1:1:nTrials
    
        % Set trial mass and update
        %zr.vehicle_mass = 184
        zr.Cl = trialClA(i);
        zr.update;
    
        % Calculate total force and Stability Index vs LLT
        [totalForceAvail,stabilityIndex] = sensitivityJawn(zr,fTire,rTire,lltArray,targetLatG,rwd_pct,rSA,velocity);
    
        maxFy = max(totalForceAvail);
        trialClAMaxG(i) = maxPredictedG(zr.mass_total,maxFy);
        maxFyStable = getMaxStableLatForce(lltArray,totalForceAvail,stabilityIndex);

        trialClAMaxGStable(i) = maxPredictedG(zr.mass_total,maxFyStable);
        trialClAStable(i) = getmaxStableLLT(lltArray,stabilityIndex);
    
    end
    
    zr.reset

%% CoP
    trialCoP = aeroBias*trialPercent;
    trialCoPMaxG = trialPercent.*0;
    trialCoPMaxGStable = trialPercent.*0;
    trialCoPStable = trialPercent.*0;
    for i = 1:1:nTrials
    
        % Set trial mass and update
        zr.center_of_pressure_distribution = trialCoP(i);
        zr.update;
    
        % Calculate total force and Stability Index vs LLT
        [totalForceAvail,stabilityIndex] = sensitivityJawn(zr,fTire,rTire,lltArray,targetLatG,rwd_pct,rSA,velocity);
    
        maxFy = max(totalForceAvail);
        trialCoPMaxG(i) = maxPredictedG(zr.mass_total,maxFy);

        maxFyStable = getMaxStableLatForce(lltArray,totalForceAvail,stabilityIndex);

        trialCoPMaxGStable(i) = maxPredictedG(zr.mass_total,maxFyStable);
        trialCoPStable(i) = getmaxStableLLT(lltArray,stabilityIndex);
    
    end
    
    zr.reset

 %% Rear Roll Rate
    trialRearRoll = springRate*trialPercent;
    trialRollMaxG = trialPercent.*0;
    trialRollMaxGStable = trialPercent.*0;
    trialRollStable = trialPercent.*0;
    for i = 1:1:nTrials
    
        % Set trial mass and update
        zr.spring_rate_rear = trialRearRoll(i);
        zr.update;
    
        % Calculate total force and Stability Index vs LLT
        [totalForceAvail,stabilityIndex] = sensitivityJawn(zr,fTire,rTire,lltArray,targetLatG,rwd_pct,rSA,velocity);
    
        maxFy = max(totalForceAvail);
        trialRollMaxG(i) = maxPredictedG(zr.mass_total,maxFy);

        maxFyStable = getMaxStableLatForce(lltArray,totalForceAvail,stabilityIndex);

        trialRollMaxGStable(i) = maxPredictedG(zr.mass_total,maxFyStable);
        trialRollStable(i) = getmaxStableLLT(lltArray,stabilityIndex);
    
    end
    
    zr.reset



%% Plot Sensitivities
% Max LatG
figure(2)
plot(trialPercent,trialMassMaxG)
hold on
plot(trialPercent,trialCGMaxG)

plot(trialPercent,trialWDMaxG)
plot(trialPercent,trialClAMaxG)
plot(trialPercent,trialCoPMaxG)

grid on
legend('Mass','CG Height','Weight Distribution','ClA','CoP')
xlabel('% change from target')
ylabel('Max Lateral G')
title('Lateral G Sensitivities (Steeper is more Sensitive)')
hold off

massSensitivity = mean(diff(trialMassMaxG))
cgSensitivity = mean(diff(trialCGMaxG))
aeroSensitivity = mean(diff(trialClAMaxG))

massSensitivity/ aeroSensitivity




% Stable Lat G
figure(3)
subplot(2,1,1)
hold on
%plot(trialPercent,trialWDMaxGStable)
%plot(trialPercent,trialCoPMaxGStable)

plot(trialPercent,trialMassMaxGStable)

plot(trialPercent,trialCGMaxGStable)

plot(trialPercent,trialWDMaxGStable)
plot(trialPercent,trialClAMaxGStable)
plot(trialPercent,trialCoPMaxGStable)
grid on

legend('Weight Distribution','CoP')
xlabel('% change from target')
ylabel('Max Lateral G')
title('Stable Lateral G Sensitivities (Steeper is more Sensitive)')
hold off


subplot(2,1,2)
hold on
plot(trialPercent,trialWDStable)
plot(trialPercent,trialCoPStable)

legend('Weight Distribution','CoP')
xlabel('% change from target')
ylabel('LLT')
title('Maximum Stable LLT')
grid on
hold off


end


%% Functions to hide mess
function [FyAvail,stbi] = sensitivityJawn(zr,fTire,rTire,lltArray,targetLatG,rwd_pct,rSA,velocity)

    [downforce,~] = zr.calculate_aero_forces(velocity);

    % roll axis to CG height at CG
    rcm = zr.roll_center_front + (1-zr.front_mass_distribution)*(zr.roll_center_rear-zr.roll_center_front);
    height_sprung_mass = zr.cg_height-rcm;
    
    % Lateral Load Transfer for current trial
    prm = zr.track_width_front*zr.total_roll_rate_front/(zr.track_width_front*zr.total_roll_rate_front+zr.track_width_rear*zr.total_roll_rate_rear);
    
    % Static Load- Fz on tires
    totalLatForce = targetLatG*zr.mass_total * zr.g;   % in N

    fTire.staticLoad = (zr.front_mass*zr.g + downforce*zr.center_of_pressure_distribution)/2;
    rTire.staticLoad = (zr.rear_mass*zr.g + downforce*(1-zr.center_of_pressure_distribution))/2;

   
    % Max Total Lateral Force
    
    FyAvail = calcTotalForceAvailable(zr,fTire,rTire,lltArray,rwd_pct,totalLatForce); 
    
    % Stability index:
    % uses bycicle model 
    % horribleness hidden in function
    stbi = getStabilityIndex(zr,fTire,rTire,lltArray,rSA,rwd_pct,totalLatForce); 

end
function deez = plotLLTvsFy(lltArray,totalForceAvail,stabilityIndex)
% plot of Total force and Stability Index vs Lat. Load Transfer %

        subplot(2,1,1)
        plot(lltArray,totalForceAvail)
        grid on
        xlabel('LLT')
        ylabel('Max Lateral Force')
        title('Lateral Load Transfer vs Max Lateral Force and Stability Index')
        
        subplot(2,1,2)
        plot(lltArray,stabilityIndex,lltArray,-.025.*ones(size(lltArray)))
        grid on
        xlabel('LLT')
        ylabel('Stability Index (sans Fx Component)')
        legend('Stability Index as fn of LLTD','Critically Stable Point (-.025)')
end

function maxLLT = getMaxLLT(LLT,FyAvail)
    maxFy = max(FyAvail);
    maxLLT = LLT(find(FyAvail==maxFy));
end


function stableMaxFy = getMaxStableLatForce(llt,FyAvail,stability)
    % returns lateral force available at critically stable point on LLTD vs
    % stb. index (-.025)
    excludeIndex = length(stability) - length(stability)*.25;
    stability = stability(1:excludeIndex);
    i=find(min(abs(stability-(-.025)))==abs(stability-(-.025)));
    stableMaxFy = FyAvail(i);
end

function stableLLT = getmaxStableLLT(LLT,stability)
excludeIndex = length(stability) - length(stability)*.25;
stability = stability(1:excludeIndex);
    i=find(min(abs(stability-(-.025)))==abs(stability-(-.025)));  
    stableLLT = LLT(i);
end

function maxG = maxPredictedG(m,Fy) 
    maxG = Fy/m/9.81;
end

function totalForceAvail = calcTotalForceAvailable(zr,fTire,rTire,lltArray,rwd_pct,totalLatForce)
    totalForceAvail = fTire.staticForce + rwd_pct*rTire.staticForce - 2*(lltArray.^2)*fTire.ka*fTire.kb*(totalLatForce*zr.cg_track_ratio)^2 ;
    totalForceAvail = totalForceAvail - 2*((1-lltArray).^2)*rwd_pct*rTire.ka*rTire.kb*(totalLatForce*zr.cg_track_ratio)^2;
end


function [frontSlipSlope,rearSlipSlope] = calcSlipSlopes(zr,fTire,rTire,lltArray,rSA,rwd_pct,totalLatForce)
    % these are from chassis guy, idk man
    % Calculate front slip angle
    fSA = rSA*(zr.b*(rwd_pct*rTire.staticForce  - 2*((1-lltArray).^2).*rwd_pct*rTire.ka*rTire.kb .* (totalLatForce*zr.cg_track_ratio).^2));
    fSA = fSA./ (zr.a*(fTire.staticForce - 2*(lltArray.^2).*fTire.ka*fTire.kb.*(totalLatForce*zr.cg_track_ratio).^2));
   
    % Slope of slip angles
    frontSlipSlope = fTire.a_1 + 2*fTire.a_2.*(fSA./57.3) + 3*fTire.a_3.*(fSA./57.3).^2 + 4*fTire.a_4.*(fSA./57.3).^3;
    % not sure what these constants do
    rearSlipSlope = 14.323 -1306.2.*(rSA./57.3).^2 * ones(1,length(lltArray));  
end

function [frontForce,rearForce] = calcFrontRearForces(zr,fTire,rTire,lltArray,totalLatForce,rwd_pct)
    % force on front and rear tires
    frontForce = (fTire.staticForce - 2*(lltArray.^2).*fTire.ka*fTire.kb .* (totalLatForce*zr.cg_track_ratio).^2);
    rearForce = rwd_pct*rTire.staticForce - 2*((1-lltArray).^2) * rwd_pct * rTire.ka*rTire.kb .* (totalLatForce*zr.cg_track_ratio).^2;
end

function stbi = getStabilityIndex(zr,fTire,rTire,lltArray,rSA,rwd_pct,totalLatForce)
    % neglects term for horizontal force

    [frontSlipSlope,rearSlipSlope] = calcSlipSlopes(zr,fTire,rTire,lltArray,rSA,rwd_pct,totalLatForce);
    [frontForce,rearForce] = calcFrontRearForces(zr,fTire,rTire,lltArray,totalLatForce,rwd_pct);

    stbi =(zr.a*frontSlipSlope.*frontForce - zr.b*rearSlipSlope.*rearForce) ./ (zr.wheelbase*(frontSlipSlope.*frontForce + rearSlipSlope.*rearForce ) );

end



