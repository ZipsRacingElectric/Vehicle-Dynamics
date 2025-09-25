% toe compliance data processer feed in forces and deflections get back
% angles num num. 

clear, clc

Excel = readtable("C:\Users\ATuck\OneDrive - The University of Akron\Documents\GitHub\Vehicle-Dynamics\MATLAB\compliance_models\Compliance Testing ZR25.xlsx")
Distancem = 1


% Positive applied toe
DistancePY = Excel{3,6}
DistancePY = str2double(DistancePY)*.0254
AppliedForcesP = Excel{1:4,1}
DisplacementDial1P = Excel{1:4,2}.*0.01.*.0254
DisplacementDial2P = Excel{1:4,3}.*0.01.*.0254
DisplacementDial3P = Excel{1:4,4}.*0.01.*.0254
DisplacementDial4P = Excel{1:4,5}.*0.01.*.0254

PositiveMoments = GetMoment(AppliedForcesP,Distancem)
PDC = GetDisplacements(DisplacementDial1P, DisplacementDial2P) % displacement at corner
PDA = GetDisplacements(DisplacementDial3P, DisplacementDial4P) % displacement at adjacent corner
PDC = PDC % convert to meters
PDA = PDA % Convert to meters
Anglea = GetCompAng(DistancePY, PDC) % corner side theta calc
Angle1 = rad2deg(Anglea)
Angleb = GetCompAng(DistancePY, PDA) % adjacent side theta calc
Angle2 = rad2deg(Angleb)
Compliance = GetToeCompliance(Angle1, PositiveMoments) % compliance at corner

% Plot
plot(Compliance)


% negative applied toe
%DistanceNY = Excel{5,6}
%AppliedForcesN = Excel{9:12,1}
%DisplacementDial1N = Excel{9:12,2}
%DisplacementDial2N = Excel{9:12,3}
%DisplacementDial3N = Excel{9:12,4}
%DisplacementDial4N = Excel{9:12,5}




%% Functions


function moments = GetMoment(FKg, Dm) % get moments back in Nm

    moments = FKg*Dm % in Newtons?
end

function r = GetDisplacements(r1in, r2in)

r = r2in-r1in

end


function theta = GetCompAng(Dispx,Dispy)

c = Dispy./Dispx


theta = atan(c)

end

function TC = GetToeCompliance(ang,Moment) % outputs degrees per Kg

    TC = ang./Moment

end