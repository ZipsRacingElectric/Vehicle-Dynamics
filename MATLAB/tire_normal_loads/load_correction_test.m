clc; clear all;

fz = [62.5, -12.5, 150, -10];
total_normal_load = sum(fz);

% set tires not in contact w/ ground to 0 N
fz(fz < 0) = 0;

% find the tires that are on the ground
grounded_tires = fz > 0;

% Redistribution is needed only if a tire is off the ground
if sum(grounded_tires) < 4

    % correct the other tire normal loads for the contact adjustment
    excess_load = sum(fz) - total_normal_load;
    
    % We will subtract load from tires that are contacting the ground based
    % on the ratios of their normal loads
    fz = fz - (excess_load .* (fz ./ (sum(fz(grounded_tires)))));

    correct_normal_load = (total_normal_load == sum(fz));

    % check for car tipping laterally edge case
    if grounded_tires(1) && grounded_tired(3) || grounded_tires(2) && grounded_tires(4)
        % car is up on two wheels laterally, and the normal loads should
        % reflect the car's static weight distribution & long load transfer
    end

    % check for car tipping longitudinally edge case
    if grounded_tires

    end
end