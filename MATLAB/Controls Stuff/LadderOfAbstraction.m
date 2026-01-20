%% LADDER OF ABSTRACTION — FSAE VEHICLE MODELING
% Auto-incrementing abstraction ladder with step-up / step-down tracking

clear; clc;

%% -------------------------------
% INITIALIZE LADDER
%% -------------------------------
ladder = struct( ...
    'level', {}, ...
    'direction', {}, ...   % "up" or "down"
    'name', {}, ...
    'assumptions', {}, ...
    'features', {}, ...
    'reason', {}, ...      % why this rung was added
    'timestamp', {}, ...
    'files', {} );

%% -------------------------------
% ADD RUNGS (AUTO LEVEL)
%% -------------------------------

ladder = addRung( ...
    ladder, ...
    "Steady-State Bicycle Model", ...
    "up", ...
    { ...
        "Constant speed", ...
        "Linear tire stiffness", ...
        "No suspension dynamics" ...
    }, ...
    { ...
        "Understeer gradient" ...
    }, ...
    "Baseline understanding of vehicle behavior", ...
    { "ControlsBicycleModel.m" } ...
);


ladder = addRung(ladder, "nonlinear tire model", 'up', 'no lateral load transfer', 'works better in nonlinear range', 'works better', 'file');
ladder = addRung(ladder, "nonlinear tire model", 'up', 'no lateral load transfer', 'works better in nonlinear range', 'works better', 'file');


%% -------------------------------
% DISPLAY LADDER SUMMARY
%% -------------------------------
fprintf('\n--- LADDER OF ABSTRACTION LOG ---\n');
for i = 1:length(ladder)
    fprintf('\nStep %d | Level %d | %s\n', i, ladder(i).level, ladder(i).direction);
    fprintf('  Model: %s\n', ladder(i).name);
    fprintf('  Reason: %s\n', ladder(i).reason);
end

%% -------------------------------
% VISUALIZE STEP UP / DOWN
%% -------------------------------
levels = [ladder.level];
times  = [ladder.timestamp];

figure
stairs(times, levels, 'LineWidth', 2)
hold on

% Mark step-downs
downIdx = find(strcmp({ladder.direction}, "down"));
plot(times(downIdx), levels(downIdx), 'ro', 'MarkerSize', 8)

xlabel('Time')
ylabel('Model Fidelity Level')
title('Vehicle Model Ladder with Step-Up and Step-Down')
legend('Abstraction Level', 'Step Down', 'Location', 'best')
grid on

%% Local function

function ladder = addRung(ladder, name, direction, assumptions, features, reason, files)

    if isempty(ladder)
        currentLevel = 0;
    else
        currentLevel = ladder(end).level;
    end

    if strcmp(direction, "up")
        newLevel = currentLevel + 1;
    elseif strcmp(direction, "down")
        newLevel = max(currentLevel - 1, 1);
    else
        error('Direction must be "up" or "down"');
    end

    rung.level = newLevel;
    rung.direction = direction;
    rung.name = name;
    rung.assumptions = assumptions;
    rung.features = features;
    rung.reason = reason;
    rung.timestamp = datetime("now");
    rung.files = files;

    ladder(end+1) = rung;
end
