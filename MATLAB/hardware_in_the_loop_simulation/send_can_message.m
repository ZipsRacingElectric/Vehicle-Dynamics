function output = send_can_message(data_vector)
    input_data = data_vector(1);
    plant_data = data_vector(2);
    time_step = data_vector(3);
    output = 0.0;
    
    % Construct the command string.
    % Make sure to use the correct path to your Python interpreter in the virtual env.
    command = sprintf('./venv/bin/python can_module.py %d %d %d', input_data, plant_data, time_step);
    
    % Call the Python script using system()
    [status, cmdout] = system(command);
    
    % Check if the command executed successfully.
    if status == 0
        % Assume the final line of output contains the actuating signal.
        % Split the output by newline and take the last non-empty line.
        lines = strsplit(strtrim(cmdout), '\n');
        actuating_signal_str = lines{end};
        actuating_signal = str2double(actuating_signal_str);
        disp(cmdout);
        fprintf('Actuating Signal: %f\n', actuating_signal);
    else
        disp('Error running Python script:');
        disp(cmdout);
    end
    
    % If the actuating signal is invalid, stop the simulation.
    if isnan(actuating_signal)
        fprintf('Bad response from DUT. Stopping Simulation.\n');
        actuating_signal = 2;
        %set_param(bdroot, 'SimulationCommand', 'stop');
    end
    
    % Assign the actuating signal to the output.
    output = double(actuating_signal);
end