function output = send_can_message(data)
    persistent py_process py_stdin py_stdout

    java.lang.Thread.sleep(10);
    
    CAN_DATA_SCALING_FACTOR = 1000;

    PYTHON_CMD_START = sprintf('./venv/bin/python can_service.py');

    output = 0.0;

    % Start Python process only once
    if isempty(py_process) || ~py_process.isAlive()
        fprintf('Starting persistent Python CAN process...\n');
        py_process = java.lang.Runtime.getRuntime().exec(PYTHON_CMD_START);
        py_stdin = py_process.getOutputStream();
        py_stdout = java.io.BufferedReader( ...
            java.io.InputStreamReader(py_process.getInputStream()));
        
        % Wait for the "READY" signal
        ready_line = char(py_stdout.readLine());
        if ~strcmp(ready_line, "READY")
            %fprintf(ready_line);
            error('CAN Python process did not start correctly.');
        end
    end

    % Format command
    if(strcmp(data,'exit'))
        data_command = sprintf('exit\n');
        py_stdin.write(uint8(data_command));  % Convert to bytes
        py_stdin.flush();
        fprintf('Closing can_service.py \n');
    else
        input_data = data(1);
        plant_data = data(2);
        time_step = data(3);
        data_command = sprintf('%d %d %d\n', input_data, plant_data, time_step);
        fprintf('Command sent was %s', data_command);
        py_stdin.write(uint8(data_command));  % Convert to bytes
        py_stdin.flush();
    
        % Read response
        response = char(py_stdout.readLine());
        fprintf('Response was %s\n', response);
        actuating_signal = str2double(strtrim(response));
    
        fprintf('CAN Response: %s (scaled: %f)\n', response, actuating_signal);
    
        % Check for bad result
        if isnan(actuating_signal)
            fprintf('Bad response from DUT. Stopping Simulation.\n');
            % Close the simulation
            %set_param(bdroot, 'SimulationCommand', 'stop');
        end

        output = double(actuating_signal);
    end
end