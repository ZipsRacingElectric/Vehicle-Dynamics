#!/usr/bin/env python

import sys
import platform
import can
import subprocess
import time
import math

def ensure_can_interface(interface='can0', bitrate=1000000):
    """
    Checks if the CAN interface is available and up.
    On Linux, uses 'ip' command.
    On macOS (Darwin), uses 'ifconfig'.
    Returns True if the interface is detected and active, otherwise False.
    """
    system = platform.system()
    if system == 'Linux':
        try:
            output = subprocess.check_output(['ip', 'link', 'show', interface],
                                             stderr=subprocess.STDOUT).decode()
        except subprocess.CalledProcessError:
            print(f"Interface {interface} not found. Please check your CAN hardware connection.")
            return False
        # On Linux, look for 'state UP'
        if 'state UP' not in output:
            try:
                print(f"{interface} is down. Attempting to bring it up at {bitrate} bit/s...")
                subprocess.check_call([
                    'sudo', 'ip', 'link', 'set', interface,
                    'up', 'type', 'can', 'bitrate', str(bitrate)
                ])
            except subprocess.CalledProcessError:
                print(f"Failed to bring up interface {interface}.")
                return False
            print(f"{interface} is now up.")
        else:
            print(f"{interface} is already up.")
        return True

    elif system == 'Darwin':
        # On macOS, use ifconfig to check for the interface.
        try:
            output = subprocess.check_output(['ifconfig', interface],
                                             stderr=subprocess.STDOUT).decode()
        except subprocess.CalledProcessError:
            print(f"Interface {interface} not found on macOS. Please check your CAN hardware connection.")
            return False
        # Check if the interface appears active.
        if "status: active" in output or "UP" in output:
            print(f"{interface} appears active on macOS.")
            return True
        else:
            print(f"{interface} is not active on macOS. Bringing up the interface may not be supported.")
            return False

    else:
        print(f"Unsupported OS: {system}")
        return False

def main():
    # Expecting exactly three arguments: input_data, plant_data, and time_step.
    if len(sys.argv) != 4:
        print("Usage: can_module.py <input_data> <plant_data> <time_step>")
        sys.exit(1)
    try:
        input_data = int(sys.argv[1])
        plant_data = int(sys.argv[2])
        time_step  = int(sys.argv[3])
    except Exception as e:
        print("Error reading arguments:", e)
        print(float('nan'))
        sys.exit(1)

    # Check if the CAN interface is available.
    if not ensure_can_interface():
        print("CAN interface not available.")
        print(float('nan'))
        sys.exit(0)
    
    try:
        # Initialize the CAN bus (using Linux socketcan or appropriate driver)
        bus = can.interface.Bus(channel='can0', bustype='socketcan')
    except Exception as e:
        print("Error initializing CAN bus:", e)
        print(float('nan'))
        sys.exit(0)
    
    try:
        # Prepare the data payload (each value is limited to one byte)
        data = [input_data & 0xFF, plant_data & 0xFF, time_step & 0xFF]
    except Exception as e:
        print("Error converting data values:", e)
        data = [0, 0, 0]
    
    send_id = 0x123
    response_id = 0x124
    msg = can.Message(arbitration_id=send_id, data=data, extended_id=False)
    
    try:
        bus.send(msg)
        print("Sent CAN message with ID {} and data: {}".format(hex(send_id), data))
    except can.CanError as e:
        print("CAN message NOT sent:", e)
        print(float('nan'))
        sys.exit(0)
    
    # Block and wait for a response (up to 1 second)
    start_time = time.time()
    timeout = 1.0
    actuating_signal = None
    while True:
        response = bus.recv(timeout=0.1)
        if response is not None and response.arbitration_id == response_id:
            actuating_signal = response.data[0]
            print("Received response with ID {}: {}".format(hex(response_id), response.data))
            break
        if time.time() - start_time > timeout:
            print("Timeout waiting for CAN response.")
            actuating_signal = float('nan')
            break
    
    try:
        actuating_signal = float(actuating_signal)
    except Exception:
        actuating_signal = float('nan')
    
    # Print the final actuating signal as the last line (MATLAB will capture this)
    print(actuating_signal)
    
    # Clean up the bus (optional)
    try:
        bus.shutdown()
    except Exception:
        pass

if __name__ == '__main__':
    main()