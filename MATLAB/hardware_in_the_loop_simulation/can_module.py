import sys
import platform
import os
import can
import subprocess
import time
import math

baudrate = 1000000                      # Baud Rate in bits/s
channel_mac = '/dev/cu.usbmodem14101'   # Path to serial device channel if using macOS
send_id = 0x123                         # Transmit message ID
response_id = 0x124                     # Recieve message ID
timeout = 0.01                           # Maximum time to wait for response

def check_can_interface(system="none"):
    """
    Checks if the CAN interface is available.
    For Linux (socketcan), uses the 'ip' command to check interface status.
    For macOS (Darwin), checks if the serial device exists.
    Returns True if the interface/device is available, else False.
    """
    if system == 'Linux':
        try:
            output = subprocess.check_output(['ip', 'link', 'show', can.rc['channel']],
                                             stderr=subprocess.STDOUT).decode()
        except subprocess.CalledProcessError:
            print(f"Channel {can.rc['channel']} not found. Please check your CAN hardware connection.")
            return False
        # Check if interface is UP
        if 'state UP' not in output:
            try:
                print(f"{can.rc['channel']} is down. Attempting to bring it up at {can.rc['bitrate']} bit/s...")
                subprocess.check_call([
                    'sudo', 'ip', 'link', 'set', can.rc['channel'],
                    'up', 'type', 'can', 'bitrate', str(can.rc['bitrate'])
                ])
            except subprocess.CalledProcessError:
                print(f"Failed to bring up interface {can.rc['channel']}.")
                return False
            print(f"{can.rc['channel']} is now up.")
        else:
            print(f"{can.rc['channel']} is already up.")
        return True

    elif system == 'Darwin':
        if not os.path.exists(can.rc['channel']):
            print(f"Serial channel {can.rc['channel']} not found. Please check your CAN hardware connection.")
            return False
        else:
            print(f"Serial channel {can.rc['channel']} found on macOS.")
            # On Darwin no need to bring the interface up beforehand.
            return True
    else:
        print(f"Unsupported OS: {system}")
        return False
    
def pack_int16(x):
    # Check that x is within the int16 range.
    if x < -32768 or x > 32767:
        raise ValueError(f"Value {x} out of range for int16")
    # Convert to 2 bytes using big-endian. Use 'little' if required.
    return list(x.to_bytes(2, byteorder='big', signed=True))

def main():
    # Expecting exactly three arguments: input_data, plant_data, and time_step.
    if len(sys.argv) != 4:
        print("Usage: can_module.py <input_data> <plant_data> <time_step>")
        sys.exit(1)
    try:
        # Data will be converted to the equivalent to int16_t in C:
        input_data = pack_int16(int(round(float(sys.argv[1]))))
        plant_data = pack_int16(int(round(float(sys.argv[2]))))
        time_step  = pack_int16(int(round(float(sys.argv[3]))))
    except Exception as e:
        print("Error reading arguments:", e)
        print(float('nan'))
        sys.exit(1)
    
    # Check the operating system
    system = platform.system()

    if system == 'Linux':
        print("Linux system detected.")
        can.rc['interface'] = 'socketcan'
        can.rc['channel'] = 'vcan0'
        can.rc['bitrate'] = baudrate
    elif system == 'Darwin':
        print("macOS system detected.")
        can.rc['interface'] = 'slcan'
        can.rc['channel'] = channel_mac
        can.rc['bitrate'] = baudrate
    else:
        print("Unsupported OS.")
        print(float('nan'))
        sys.exit(0)
    
    # Check if the CAN interface is available.
    if not check_can_interface(system):
        print("CAN interface not available.")
        print(float('nan'))
        sys.exit(0)
    
    # Initialize the CAN bus with the appropriate backend.
    try:
        bus = can.interface.Bus()
    except Exception as e:
        print("Error initializing CAN bus:", e)
        print(float('nan'))
        sys.exit(0)
    
    # Create CAN message
    try:
        # Prepare the data payload (each value constrained to two byte)
        data = input_data + plant_data + time_step
    except Exception as e:
        print("Error converting data values:", e)
        data = [0, 0, 0]
    
    msg = can.Message(arbitration_id=send_id, data=data, is_extended_id=False)
    
    # Send CAN message
    try:
        bus.send(msg)
        print("Sent CAN message with ID {} and data: {}".format(hex(send_id), data))
    except can.CanError as e:
        print("CAN message NOT sent:", e)
        print(float('nan'))
        sys.exit(0)
    
    # Block and wait for a response (timeout: 1 second)
    start_time = time.time()
    actuating_signal = None
    while True:
        response = bus.recv(timeout)
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
    
    # Print the final actuating signal (MATLAB can capture the last printed line)
    print(actuating_signal)
    
    # Clean up the bus (optional)
    try:
        bus.shutdown()
    except Exception:
        pass

if __name__ == '__main__':
    main()
