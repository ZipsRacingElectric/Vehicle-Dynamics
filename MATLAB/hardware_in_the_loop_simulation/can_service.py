import sys
import platform
import os
import can
import subprocess
import time
import math
import struct

baudrate = 1000000                      # Baud Rate in bits/s
channel_mac = '/dev/cu.usbmodem14201'   # Path to serial device channel if using macOS
send_id = 0x123                         # Transmit message ID
response_id = 0x124                     # Recieve message ID
timeout = 0.01                          # Maximum time to wait for response
actuating_signal = None

#--------------------------------- FUNCTIONS ------------------------------------

# Check if bus device is there
def check_can_interface(system):
    if system == 'Darwin':
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

# Packs data into int16 types to send in can message
def pack_int16(x):
    # Check that x is within the int16 range.
    if x < -32768 or x > 32767:
        raise ValueError(f"Value {x} out of range for int16")
    # Convert to 2 bytes using big-endian. Use 'little' if required.
    return list(x.to_bytes(2, byteorder='big', signed=True))

#--------------------------------- PROGRAM START ------------------------------------

# Check the operating system
system = platform.system()

# CAN interface configuration
if system == 'Darwin':
    can.rc['interface'] = 'slcan'
    can.rc['channel'] = channel_mac
    can.rc['bitrate'] = baudrate
else:
    sys.exit(0)

# Initialize the CAN bus with the appropriate backend.
try:
    bus = can.interface.Bus()
except Exception as e:
    print("Error initializing CAN bus:", e)
    sys.exit(0)

print("READY", flush=True)  # Let MATLAB know it's ready

while True:
    try:
        line = sys.stdin.readline()
        if not line:
            break  # EOF

        line = line.strip()

        if line == "exit":
            #print("Closing CAN service.", flush=True)
            bus.shutdown()
            break

        parts = line.split()
        if len(parts) != 3:
            print("Error - invalid input", flush=True)
            continue

        # Parse and pack input values
        input_val = pack_int16(int(round(float(parts[0]))))
        plant_val = pack_int16(int(round(float(parts[1]))))
        timestep_val = pack_int16(int(round(float(parts[2]))))
        data = input_val + plant_val + timestep_val

        # Send the CAN message
        msg = can.Message(arbitration_id=send_id, data=data, is_extended_id=False)
        bus.send(msg)
        # print(f"Data sent: {hex(send_id)} {list(data)}", flush=True)

        # Wait for a response
        response = bus.recv(timeout=0.05)
        if response is not None and response.arbitration_id == response_id:
            # Unpack little-endian signed 16-bit int
            actuating_signal = struct.unpack_from('<h', response.data, 0)[0]
            print(actuating_signal, flush=True)
        else:
            print("NaN", flush=True)

    except Exception as e:
        print("NaN", e, flush=True)
        continue