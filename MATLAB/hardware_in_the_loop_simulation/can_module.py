import can
import subprocess
import time
import math

try:
    input_data = int(input_data)
    plant_data = int(plant_data)
    time_step = int(time_step)
except Exception as e:
    print("Error reading inputs:", e)
    actuating_signal = float('nan')
else:
    actuating_signal = 42.0  # just a test value