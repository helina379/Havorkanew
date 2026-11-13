"""Logs the Hovorka board's serial output to a CSV, skipping garbled lines.

Setup (once):   pip install pyserial
Run:            python run_log.py
Stop:           Ctrl+C

Close the Arduino Serial Monitor first -- only one program can use the port.
Opening the port resets the Hovorka board, so the simulation restarts at t=0.
"""
import datetime
import serial

PORT = "COM12"      # change if the Hovorka board shows up on a different port
BAUD = 115200

fname = datetime.datetime.now().strftime("hovorka_run_%Y%m%d_%H%M.csv")


def valid(line):
    parts = line.split(",")
    if len(parts) != 3:
        return False
    try:
        [float(p) for p in parts]
        return True
    except ValueError:
        return False


with serial.Serial(PORT, BAUD, timeout=2) as ser, open(fname, "w", newline="") as f:
    f.write("time_min,glucose_mgdl,insulin_mUmin\n")
    print("Logging to", fname, "- press Ctrl+C to stop")
    try:
        while True:
            line = ser.readline().decode("ascii", errors="ignore").strip()
            if valid(line):
                f.write(line + "\n")
                f.flush()
                print(line)
    except KeyboardInterrupt:
        print("Stopped. Saved", fname)
