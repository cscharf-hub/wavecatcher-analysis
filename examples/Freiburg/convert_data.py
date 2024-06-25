import sys, argparse
sys.path.append('./examples/Freiburg/gandalf_usb/')
sys.path.append('./examples/Freiburg/gandalf_usb/amc_tools/')
sys.path.append('./examples/Freiburg/gandalf_usb/gandalf_env')
#import reader, decode
import amc_hax

import gzip
import numpy as np

# convert Freiburg GANDALF DAQ binary data to simple binary list
def convert_data(input_file, output_file, end=None):
    print(input_file)
    print(output_file)

    measurement_data = []
    try:
        for k, event in enumerate(amc_hax.frame_events(gzip.open(input_file, 'rb'), pattern=amc_hax.PATTERN_SLINK)):
            for ch, samples in event.frames.items():
                measurement_data.append([k, ch] + samples)
            if end is not None and k == end:
                break
    except Exception as e:
        print('Exception in amc_hax:')
        print(k, e)

    measurement_data = np.array(measurement_data, dtype=np.int32)
    rows, cols = measurement_data.shape

    with open(output_file, "wb") as file:
        file.write(np.array([rows, cols], dtype=np.int32).tobytes())
        measurement_data.tofile(file)