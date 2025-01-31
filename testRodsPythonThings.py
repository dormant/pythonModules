#!/usr/bin/env python
# testRodsPythonThings.py
#
# R.C. Stewart, 18 Aug 2023
#



############  Imports
import os
import sys
import glob
import subprocess
import argparse
import re
import obspy
import warnings
import matplotlib.pyplot as plt
import numpy as np
import socket
from datetime import datetime, date, timedelta, time
from dateutil import parser as dparser
from dateutil.rrule import rrule, DAILY
from obspy.clients.earthworm import Client
from obspy.core import UTCDateTime, Stream, Trace
from pathlib import Path
from obspy.signal.tf_misfit import plot_tfr
from obspy.clients.filesystem.sds import Client as sdsClient


############  Location of python modules is different on opsproc2
hostname = socket.gethostname()
if hostname == 'opsproc2':
    sys.path.append( '/'.join( [os.path.expanduser('~'), 'src/pythonModules' ] ) )
elif hostname == 'opsproc3':
    sys.path.append( '/'.join( [os.path.expanduser('~'), 'src/pythonModules' ] ) )
else:
    sys.path.append( '/'.join( [os.path.expanduser('~'), 'STUFF/src/pythonModules' ] ) )

import rodsPythonThings as rpt
import rodsPlotTfr as rodstfr

dirEvents = "/mnt/mvofls2/Seismic_Data/monitoring_data/events"


for file in os.listdir(dirEvents):
     filename = os.fsdecode(file)
     if filename.startswith("2023"):
        fileEvent = os.path.join(dirEvents, filename)
        st = obspy.read( fileEvent )
        st1 = st.select(station="MSS1")
        tr = st1[0]
        clipped = rpt.isClipped( tr )
        if clipped:
            print( filename, 'clipped' )
     else:
        continue

