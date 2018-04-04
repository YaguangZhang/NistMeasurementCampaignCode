# Get the mounted GPS on the USRP and log the GPS info into a file (Linux new
# line is used: \n).
#
# Yaguang Zhang, Purdue University, 2017-06-15

import time, datetime
from gnuradio import uhd

DEBUG_DIV_MARKER = '>>>>>>>>>>>>>>>>>>'
class gpsdo:
    def __init__(self, usrp, LOG_FILE_PATH):
      mboardSersorNames = usrp.get_mboard_sensor_names()
      self.usrp = usrp
      self.logFilePath = LOG_FILE_PATH

      print(DEBUG_DIV_MARKER)
      print('Initialized GPSDO.')
      print('')
      print('    Available sersors: ')
      print(mboardSersorNames)
      print('')
      print('    Log file path:')
      print(self.logFilePath)
      print('')
      return

    # Log the GPS info together with channel gain.
    def log(self, rxChannelGain):
        print('<<<<<<<<<<<<<<<<<<')
        print('GPSDO: logging...')
        print('')

        # Get all the information needed to be logged.
        sysEpochTime = time.time()
        sysDateAndTime = datetime.datetime.fromtimestamp(sysEpochTime).strftime('%Y-%m-%d %H:%M:%S')
        refLocked = 1 if self.usrp.get_mboard_sensor('ref_locked').value.strip()=='true' else 0
        gpsLocked = 1 if self.usrp.get_mboard_sensor('gps_locked').value.strip()=='true' else 0
        gpsTime = self.usrp.get_mboard_sensor('gps_time').value.strip()
        gpsLocation = self.usrp.get_mboard_sensor('gps_gpgga').value.strip()

        # Construct the strings to write into the file.
        sampleLine = ', '.join([sysDateAndTime, str(sysEpochTime), str(rxChannelGain), \
            str(refLocked), str(gpsLocked), gpsTime, \
            '{'+gpsLocation+'}'])
        linesToWrite = ['sysDateAndTime, sysEpochTime, rxChannelGain, refLocked, gpsLocked, gpsTime, {gpsLocation}', \
            sampleLine]

        # Write lines to file.
        print('    Fetched data: ')
        print(linesToWrite[0])
        print(linesToWrite[1])
        print('')

        hLogFile = open(self.logFilePath, 'w')
        hLogFile.writelines( (line+'\n') for line in linesToWrite)
        hLogFile.close()

        print('Done!')
        print(DEBUG_DIV_MARKER)