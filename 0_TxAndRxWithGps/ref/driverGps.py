from sequential_pixhawk_process import DroneMetadataGenerator
from time import sleep
from config.config import get_config
from config.increment_flight_number import  increment_flight_number
from sigmf_utils import epoch_time, update_and_save_sigmf_global
from grc.ReceiveBlock import ReceiveBlock
from grc.TransmitBlock import TransmitBlock
import argparse
import json


# save GPS time and CPU time to factor in offset when post processing
def create_gps_offset_file(config, usrp):
    usrp = rx_block.get_usrp_block()

    gps_offset = {
        'gps_locked': False
    }

    with open(config['sdr']['gps_offset_file'], 'a') as gps_offset_file:
        try:
            gps_locked = usrp.get_mboard_sensor('gps_locked').to_bool()

            if gps_locked:
                gps_offset['gps_locked'] = True
                gps_time = usrp.get_mboard_sensor('gps_time').to_real()
                cpu_time = epoch_time()

                gps_offset['gps_time'] = gps_time
                gps_offset['cpu_time'] = cpu_time

        except Exception as exception:
            gps_offset['exception_message'] = str(exception)

        json.dump(gps_offset, gps_offset_file, indent=1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Driver to run Pixhawk process and GNU radio flow graph')

    parser.add_argument('delay_seconds', type=float,
                        help='Number of seconds to delay starting the flow graph and pixhawk metadata generator process')

    parser.add_argument('run_seconds', type=float,
                        help='Number of seconds to run the flow graph and pixhawk metadata generator process')

    parser.add_argument('description', help='description that will be set in SigMf global object')
    args = parser.parse_args()

    config = get_config('./config/sdr_config.json')

    rx_block = ReceiveBlock(config)
    tx_block = TransmitBlock(config)

    create_gps_offset_file(config, rx_block)

    # pixhawk_process = DroneMetadataGenerator(args.run_seconds, config)

    # update and save sig-mf global JSON file
    update_and_save_sigmf_global(config['sigmf'], rx_block.get_samp_rate(),
                                 args.description, config['sdr']['global_file'])

    tx_block.start()
    sleep(args.delay_seconds)

    rx_block.start()
    # pixhawk_process.start()
    # pixhawk_process.join() # will run for args.run_seconds

    sleep(args.run_seconds)  # used when not testing with pixhawk process, thus join is not being called

    tx_block.stop()
    rx_block.stop()

    increment_flight_number()
