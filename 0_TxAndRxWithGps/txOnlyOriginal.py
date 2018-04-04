#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: Txonlyoriginal
# Generated: Sat Mar 31 09:58:14 2018
##################################################

from gnuradio import analog
from gnuradio import eng_notation
from gnuradio import gr
from gnuradio import uhd
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from optparse import OptionParser
import time


class txOnlyOriginal(gr.top_block):

    def __init__(self):
        gr.top_block.__init__(self, "Txonlyoriginal")

        ##################################################
        # Variables
        ##################################################
        self.samp_rate = samp_rate = 2e6
        self.TX_freq = TX_freq = 400e6
        self.TX_Gain = TX_Gain = 65

        ##################################################
        # Blocks
        ##################################################
        self.uhd_usrp_sink_0 = uhd.usrp_sink(
        	",".join(("", "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(1),
        	),
        )
        self.uhd_usrp_sink_0.set_samp_rate(samp_rate)
        self.uhd_usrp_sink_0.set_center_freq(TX_freq, 0)
        self.uhd_usrp_sink_0.set_gain(TX_Gain, 0)
        self.analog_const_source_x_0 = analog.sig_source_c(0, analog.GR_CONST_WAVE, 0, 0, 1)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_const_source_x_0, 0), (self.uhd_usrp_sink_0, 0))

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.uhd_usrp_sink_0.set_samp_rate(self.samp_rate)

    def get_TX_freq(self):
        return self.TX_freq

    def set_TX_freq(self, TX_freq):
        self.TX_freq = TX_freq
        self.uhd_usrp_sink_0.set_center_freq(self.TX_freq, 0)

    def get_TX_Gain(self):
        return self.TX_Gain

    def set_TX_Gain(self, TX_Gain):
        self.TX_Gain = TX_Gain
        self.uhd_usrp_sink_0.set_gain(self.TX_Gain, 0)



def main(top_block_cls=txOnlyOriginal, options=None):

    tb = top_block_cls()
    tb.start()
    try:
        raw_input('Press Enter to quit: ')
    except EOFError:
        pass
    tb.stop()
    tb.wait()


if __name__ == '__main__':
    main()
