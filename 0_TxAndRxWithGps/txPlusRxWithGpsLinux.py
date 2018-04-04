#!/usr/bin/env python2
# -*- coding: utf-8 -*-
##################################################
# GNU Radio Python Flow Graph
# Title: txPlusRxWithGpsOriginal
# Generated: Fri Mar 30 20:46:54 2018
##################################################

if __name__ == '__main__':
    import ctypes
    import sys
    if sys.platform.startswith('linux'):
        try:
            x11 = ctypes.cdll.LoadLibrary('libX11.so')
            x11.XInitThreads()
        except:
            print "Warning: failed to XInitThreads()"

from PyQt4 import Qt
from gnuradio import analog
from gnuradio import blocks
from gnuradio import eng_notation
from gnuradio import filter
from gnuradio import gr
from gnuradio import qtgui
from gnuradio import uhd
from gnuradio.eng_option import eng_option
from gnuradio.filter import firdes
from gnuradio.qtgui import Range, RangeWidget
from optparse import OptionParser
import sip
import sys
import time

# ZYG
import os
from threading import Timer
from lib import gpsdo
curFileDir = os.path.dirname(os.path.realpath(__file__))
FLAG_TRY_TO_LOCK_GPSDO = False

class txPlusRxWithGpsOriginal(gr.top_block, Qt.QWidget):

    def __init__(self):
        gr.top_block.__init__(self, "txPlusRxWithGpsOriginal")
        Qt.QWidget.__init__(self)
        self.setWindowTitle("txPlusRxWithGpsOriginal")
        try:
            self.setWindowIcon(Qt.QIcon.fromTheme('gnuradio-grc'))
        except:
            pass
        self.top_scroll_layout = Qt.QVBoxLayout()
        self.setLayout(self.top_scroll_layout)
        self.top_scroll = Qt.QScrollArea()
        self.top_scroll.setFrameStyle(Qt.QFrame.NoFrame)
        self.top_scroll_layout.addWidget(self.top_scroll)
        self.top_scroll.setWidgetResizable(True)
        self.top_widget = Qt.QWidget()
        self.top_scroll.setWidget(self.top_widget)
        self.top_layout = Qt.QVBoxLayout(self.top_widget)
        self.top_grid_layout = Qt.QGridLayout()
        self.top_layout.addLayout(self.top_grid_layout)

        self.settings = Qt.QSettings("GNU Radio", "txPlusRxWithGpsOriginal")
        self.restoreGeometry(self.settings.value("geometry").toByteArray())

        ##################################################
        # Variables
        ##################################################
        self.time_sink_trigger_level = time_sink_trigger_level = 0.002
        self.samp_rate = samp_rate = 2.0e6
        self.TX_freq = TX_freq = 400e6
        self.TX_Gain = TX_Gain = 65
        self.RX_freq = RX_freq = 2.5e9
        self.RX_Gain = RX_Gain = 50

        ##################################################
        # Blocks
        ##################################################
        self._time_sink_trigger_level_range = Range(0.0001, 1, 0.0001, 0.002, 200)
        self._time_sink_trigger_level_win = RangeWidget(self._time_sink_trigger_level_range, self.set_time_sink_trigger_level, "time_sink_trigger_level", "counter_slider", float)
        self.top_layout.addWidget(self._time_sink_trigger_level_win)
        self._TX_Gain_range = Range(0, 100, 1, 65, 200)
        self._TX_Gain_win = RangeWidget(self._TX_Gain_range, self.set_TX_Gain, "TX_Gain", "counter_slider", float)
        self.top_layout.addWidget(self._TX_Gain_win)
        self._RX_Gain_range = Range(0, 100, 1, 50, 200)
        self._RX_Gain_win = RangeWidget(self._RX_Gain_range, self.set_RX_Gain, "RX_Gain", "counter_slider", float)
        self.top_layout.addWidget(self._RX_Gain_win)
        self.uhd_usrp_source_0_0 = uhd.usrp_source(
        	",".join(("", "")),
        	uhd.stream_args(
        		cpu_format="fc32",
        		channels=range(1),
        	),
        )
        self.uhd_usrp_source_0_0.set_samp_rate(samp_rate)
        self.uhd_usrp_source_0_0.set_center_freq(RX_freq, 0)
        self.uhd_usrp_source_0_0.set_gain(RX_Gain, 0)
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
        self.qtgui_time_sink_x_0_0 = qtgui.time_sink_c(
        	1024, #size
        	samp_rate, #samp_rate
        	"After LPF", #name
        	1 #number of inputs
        )
        self.qtgui_time_sink_x_0_0.set_update_time(0.10)
        self.qtgui_time_sink_x_0_0.set_y_axis(-1, 1)

        self.qtgui_time_sink_x_0_0.set_y_label("Amplitude", "")

        self.qtgui_time_sink_x_0_0.enable_tags(-1, True)
        self.qtgui_time_sink_x_0_0.set_trigger_mode(qtgui.TRIG_MODE_NORM, qtgui.TRIG_SLOPE_POS, time_sink_trigger_level, 0, 0, "")
        self.qtgui_time_sink_x_0_0.enable_autoscale(False)
        self.qtgui_time_sink_x_0_0.enable_grid(False)
        self.qtgui_time_sink_x_0_0.enable_control_panel(True)

        if not True:
          self.qtgui_time_sink_x_0_0.disable_legend()

        labels = ["", "", "", "", "",
                  "", "", "", "", ""]
        widths = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        colors = ["blue", "red", "green", "black", "cyan",
                  "magenta", "yellow", "dark red", "dark green", "blue"]
        styles = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        markers = [-1, -1, -1, -1, -1,
                   -1, -1, -1, -1, -1]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]

        for i in xrange(2*1):
            if len(labels[i]) == 0:
                if(i % 2 == 0):
                    self.qtgui_time_sink_x_0_0.set_line_label(i, "Re{{Data {0}}}".format(i/2))
                else:
                    self.qtgui_time_sink_x_0_0.set_line_label(i, "Im{{Data {0}}}".format(i/2))
            else:
                self.qtgui_time_sink_x_0_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_0_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_0_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_0_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_0_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_0_0_win = sip.wrapinstance(self.qtgui_time_sink_x_0_0.pyqwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_0_0_win)
        self.qtgui_time_sink_x_0 = qtgui.time_sink_c(
        	1024, #size
        	samp_rate, #samp_rate
        	"Signal recevied", #name
        	1 #number of inputs
        )
        self.qtgui_time_sink_x_0.set_update_time(0.10)
        self.qtgui_time_sink_x_0.set_y_axis(-1, 1)

        self.qtgui_time_sink_x_0.set_y_label("Amplitude", "")

        self.qtgui_time_sink_x_0.enable_tags(-1, True)
        self.qtgui_time_sink_x_0.set_trigger_mode(qtgui.TRIG_MODE_NORM, qtgui.TRIG_SLOPE_POS, time_sink_trigger_level, 0, 0, "")
        self.qtgui_time_sink_x_0.enable_autoscale(False)
        self.qtgui_time_sink_x_0.enable_grid(False)
        self.qtgui_time_sink_x_0.enable_control_panel(True)

        if not True:
          self.qtgui_time_sink_x_0.disable_legend()

        labels = ["", "", "", "", "",
                  "", "", "", "", ""]
        widths = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        colors = ["blue", "red", "green", "black", "cyan",
                  "magenta", "yellow", "dark red", "dark green", "blue"]
        styles = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        markers = [-1, -1, -1, -1, -1,
                   -1, -1, -1, -1, -1]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]

        for i in xrange(2*1):
            if len(labels[i]) == 0:
                if(i % 2 == 0):
                    self.qtgui_time_sink_x_0.set_line_label(i, "Re{{Data {0}}}".format(i/2))
                else:
                    self.qtgui_time_sink_x_0.set_line_label(i, "Im{{Data {0}}}".format(i/2))
            else:
                self.qtgui_time_sink_x_0.set_line_label(i, labels[i])
            self.qtgui_time_sink_x_0.set_line_width(i, widths[i])
            self.qtgui_time_sink_x_0.set_line_color(i, colors[i])
            self.qtgui_time_sink_x_0.set_line_style(i, styles[i])
            self.qtgui_time_sink_x_0.set_line_marker(i, markers[i])
            self.qtgui_time_sink_x_0.set_line_alpha(i, alphas[i])

        self._qtgui_time_sink_x_0_win = sip.wrapinstance(self.qtgui_time_sink_x_0.pyqwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_time_sink_x_0_win)
        self.qtgui_freq_sink_x_0_0 = qtgui.freq_sink_c(
        	4096, #size
        	firdes.WIN_BLACKMAN_hARRIS, #wintype
        	RX_freq, #fc
        	samp_rate, #bw
        	"After LPF", #name
        	1 #number of inputs
        )
        self.qtgui_freq_sink_x_0_0.set_update_time(0.10)
        self.qtgui_freq_sink_x_0_0.set_y_axis(-140, 10)
        self.qtgui_freq_sink_x_0_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, 0.0, 0, "")
        self.qtgui_freq_sink_x_0_0.enable_autoscale(False)
        self.qtgui_freq_sink_x_0_0.enable_grid(False)
        self.qtgui_freq_sink_x_0_0.set_fft_average(1.0)
        self.qtgui_freq_sink_x_0_0.enable_control_panel(True)

        if not True:
          self.qtgui_freq_sink_x_0_0.disable_legend()

        if "complex" == "float" or "complex" == "msg_float":
          self.qtgui_freq_sink_x_0_0.set_plot_pos_half(not True)

        labels = ["", "", "", "", "",
                  "", "", "", "", ""]
        widths = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        colors = ["blue", "red", "green", "black", "cyan",
                  "magenta", "yellow", "dark red", "dark green", "dark blue"]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]
        for i in xrange(1):
            if len(labels[i]) == 0:
                self.qtgui_freq_sink_x_0_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_freq_sink_x_0_0.set_line_label(i, labels[i])
            self.qtgui_freq_sink_x_0_0.set_line_width(i, widths[i])
            self.qtgui_freq_sink_x_0_0.set_line_color(i, colors[i])
            self.qtgui_freq_sink_x_0_0.set_line_alpha(i, alphas[i])

        self._qtgui_freq_sink_x_0_0_win = sip.wrapinstance(self.qtgui_freq_sink_x_0_0.pyqwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_freq_sink_x_0_0_win)
        self.qtgui_freq_sink_x_0 = qtgui.freq_sink_c(
        	4096, #size
        	firdes.WIN_BLACKMAN_hARRIS, #wintype
        	RX_freq, #fc
        	samp_rate, #bw
        	"Signal recevied", #name
        	1 #number of inputs
        )
        self.qtgui_freq_sink_x_0.set_update_time(0.10)
        self.qtgui_freq_sink_x_0.set_y_axis(-140, 10)
        self.qtgui_freq_sink_x_0.set_trigger_mode(qtgui.TRIG_MODE_FREE, 0.0, 0, "")
        self.qtgui_freq_sink_x_0.enable_autoscale(False)
        self.qtgui_freq_sink_x_0.enable_grid(False)
        self.qtgui_freq_sink_x_0.set_fft_average(1.0)
        self.qtgui_freq_sink_x_0.enable_control_panel(True)

        if not True:
          self.qtgui_freq_sink_x_0.disable_legend()

        if "complex" == "float" or "complex" == "msg_float":
          self.qtgui_freq_sink_x_0.set_plot_pos_half(not True)

        labels = ["", "", "", "", "",
                  "", "", "", "", ""]
        widths = [1, 1, 1, 1, 1,
                  1, 1, 1, 1, 1]
        colors = ["blue", "red", "green", "black", "cyan",
                  "magenta", "yellow", "dark red", "dark green", "dark blue"]
        alphas = [1.0, 1.0, 1.0, 1.0, 1.0,
                  1.0, 1.0, 1.0, 1.0, 1.0]
        for i in xrange(1):
            if len(labels[i]) == 0:
                self.qtgui_freq_sink_x_0.set_line_label(i, "Data {0}".format(i))
            else:
                self.qtgui_freq_sink_x_0.set_line_label(i, labels[i])
            self.qtgui_freq_sink_x_0.set_line_width(i, widths[i])
            self.qtgui_freq_sink_x_0.set_line_color(i, colors[i])
            self.qtgui_freq_sink_x_0.set_line_alpha(i, alphas[i])

        self._qtgui_freq_sink_x_0_win = sip.wrapinstance(self.qtgui_freq_sink_x_0.pyqwidget(), Qt.QWidget)
        self.top_layout.addWidget(self._qtgui_freq_sink_x_0_win)
        self.low_pass_filter_0 = filter.fir_filter_ccf(1, firdes.low_pass(
        	1, samp_rate, 60e3, 5e3, firdes.WIN_HAMMING, 6.76))

        # ZYG
        self.epochTimeStrForLogFile = str(int(time.time()))
        outFileName = 'measureSignalCont_'+self.epochTimeStrForLogFile+'.out'
        outFileFilteredName = 'measureSignalCont_'+self.epochTimeStrForLogFile+'_filtered.out'
        self.outFilesPath = os.path.join(curFileDir, 'measureSignalOutput')
        outFilePath = os.path.join(self.outFilesPath, outFileName)
        outFileFilteredPath = os.path.join(self.outFilesPath, outFileFilteredName)

        # ZYG
        self.blocks_file_sink_0_0 = blocks.file_sink(gr.sizeof_gr_complex*1, outFileFilteredPath, False)
        self.blocks_file_sink_0_0.set_unbuffered(False)
        self.blocks_file_sink_0 = blocks.file_sink(gr.sizeof_gr_complex*1, outFilePath, False)
        self.blocks_file_sink_0.set_unbuffered(False)

        self.analog_const_source_x_0 = analog.sig_source_c(0, analog.GR_CONST_WAVE, 0, 0, 1)

        ##################################################
        # Connections
        ##################################################
        self.connect((self.analog_const_source_x_0, 0), (self.uhd_usrp_sink_0, 0))
        self.connect((self.low_pass_filter_0, 0), (self.blocks_file_sink_0_0, 0))
        self.connect((self.low_pass_filter_0, 0), (self.qtgui_freq_sink_x_0_0, 0))
        self.connect((self.low_pass_filter_0, 0), (self.qtgui_time_sink_x_0_0, 0))
        self.connect((self.uhd_usrp_source_0_0, 0), (self.blocks_file_sink_0, 0))
        self.connect((self.uhd_usrp_source_0_0, 0), (self.low_pass_filter_0, 0))
        self.connect((self.uhd_usrp_source_0_0, 0), (self.qtgui_freq_sink_x_0, 0))
        self.connect((self.uhd_usrp_source_0_0, 0), (self.qtgui_time_sink_x_0, 0))

    def closeEvent(self, event):
        self.settings = Qt.QSettings("GNU Radio", "txPlusRxWithGpsOriginal")
        self.settings.setValue("geometry", self.saveGeometry())
        event.accept()


    def get_time_sink_trigger_level(self):
        return self.time_sink_trigger_level

    def set_time_sink_trigger_level(self, time_sink_trigger_level):
        self.time_sink_trigger_level = time_sink_trigger_level
        self.qtgui_time_sink_x_0.set_trigger_mode(qtgui.TRIG_MODE_NORM, qtgui.TRIG_SLOPE_POS, self.time_sink_trigger_level, 0, 0, "")
        self.qtgui_time_sink_x_0_0.set_trigger_mode(qtgui.TRIG_MODE_NORM, qtgui.TRIG_SLOPE_POS, self.time_sink_trigger_level, 0, 0, "")

    def get_samp_rate(self):
        return self.samp_rate

    def set_samp_rate(self, samp_rate):
        self.samp_rate = samp_rate
        self.low_pass_filter_0.set_taps(firdes.low_pass(1, self.samp_rate, 60e3, 5e3, firdes.WIN_HAMMING, 6.76))
        self.qtgui_freq_sink_x_0.set_frequency_range(self.RX_freq, self.samp_rate)
        self.qtgui_freq_sink_x_0_0.set_frequency_range(self.RX_freq, self.samp_rate)
        self.qtgui_time_sink_x_0.set_samp_rate(self.samp_rate)
        self.qtgui_time_sink_x_0_0.set_samp_rate(self.samp_rate)
        self.uhd_usrp_sink_0.set_samp_rate(self.samp_rate)
        self.uhd_usrp_source_0_0.set_samp_rate(self.samp_rate)

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


    def get_RX_freq(self):
        return self.RX_freq

    def set_RX_freq(self, RX_freq):
        self.RX_freq = RX_freq
        self.qtgui_freq_sink_x_0.set_frequency_range(self.RX_freq, self.samp_rate)
        self.qtgui_freq_sink_x_0_0.set_frequency_range(self.RX_freq, self.samp_rate)
        self.uhd_usrp_source_0_0.set_center_freq(self.RX_freq, 0)

    def get_RX_Gain(self):
        return self.RX_Gain

    def set_RX_Gain(self, RX_Gain):
        self.RX_Gain = RX_Gain
        self.uhd_usrp_source_0_0.set_gain(self.RX_Gain, 0)



def main(top_block_cls=txPlusRxWithGpsOriginal, options=None):

    from distutils.version import StrictVersion
    if StrictVersion(Qt.qVersion()) >= StrictVersion("4.5.0"):
        style = gr.prefs().get_string('qtgui', 'style', 'raster')
        Qt.QApplication.setGraphicsSystem(style)
    qapp = Qt.QApplication(sys.argv)

    tb = top_block_cls()

    # ZYG: Move quitting here and add new one for debugging.
    def quitting():
            tb.stop()
            tb.wait()
    def logGps():
        gpsDo = gpsdo.gpsdo(tb.uhd_usrp_source_0_0, \
            os.path.join(tb.outFilesPath, \
            'measureSignalCont_'+tb.epochTimeStrForLogFile+'_GPS.log'))
        gpsDo.log(tb.get_RX_Gain())
    def logGpsCont():
        gpsDo = gpsdo.gpsdo(tb.uhd_usrp_source_0_0, \
            os.path.join(tb.outFilesPath, \
            'measureSignalCont_'+str(int(time.time()))+'_GPS.log'))
        gpsDo.log(tb.get_RX_Gain())
        t = Timer(1.0, logGpsCont)
        t.start()
    keepProbingGps = True
    def probeGps():
        logGps()
        time.sleep(3)
        if keepProbingGps:
            probeGps()
    def quiteFromTimer():
        print(tb.uhd_usrp_source_0_0)
        logGps()
        print('Quitting by Timer...')
        keepProbingGps = False
        quitting()
        qapp.quit()
    if FLAG_TRY_TO_LOCK_GPSDO:
        print(' ')
        print('Trying to lock GPSDO to satellite...')
        print(' ')
        print('    Available sensors:')
        print(tb.uhd_usrp_source_0.get_mboard_sensor_names())
        print(' ')
        time.sleep(3)
        # t = Timer(240.0, quiteFromTimer)
        probeGps()
    # else:
        # t = Timer(3.0, logGpsNewFile)

    tb.start()
    tb.show()

    # ZYG
    logGpsCont()
    # t.start()

    # ZYG
    # def quitting():
    #     tb.stop()
    #     tb.wait()

    qapp.connect(qapp, Qt.SIGNAL("aboutToQuit()"), quitting)
    qapp.exec_()


if __name__ == '__main__':
    main()
