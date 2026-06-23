# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 10:09:08 2026

@author: giorg
"""

 
import matplotlib.pyplot as plt
import numpy as np
import sequence.components.optical_channel as optch
import sequence.components.detector as det
import sequence.components.beam_splitter as BS
import sequence.components.light_source as ls

from sequence.kernel.event import Event
from sequence.kernel.process import Process
from sequence.kernel.timeline import Timeline
from sequence.qkd.BB84 import pair_bb84_protocols
from sequence.topology.node import QKDNode
import sequence.utils.log as log




class BBM92_receiver:
  #class that represents a BBM92 receiver
  
  def __init__(self,det_eff=0.3,dark_counts=5000,count_rate_max=3*10**(6), tl = Timeline(2e10)):
    ###### detector definition ##############
    self.det_1= det.Detector(timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max)
    self.det_2= det.Detector(timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max)
    self.det_3= det.Detector(timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max)
    self.det_4= det.Detector(timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max)
    ######### Polarization beam splitter definition ############
    self.PBS = BS.BeamSplitter(timeline=tl, fidelity=0.98)
    ######### Normal beam splitter definition ################
    self.BS_A =BS.FockBeamSplitter2(timeline=tl, fidelity=0.98)
    self.BS_B =BS.FockBeamSplitter2(timeline=tl, fidelity=0.98)
    ######## Fiber optic connecting the several components definition ###########
    self.qcin = optch.QuantumChannel("qcin", tl, distance=5, polarization_fidelity=0.99, attenuation=0.0002)
    self.qcout1 = optch.QuantumChannel("qcout1", tl, distance=5, polarization_fidelity=0.99, attenuation=0.0002)
    self.qcout2 = optch.QuantumChannel("qcout2", tl, distance=5, polarization_fidelity=0.99, attenuation=0.0002)
    self.qcoutA = optch.QuantumChannel("qcoutA", tl, distance=5, polarization_fidelity=0.99, attenuation=0.0002)
    self.qcoutB = optch.QuantumChannel("qcoutB", tl, distance=5, polarization_fidelity=0.99, attenuation=0.0002)
    self.qcoutC = optch.QuantumChannel("qcoutC", tl, distance=5, polarization_fidelity=0.99, attenuation=0.0002)
    self.qcoutD = optch.QuantumChannel("qcoutD", tl, distance=5, polarization_fidelity=0.98, attenuation=0.0002)



#  def modulus(self):
#    return sqrt(pow(self.re,2)+pow(self.im,2))


class BBM92_SPDC_source:
  #class that represents a BBM92 source
  
  def __init__(self,freq=10000000.0, tl = Timeline(2e10)):
      self.SPDC = ls.SPDCSource(timeline=tl, wavelengths=[1550,1550], frequency= freq, mean_photon_num=0.1, encoding_type={'bases': None, 'name': 'fock'}, phase_error=0.1, linewidth=5)




runtime = 2e10
distance = 1e3

tl = Timeline(runtime)

qc0 = optch.QuantumChannel("qc0", tl, distance=distance, polarization_fidelity=0.97, attenuation=0.0002)
qc1 = optch.QuantumChannel("qc1", tl, distance=distance, polarization_fidelity=0.97, attenuation=0.0002)
cc0 = optch.ClassicalChannel("cc0", tl, distance=distance)
cc1 = optch.ClassicalChannel("cc1", tl, distance=distance)




