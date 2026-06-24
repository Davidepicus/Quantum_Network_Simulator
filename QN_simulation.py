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
import sequence.topology.node as nd
import sequence.topology.topology as tp
import sequence.topology.qkd_topo as qkdtp
from sequence.kernel.timeline import Timeline


class Counter:
    def __init__(self):
        self.count = 0

    def trigger(self, detector, info):
        self.count += 1



class BBM92_receiver(nd.Node):
  #class that represents a BBM92 receiver
  
  def __init__(self,name="unnamed",det_eff=0.3,dark_counts=5000,count_rate_max=3*10**(6), tl = Timeline(2e10)):
    ###### detector definition ##############
    super().__init__(name, tl)
    self.name=name
    dect0_name= name+ "det0"
    dect1_name= name+ "det1"
    dectA_name= name+ "det+"
    dectB_name= name+ "det-"
    det_1= det.Detector(name= dect0_name, timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max)
    det_2= det.Detector(name= dect1_name, timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max)
    det_3= det.Detector(name= dectA_name, timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max)
    det_4= det.Detector(name= dectB_name, timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max)
    self.add_component(det_1)
    self.add_component(det_2)
    self.add_component(det_3)
    self.add_component(det_4)
    ######### Polarization beam splitter definition ############
    PBS_name= name + "Polarizing_Beam_splitter"
    PBS=BS.BeamSplitter(name=PBS_name, timeline=tl, fidelity=0.98)
    self.add_component(PBS)
    ######### Normal beam splitter definition ################
    BSA_name= name + "Beam_splitter_comp_basis"
    BSB_name= name + "Beam_splitter_x_basis"
    self.BSA_count=0
    self.BSB_count=0
    BS_A =BS.FockBeamSplitter2(name= BSA_name,owner=self.name,efficiency=0.99,timeline=tl,photon_counter=self.BSA_count, src_list=[PBS])
    BS_B =BS.FockBeamSplitter2(name= BSB_name,owner=self.name,efficiency=0.99,timeline=tl,photon_counter=self.BSB_count, src_list=[PBS]) 
    self.add_component(BS_A)
    self.add_component(BS_B)
    ######## Definition of fiber optic cables connecting the several components ###########
    #name_qcin= name + "qcin"
    #name_qcout1= name + "qcout1"
    #name_qcout2= name + "qcout2"
    #name_qcoutA= name + "qcoutA"
    #name_qcoutB= name + "qcoutB"
    #name_qcoutC= name + "qcoutC"
    #name_qcoutD= name + "qcoutD"
    #qcin = optch.QuantumChannel(name=name_qcin, timeline=tl, distance=5, polarization_fidelity=0.99, attenuation=0.0002)
    #qcout1 = optch.QuantumChannel(name=name_qcout1, timeline=tl, distance=5, polarization_fidelity=0.99, attenuation=0.0002)
    #qcout2 = optch.QuantumChannel(name=name_qcout2, timeline=tl, distance=5, polarization_fidelity=0.99, attenuation=0.0002)
    #qcoutA = optch.QuantumChannel(name=name_qcoutA, timeline=tl, distance=5, polarization_fidelity=0.99, attenuation=0.0002)
    #qcoutB = optch.QuantumChannel(name=name_qcoutB, timeline=tl, distance=5, polarization_fidelity=0.99, attenuation=0.0002)
    #qcoutC = optch.QuantumChannel(name=name_qcoutC, timeline=tl, distance=5, polarization_fidelity=0.99, attenuation=0.0002)
    #qcoutD = optch.QuantumChannel(name=name_qcoutD, timeline=tl, distance=5, polarization_fidelity=0.98, attenuation=0.0002)
    #self.add_component(qcin)
    #self.add_component(qcout1)
    #self.add_component(qcout2)
    #self.add_component(qcoutA)
    #self.add_component(qcoutB)
    #self.add_component(qcoutC)
    #self.add_component(qcoutD)
    print("BBM92_receiver named: ",  self.name ," has been set")
    ######## Connecting the hardware ##################
    self.set_first_component(PBS)
    BS_A.add_outputs(outputs=[det_1,det_2])
    BS_B.add_outputs(outputs=[det_3,det_4])

#  def modulus(self):
#    return sqrt(pow(self.re,2)+pow(self.im,2))


class BBM92_SPDC_source(nd.Node):
  #class that represents a BBM92 source
  
  def __init__(self,name="unnamed", freq=10000000.0, tl = Timeline(2e10)):
    super().__init__(name, tl)
    self.name=name
    SPDC_name= name +" SPDC_source"
    SPDC_source = ls.SPDCSource(name=SPDC_name, timeline=tl, wavelengths=[1550,1550], frequency= freq, mean_photon_num=0.1, encoding_type={'bases': None, 'name': 'fock'}, phase_error=0.1)
    self.add_component(SPDC_source)
    print("BBM92_SPDC_source named: ",  self.name ," has been set")


##################### main program definition ################

def main():

  runtime = 2e10
  distance = 1e3
  tl = Timeline(runtime)


#### Defining nodes, Alice and Bob (the receivers) and Charlie (the sender)
  Alice = BBM92_receiver(name="Alice")
  Bob = BBM92_receiver(name="Bob")
  Charlie = BBM92_SPDC_source(name="Charlie")
  Alice.set_seed(0)
  Bob.set_seed(1)
  Charlie.set_seed(2)


#### Defining channels and connecting the nodes
  qc0 = optch.QuantumChannel("qc0", tl, distance=distance, polarization_fidelity=0.97, attenuation=0.0002)
  qc1 = optch.QuantumChannel("qc1", tl, distance=distance, polarization_fidelity=0.97, attenuation=0.0002)
  qc0.set_ends(Charlie, Alice.name)
  qc1.set_ends(Charlie, Bob.name)

  cc0 = optch.ClassicalChannel("cc0", tl, distance=distance)
  cc1 = optch.ClassicalChannel("cc1", tl, distance=distance)
  cc0.set_ends(Alice, Bob.name)
  cc1.set_ends(Bob, Alice.name)


  #qkdtp.QKDTopo(nodes=[Alice, Bob, Charlie], qchannels=[qc0, qc1], cchannels=[cc0,cc1], tl=tl)


if __name__ == "__main__":
  main()



