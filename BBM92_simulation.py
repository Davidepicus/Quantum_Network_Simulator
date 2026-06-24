
import math as mt
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
from sequence.kernel.event import Event
from sequence.kernel.process import Process

################# Class that represents a BBM92 receiver ################################
class BBM92_receiver(nd.Node):
  
  def __init__(self,name="unnamed",det_eff=0.3,dark_counts=5000,count_rate_max=3*10**(6), tl = Timeline(2e10)):
    ###### detector definition ##############
    super().__init__(name, tl)
    self.name=name
    dect0_name= name+ "det0"
    dect1_name= name+ "det1"
    self.det_1= det.Detector(name= dect0_name, timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max)
    self.det_2= det.Detector(name= dect1_name, timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max)
    self.add_component(self.det_1)
    self.add_component(self.det_2)

    ######### Polarization beam splitter definition ############
    PBS_name= name + "Polarizing_Beam_splitter"
    self.PBS=BS.BeamSplitter(name=PBS_name, timeline=tl, fidelity=0.98)
    self.add_component(self.PBS)
    self.PBS.set_basis_list(basis_list=[0,1], start_time=0, frequency=0)

    print("BBM92_receiver named: ",  self.name ," has been set")
    ######## Connecting the hardware ##################
    self.set_first_component(self.PBS)
    self.PBS.add_receiver(self.det_1)
    self.PBS.add_receiver(self.det_2)

  def get_counts_D1(self):
    return self.det_1.photon_counter

  def get_counts_D2(self):
    return self.det_2.photon_counter

############### class that represents a BBM92 source ###################################

class BBM92_SPDC_source(nd.Node):
  
  def __init__(self,name="unnamed", freq=10000000.0, tl = Timeline(2e10)):
  ###### source definition ##############
    super().__init__(name, tl)
    self.name=name
    SPDC_name= name +" SPDC_source"
    self.SPDC_source = ls.SPDCSource(name=SPDC_name, timeline=tl, wavelengths=[1550,1550], frequency= freq, mean_photon_num=0.1, encoding_type={'bases': [((1 + 0j, 0j), (0j, 1 + 0j)), ((0.7071067811865476 + 0j, 0.7071067811865476 + 0j), (-0.7071067811865476 + 0j, 0.7071067811865476 + 0j))], 'name': 'polarization'}, phase_error=0.1)
    self.add_component(self.SPDC_source)
    print("BBM92_SPDC_source named: ",  self.name ," has been set")

  def emit_photon(self):
         c = 1 / mt.sqrt(2)
         self.SPDC_source.emit([[c+0j,c+0j],[c+0j,-c+0j]])


##################### main program definition ############################

def main():

  runtime = 2e10
  distance = 1e3
  tl = Timeline(runtime)
  tl.show_progress = True

#### Defining nodes, Alice and Bob (the receivers) and Charlie (the sender)
  Alice = BBM92_receiver(name="Alice")
  Bob = BBM92_receiver(name="Bob")
  Charlie = BBM92_SPDC_source(name="Charlie")



#### Defining channels and connecting the nodes
  qc0 = optch.QuantumChannel("qc0", tl, distance=distance, polarization_fidelity=0.97, attenuation=0.0002)
  qc1 = optch.QuantumChannel("qc1", tl, distance=distance, polarization_fidelity=0.97, attenuation=0.0002)
  qc0.set_ends(Charlie, Alice.name)
  qc1.set_ends(Charlie, Bob.name)

  cc0 = optch.ClassicalChannel("cc0", tl, distance=distance)
  cc1 = optch.ClassicalChannel("cc1", tl, distance=distance)
  cc0.set_ends(Alice, Bob.name)
  cc1.set_ends(Bob, Alice.name)


  for i in range (1000):
    Charlie.emit_photon()

  print("Counter of alice D1 is: ", Alice.get_counts_D1())
  print("Counter of alice D2 is: ", Alice.get_counts_D1())

  print("Counter of Bob D1 is: ", Bob.get_counts_D1())
  print("Counter of Bob D2 is: ", Bob.get_counts_D1())


if __name__ == "__main__":
  main()
