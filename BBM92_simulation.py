
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
from sequence.components.memory import Memory




################# Class that represents the counters for the photon detectors #########################


class det_counter:
  def __init__(self):
    self.count = 0
    self.time = 0
    self.det_time=[]

  def trigger(self, detector, info):
    self.count += 1
    self.time = info['time']
    self.det_time.append(self.time)



################# Class that represents a BBM92 receiver ################################
class BBM92_receiver(nd.Node):
  
  def __init__(self,name="unnamed",det_eff=0.3,dark_counts=5000,count_rate_max=3*10**(6), tl = Timeline(2e10)):
    ###### detector definition ##############
    super().__init__(name, tl)
    self.name=name
    dect1_name= name+ "det1"
    dect2_name= name+ "det2"
    dect3_name= name+ "det3"
    dect4_name= name+ "det4"
    self.det_1= det.Detector(name= dect1_name, timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max)
    self.det_2= det.Detector(name= dect2_name, timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max)
    self.det_3= det.Detector(name= dect3_name, timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max)
    self.det_4= det.Detector(name= dect4_name, timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max)
    self.add_component(self.det_1)
    self.add_component(self.det_2)
    self.add_component(self.det_3)
    self.add_component(self.det_4)

    ######## Beam splitter definition ###############
    bs_name= name+ "initial_beam_splitter"
    self.bs=BS.BeamSplitter(name=bs_name, timeline=tl, fidelity=0.99)
    self.add_component(self.bs)

    ######### Polarization beam splitter definition, measurement on the Z basis ############
    PBSA_name= name + "Polarizing_Beam_splitter_A"
    self.PBSA=BS.BeamSplitter(name=PBSA_name, timeline=tl, fidelity=0.99)
    self.add_component(self.PBSA)
    self.PBSA.set_basis_list(basis_list=[0], start_time=0, frequency=0)

    ######### Polarization beam splitter definition, measurement on the Z basis ############
    PBSB_name= name + "Polarizing_Beam_splitter_B"
    self.PBSB=BS.BeamSplitter(name=PBSB_name, timeline=tl, fidelity=0.99)
    self.add_component(self.PBSB)
    self.PBSB.set_basis_list(basis_list=[1], start_time=0, frequency=0)

    ######## Connecting the hardware ##################
    self.set_first_component(self.bs)
    self.bs.add_receiver(self.PBSA)
    self.bs.add_receiver(self.PBSA)
    self.PBSA.add_receiver(self.det_1)
    self.PBSA.add_receiver(self.det_2)
    self.PBSB.add_receiver(self.det_3)
    self.PBSB.add_receiver(self.det_4)
  
    ######### Adding a counter #############
    self.counter_0= det_counter()
    self.det_1.attach(self.counter_0)
    self.counter_1= det_counter()
    self.det_2.attach(self.counter_1)
    self.counter_plus= det_counter()
    self.det_3.attach(self.counter_plus)
    self.counter_minus= det_counter()
    self.det_4.attach(self.counter_minus)

    ######### Defining list of detected basis ##########
    

    print("BBM92_receiver named: ",  self.name ," has been set")

#  def emission_flag(self):
#    print(self.name, "has been informed of an emission event")


############### class that represents a BBM92 source ###################################

class BBM92_SPDC_source(nd.Node):
  
  def __init__(self,name="unnamed", freq=100000.0, tl = Timeline(2e10), receiver_1=None, receiver_2=None):
  ###### source definition ##############
    super().__init__(name, tl)
    self.name=name
    SPDC_name= name +" SPDC_source"
    self.SPDC_source = ls.SPDCSource(name=SPDC_name, timeline=tl, wavelengths=[1550,1550], frequency= freq, mean_photon_num=0.1, encoding_type={'bases': [((1 + 0j, 0j), (0j, 1 + 0j)), ((0.7071067811865476 + 0j, 0.7071067811865476 + 0j), (-0.7071067811865476 + 0j, 0.7071067811865476 + 0j))], 'name': 'polarization'}, phase_error=0.1)
    self.add_component(self.SPDC_source)
    self.receiver_1=receiver_1
    self.receiver_2=receiver_2
    self.SPDC_source.add_receiver(self.receiver_1)
    self.SPDC_source.add_receiver(self.receiver_2)
    print("BBM92_SPDC_source named: ",  self.name ," has been set")



  def emit_photon(self):
         c = 1 / mt.sqrt(2)
         self.SPDC_source.emit([[c+0j,c+0j]])
#         self.receiver_1.emission_flag()
#         self.receiver_2.emission_flag()
         



##################### main program definition ############################

def main():

  runtime = 4e11
  distance = 1e3
  tl = Timeline(runtime)
  tl.show_progress = True

#### Defining nodes, Alice and Bob (the receivers) and Charlie (the sender)
  Alice = BBM92_receiver(name="Alice", tl=tl)
  Bob = BBM92_receiver(name="Bob", tl=tl)
  Charlie = BBM92_SPDC_source(name="Charlie", tl=tl, receiver_1=Alice, receiver_2=Bob)
  Alice.set_seed(0)
  Bob.set_seed(1)
  Charlie.set_seed(0)



#### Defining channels and connecting the nodes
  qc0 = optch.QuantumChannel("qc0", tl, distance=distance, polarization_fidelity=0.97, attenuation=0.0002)
  qc1 = optch.QuantumChannel("qc1", tl, distance=distance, polarization_fidelity=0.97, attenuation=0.0002)
  qc0.set_ends(Charlie, Alice.name)
  qc1.set_ends(Charlie, Bob.name)

  cc0 = optch.ClassicalChannel("cc0", tl, distance=distance)
  cc1 = optch.ClassicalChannel("cc1", tl, distance=distance)
  cc0.set_ends(Alice, Bob.name)
  cc1.set_ends(Bob, Alice.name)




  process = Process(Charlie, "emit_photon", [])
  emission_event = Event(0, process)
  tl.schedule(emission_event)
  


  tl.init()
  tl.run()

  Alice_dection_times=[]
  Bob_dection_times=[]


  Alice_dection_times.append(Alice.counter_plus.det_time)
  Alice_dection_times.append(Alice.counter_minus.det_time)
  Alice_dection_times.append(Alice.counter_0.det_time )
  Alice_dection_times.append(Alice.counter_1.det_time)

  Bob_dection_times.append(Alice.counter_plus.det_time)
  Bob_dection_times.append(Alice.counter_minus.det_time)
  Bob_dection_times.append(Alice.counter_0.det_time )
  Bob_dection_times.append(Bob.counter_1.det_time)


  for A_det_t in Alice_dection_times:
    for B_det_t in Bob_dection_times:
      if A_det_t==B_det_t:
        print("Siumm" )

if __name__ == "__main__":
  main()
