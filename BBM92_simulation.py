
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
import sequence.utils.log as log



################## Defining register class ##############################


class Register():
   def __init__(self,name=""):
    self.n_of_events=0
    self.name="Register named "+ name
    self.detection_events = {
    "Basis": [],
    "Time": [],
    "Det_result": []
    }


   ######### Functions which is called by the detector counters and  updatwes the event list ################
   def register_event(self, owner_name, time, det_result):
    self.n_of_events=self.n_of_events+1
    self.detection_events["Time"].append(time)
    if (det_result=="+" or det_result=="-"):
      self.detection_events["Basis"].append("+-")
    elif (det_result=="0" or det_result=="1"):
      self.detection_events["Basis"].append("01")

    if (det_result=="+" or det_result=="0"):
      self.detection_events["Det_result"].append(0)
    elif (det_result=="-" or det_result=="1"):
      self.detection_events["Det_result"].append(1)    



   ######### Function which allows to obtain the dictionary with the measurement results #################
   def get_registered_events(self):
    return self.detection_events

      

################# Class that represents the counters for the photon detectors #########################


class det_counter():
  def __init__(self,det_name="none", owner="none", reg=None):
    self.name=det_name
    self.owner=owner
    self.reg=reg
    self.count = 0
    self.time = 0
    self.det_time=[]
    self.name="counter of " + det_name


  ################# Trigger funcion which is called when a photon is detected. It registers its result in the register objecgt which is passed to the node ##############
  def trigger(self, detector, info):
    self.count += 1
    self.time = info['time']
    self.det_time.append(self.time)
    result="none"

    if (self.name=="counter of det_zero"):
      result="0"
    elif (self.name=="counter of det_one"):
      result="1"
    elif (self.name=="counter of det_plus"):
      result="+"
    elif (self.name=="counter of det_minus"):
      result="-"

    self.reg.register_event(self.owner, self.time, result )




################# Class that represents a BBM92 receiver ################################
class BBM92_receiver(nd.Node):
  
  def __init__(self,name="unnamed",det_eff=0.3,dark_counts=5000,count_rate_max=3*10**(6), tl = Timeline(2e10), reg=None):
    ###### detector definition ##############
    super().__init__(name, tl)
    self.reg=reg
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
    self.counter_0= det_counter(det_name="det_zero", owner=self.name, reg=self.reg)
    self.det_1.attach(self.counter_0)
    self.counter_1= det_counter(det_name="det_one", owner=self.name, reg=self.reg)
    self.det_2.attach(self.counter_1)
    self.counter_plus= det_counter (det_name="det_plus", owner=self.name, reg=self.reg)
    self.det_3.attach(self.counter_plus)
    self.counter_minus= det_counter(det_name="det_minus", owner=self.name, reg=self.reg)
    self.det_4.attach(self.counter_minus)

    ######### Informing user about creation of the BBM92_receiver  ##########
    print("BBM92_receiver named: ",  self.name ," has been set")



############### class that represents a BBM92 source ###################################

class BBM92_SPDC_source(nd.Node):

  def __init__(self,name="unnamed", freq=10000000.0, tl = Timeline(2e10), receiver_1=None, receiver_2=None):
    super().__init__(name, tl)
    self.frequency=freq
    self.name=name
    SPDC_name= name +" SPDC_source"
    encoding_polarizzazione = {'name': 'polarization', 'bases': [ ((1 + 0j, 0j), (0j, 1 + 0j)),((0.707+0j, 0.707+0j), (-0.707+0j, 0.707+0j)) ]}
    self.SPDC_source = ls.SPDCSource(name=SPDC_name, timeline=tl, frequency=freq, wavelengths=[1550,1550], mean_photon_num=0.1, encoding_type=encoding_polarizzazione, phase_error=0.1)


#encoding_type={'bases': [((1 + 0j, 0j), (0j, 1 + 0j)), ((0.7071067811865476 + 0j, 0.7071067811865476 + 0j), (-0.7071067811865476 + 0j, 0.7071067811865476 + 0j))], 'name': 'polarization'}

    self.add_component(self.SPDC_source)
    self.receiver_1=receiver_1
    self.receiver_2=receiver_2
    self.SPDC_source.add_receiver(self.receiver_1)
    self.SPDC_source.add_receiver(self.receiver_2)

    ######### Informing user about creation of the BBM92_SPDC_emitter  ##########
    print("BBM92_SPDC_source named: ",  self.name ," has been set")


  def emit_photon(self):
         c = 1 / mt.sqrt(2)
         self.SPDC_source.emit([[c + 0j,0 +0j,0 +0j, c + 0j]])

  def get_frequency(self):
    return self.frequency

         



##################### main program definition ############################

def main():

  Reg_Alice=Register("BBM92_Reg_Alice")
  Reg_Bob=Register("BBM92_Reg_Bob")
  runtime = 1e12
  sim_time=runtime*10**(-12)
  distance = 1e3
  tl = Timeline(runtime)
  tl.show_progress = True

#### Defining nodes, Alice and Bob (the receivers) and Charlie (the sender)
  Alice = BBM92_receiver(name="Alice", tl=tl, reg=Reg_Alice)
  Bob = BBM92_receiver(name="Bob", tl=tl, reg=Reg_Bob)
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



  ############## Defining the processes and events ############# 
  process = Process(Charlie, "emit_photon", [])
  emission_event = Event(0, process)
  tl.schedule(emission_event)



  ###### Running the simulation ###############
  tl.init()
  tl.run()

  #######  Getting the dictionary with the result ###########
  results_Alice=Reg_Alice.get_registered_events()
  results_Bob=Reg_Bob.get_registered_events()

  n_Alice_measurements=len(results_Alice["Time"])
  n_Bob_measurements=len(results_Bob["Time"])
  offset=3000
  coincident_counts=0 
  Alice_raw_key=[]
  Bob_raw_key=[]


  print("\n\n")
  print("Number of Alice counts is: ",n_Alice_measurements )
  print("Number of Bob counts is: ",n_Bob_measurements )

  for i in range(n_Alice_measurements):
    for j in range(n_Bob_measurements):
      if np.abs(results_Alice["Time"][i]-results_Bob["Time"][j]-offset)<=300000 and results_Alice["Basis"][i]==results_Bob["Basis"][j]:
        coincident_counts=coincident_counts+1
        Alice_raw_key.append(results_Alice["Det_result"][i])
        Bob_raw_key.append(results_Bob["Det_result"][j])
        print("ess", i,j, "time difference= ", results_Alice["Time"][i]-results_Bob["Time"][j]) 
        break



  print("Number of coincident counts is: ",coincident_counts )


  correct_counts=0
  different_counts=0
  for k in range(coincident_counts)  :
    if Bob_raw_key[k]==Alice_raw_key[k]:
      correct_counts=correct_counts+1
    else:
      different_counts=different_counts+1



  QBER=different_counts/coincident_counts
  print("The QBER is: ", QBER)



if __name__ == "__main__":
  main()
