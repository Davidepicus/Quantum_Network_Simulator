
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
import random
import json



################ Redefining quantum channel ################################

#optch.QuantumChannel("qc0", tl, distance=distance, polarization_fidelity=0.99, attenuation=0.0002)

class Qchannel(optch.QuantumChannel):
  def __init__(self, name="none", timeline=None, attenuation=0, distance=1*10**(3), polarization_fidelity= 1.0, light_speed= 3*10**(-4), frequency= 8e7):
    super().__init__( name=name, timeline=timeline, attenuation=attenuation, distance=distance, polarization_fidelity= polarization_fidelity, light_speed=light_speed, frequency= frequency)



################ redefining correct Polarizing beam splitter class ###########################

class BeamSplitter(BS.FockBeamSplitter2):
   def __init__(self, name_i="none", owner_i=None, tml=None, eff=1, ph_counter=0, source_list=["none"]):
    super().__init__( name=name_i, owner=owner_i, timeline=tml, efficiency=eff, photon_counter=ph_counter, src_list=source_list)


   def get(self, source, photon) -> None:
        """Receive photon from one end nodes"""

        self.photon_counter += 1
            
        #print(f"Photon counter BEAM SPLITTER: {self.photon_counter} at time {self.timeline.now()}")

        selected_receiver = random.choice(self._receivers)
        #print(selected_receiver)

        if self.photon_counter == 1:
            
            selected_receiver.get(photon)
            #selected_receiver.get_2(photon)
            self.photon_counter=0 

#        elif self.photon_counter == 2:
            
#            self._receivers[0].photon_counter = 0
#            self._receivers[1].photon_counter = 0
#            self._receivers[0].photon_counter2 = 0
#            self._receivers[1].photon_counter2 = 0
#            selected_receiver.getx2(photon)
#            selected_receiver.get_2x2(photon)


################# redifining correct beam splitter class ###################


class PBeamSplitter(BS.BeamSplitter):
  def __init__(self, name="none", timeline=None, fidelity=1):
    super().__init__(name=name, timeline=timeline, fidelity=fidelity)
    self.con=0

  def get_2(self, photon):
    #self.con=self.con+1
    #print("get_2 called")
    self.get(photon)

  def getx2(self, photon):
    self.con=self.con+1
    print("getx2 called")
    #self.get(photon)

  def get_2x2(self, photon):
    self.con=self.con+1
    print("get_2x2 called")
    #self.get(photon)

  def get(self, Photon, **kwargs) -> None:
    """Method to receive a photon for measurement.

        Args:
            photon (Photon): photon to measure (must have polarization encoding)

        Side Effects:
            May call get method of one receiver.
        """

    assert Photon.encoding_type["name"] == "polarization", "Beamsplitter should only be used with polarization."

    if self.get_generator().random() < self.fidelity: 
      #index = int((self.timeline.now() - self.start_time) * self.frequency * 1e-12)
      index=1
      if 0 > index or index >= len(self.basis_list):
          return
      else:
        res = Photon.measure(Photon.encoding_type["bases"][self.basis_list[index]],Photon, self.get_generator())
        #print(self._receivers)
        self._receivers[res].get(Photon)



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
  def __init__(self,det_name="", owner="none", reg=None, timeline=None):
    self.timeline=timeline
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
    self.time = info['time'] #self.timeline.now()
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
  
  def __init__(self,name="unnamed",det_eff=0.35,dark_counts=5000,count_rate_max=3*10**(18),time_resolution= 150, tl = Timeline(2e10), reg=None):
    ###### detector definition ##############
    super().__init__(name, tl)
    self.reg=reg
    self.name=name
    dect1_name= name+ "det1"
    dect2_name= name+ "det2"
    dect3_name= name+ "det3"
    dect4_name= name+ "det4"
    self.det_1= det.Detector(name= dect1_name, timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max, time_resolution=time_resolution)
    self.det_2= det.Detector(name= dect2_name, timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max, time_resolution=time_resolution)
    self.det_3= det.Detector(name= dect3_name, timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max, time_resolution=time_resolution)
    self.det_4= det.Detector(name= dect4_name, timeline=tl, efficiency=det_eff, dark_count=dark_counts, count_rate=count_rate_max, time_resolution=time_resolution)
    self.add_component(self.det_1)
    self.add_component(self.det_2)
    self.add_component(self.det_3)
    self.add_component(self.det_4)

    ######## Beam splitter definition ###############
    self.bs_counter=0
    bs_name= name+ "initial_beam_splitter"
    self.bs=BeamSplitter(name_i=bs_name, owner_i=self.name,  tml=tl, eff=0.99, ph_counter=self.bs_counter, source_list=["Charlie"])
    self.add_component(self.bs)

    ######### Polarization beam splitter definition, measurement on the Z basis ############
    # polarization_Z  =[complex(1), complex(0), complex(0), complex(1)]
    PBSA_name= name + "Polarizing_Beam_splitter_A"
    self.PBSA=PBeamSplitter(name=PBSA_name, timeline=tl, fidelity=0.99)
    self.add_component(self.PBSA)

    

    self.PBSA.set_basis_list([0,0], start_time=0, frequency=1*10**12)

    ######### Polarization beam splitter definition, measurement on the x basis ############
    #polarization_X = [complex(np.sqrt(1 / 2)), complex(np.sqrt(1 / 2)), complex(-np.sqrt(1 / 2)), complex(np.sqrt(1 / 2))] 
    PBSB_name= name + "Polarizing_Beam_splitter_B"
    self.PBSB=PBeamSplitter(name=PBSB_name, timeline=tl, fidelity=0.99)
    self.add_component(self.PBSB)


    self.PBSB.set_basis_list([1,1], start_time=0, frequency=1*10**12)
    ######## Connecting the hardware ##################
    self.set_first_component(self.bs)
    self.bs.add_receiver(self.PBSA)
    self.bs.add_receiver(self.PBSB)
    self.PBSA.add_receiver(self.det_1)
    self.PBSA.add_receiver(self.det_2)
    self.PBSB.add_receiver(self.det_3)
    self.PBSB.add_receiver(self.det_4)
  
    ######### Adding a counter #############
    self.counter_0= det_counter(det_name="det_zero", owner=self.name, reg=self.reg, timeline=tl)
    self.det_1.attach(self.counter_0)
    self.counter_1= det_counter(det_name="det_one", owner=self.name, reg=self.reg, timeline=tl)
    self.det_2.attach(self.counter_1)
    self.counter_plus= det_counter (det_name="det_plus", owner=self.name, reg=self.reg, timeline=tl)
    self.det_3.attach(self.counter_plus)
    self.counter_minus= det_counter(det_name="det_minus", owner=self.name, reg=self.reg, timeline=tl)
    self.det_4.attach(self.counter_minus)

    ######### Informing user about creation of the BBM92_receiver  ##########
    print("BBM92_receiver named: ",  self.name ," has been set")

  def get(self, photon):
    self.bs.get("Charlie", photon)

############### class that represents a BBM92 source ###################################

class BBM92_SPDC_source(nd.Node):

  def __init__(self,name="unnamed", freq=1*10**(8), tl = Timeline(2e10), receiver_1=None, receiver_2=None):
    super().__init__(name, tl)
    self.frequency=freq
    self.name=name
    SPDC_name= name +" SPDC_source"
    polarization ={"name": "polarization","bases": [((complex(1), complex(0)), (complex(0), complex(1))),((complex(np.sqrt(1 / 2)), complex(np.sqrt(1 / 2))), (complex(-np.sqrt(1 / 2)), complex(np.sqrt(1 / 2))))]}

    self.SPDC_source = ls.SPDCSource(name=SPDC_name, timeline=tl, frequency=freq, wavelengths=[1550,1550], mean_photon_num=0.01,encoding_type=polarization, phase_error=0)


#encoding_type={'bases': [((1 + 0j, 0j), (0j, 1 + 0j)), ((0.7071067811865476 + 0j, 0.7071067811865476 + 0j), (-0.7071067811865476 + 0j, 0.7071067811865476 + 0j))], 'name': 'polarization'}

    self.add_component(self.SPDC_source)
    self.receiver_1=receiver_1
    self.receiver_2=receiver_2
    self.SPDC_source.add_receiver(self.receiver_1)
    self.SPDC_source.add_receiver(self.receiver_2)

    ######### Informing user about creation of the BBM92_SPDC_emitter  ##########
    print("BBM92_SPDC_source named: ",  self.name ," has been set")


  def emit_photons(self):
    c = 1 / mt.sqrt(2)
    phi_plus = [c + 0j, c + 0j]
    self.SPDC_source.emit( [phi_plus])


  def get_frequency(self):
    return self.frequency

         



##################### main program definition ############################

def main():
  ######### importing parameters from config.json #######
  with open('config.json', 'r', encoding='utf-8') as file:
    parameters= json.load(file)

  ######### defining registers for Alice and Bob 
  Reg_Alice=Register("BBM92_Reg_Alice")
  Reg_Bob=Register("BBM92_Reg_Bob")

  runtime =parameters["general_parameters"]["runtime"] #*(10)**(10)
  emission_frequency=parameters["emitter_parameters"]["emitter_frequency"] ###30*10**(-6) ###in numbers of emission per picoseconds
  number_of_emissions=int(emission_frequency*runtime)-1
  emission_period=(int(1/emission_frequency))
  sim_time=runtime*10**(-12)
  distance = 1e3
  tl = Timeline(runtime)
  tl.show_progress = True

#### Defining nodes, Alice and Bob (the receivers) and Charlie (the sender), and quantum channels
  Alice = BBM92_receiver(name="Alice", tl=tl,det_eff=parameters["Detector_parameters"]["det_eff"],dark_counts=parameters["Detector_parameters"]["dark_counts"],count_rate_max=parameters["Detector_parameters"]["count_rate_max"],time_resolution=parameters["Detector_parameters"]["time_resolution"],reg=Reg_Alice)
  Bob = BBM92_receiver(name="Bob", tl=tl, reg=Reg_Bob)

  qc0 = optch.QuantumChannel("qc0", tl, distance=distance, polarization_fidelity=0.99, attenuation=0.0002)
  qc1 = optch.QuantumChannel("qc1", tl, distance=distance, polarization_fidelity=0.99, attenuation=0.0002)
  Charlie = BBM92_SPDC_source(name="Charlie", tl=tl, receiver_1=Alice, receiver_2=Bob)
  Alice.set_seed(0)
  Bob.set_seed(1)
  Charlie.set_seed(2)



#### Setting Bob and Alice as the receivers

  qc0.set_ends(Charlie, Alice.name)
  qc1.set_ends(Charlie, Bob.name)

#  cc0 = optch.ClassicalChannel("cc0", tl, distance=distance)
#  cc1 = optch.ClassicalChannel("cc1", tl, distance=distance)
#  cc0.set_ends(Alice, Bob.name)
#  cc1.set_ends(Bob, Alice.name)



  ############## Defining the processes and events #############
  process = Process(Charlie, "emit_photons", [])
  emission_event=[]
  for t in range(number_of_emissions):
    emission_event.append(Event(emission_period*t, process))

  for h in range (len(emission_event)):
    tl.schedule(emission_event[h])



  ###### Running the simulation ###############
  tl.init()
  tl.run()


  #######  Getting the dictionary with the result ###########
  results_Alice=Reg_Alice.get_registered_events()
  results_Bob=Reg_Bob.get_registered_events()


  n_Alice_measurements=len(results_Alice["Time"])
  n_Bob_measurements=len(results_Bob["Time"])
  offset=0
  coincident_counts=0 
  Alice_raw_key=[]
  Bob_raw_key=[]


  print("\n\n")
  print("Number of Alice detections is: ",n_Alice_measurements )
  print("Number of Bob detections is: ",n_Bob_measurements )

  print("### Starting to post processing analysis ###")

  start=0
  for i in range (n_Alice_measurements): 
    for j in range(start,n_Bob_measurements,1):
      if np.abs(results_Alice["Time"][i]-results_Bob["Time"][j]-offset)<=parameters["Detector_parameters"]["time_resolution"] and results_Alice["Basis"][i]==results_Bob["Basis"][j]:
        coincident_counts=coincident_counts+1
        Alice_raw_key.append(results_Alice["Det_result"][i])
        Bob_raw_key.append(results_Bob["Det_result"][j])
        #print("ess", i,j, "time difference= ", results_Alice["Time"][i]-results_Bob["Time"][j])
        start=j
        break

      
    

  print("Number of coincident counts is: ",coincident_counts )
  key_rate=int(coincident_counts/(runtime*10**(-12)))
  print("Key transmission rate is: " ,key_rate, "Hz")

  correct_counts=0
  different_counts=0
  for k in range(coincident_counts)  :
    if Bob_raw_key[k]==Alice_raw_key[k]:
      correct_counts=correct_counts+1
    else:
      different_counts=different_counts+1

  if coincident_counts==0:
    print("No coincident counts have been detected, therefore QBER cannot be estimated")
  else:
    QBER=(different_counts/coincident_counts)
    print(f"The QBER is: {QBER:.3f}")
  



if __name__ == "__main__":
  main()
