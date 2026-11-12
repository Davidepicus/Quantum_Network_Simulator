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
from sequence.components.photon import Photon
from scipy import constants





### redefining Fock Beam splitter class ####

class BeamSplitter(BS.FockBeamSplitter2):
   def __init__(self, name_i="none", owner_i=None, tml=None, eff=1, ph_counter=0, source_list=["none"]):
    super().__init__( name=name_i, owner=owner_i, timeline=tml, efficiency=eff, photon_counter=ph_counter, src_list=source_list)

     ######## get funcion of the beam splitter which overrides the one present in the Sequence files ##########
   def get(self, photon):
        """Receive photon from one end nodes"""
#        print("photon in the Beam splitter", self.name)
        self.photon_counter += 1
            
        #print(f"Photon counter BEAM SPLITTER: {self.photon_counter} at time {self.timeline.now()}")

        selected_receiver = random.choice(self._receivers)
        #print(selected_receiver)

        if self.photon_counter == 1:
            
            selected_receiver.get(photon)
            #selected_receiver.get_2(photon)
            self.photon_counter=0 


    ######function that receives qubit from a channel and calls the get method for the beam splitter###############
#   def receive_qubit(self, source_name, qubit):
#      self.get(source=source_name, photon=qubit)

################# redifining correct beam splitter class ###################


class PBeamSplitter(BS.BeamSplitter):
  def __init__(self, name="none", timeline=None, fidelity=1):
    super().__init__(name=name, timeline=timeline, fidelity=fidelity)
    self.con=0

  def get_2(self, photon):
    self.con=self.con+1
    #print("get_2 called")
    #self.get(photon)

  def getx2(self, photon):
    self.con=self.con+1
    print("getx2 called")
    #self.get(photon)

  def get_2x2(self, photon):
    self.con=self.con+1
    print("get_2x2 called")
    #self.get(photon)


  ######## get funcion of the  polarizing beam splitter which overrides the one present in the Sequence files ##########
  def get(self, photon):
#    print("photon in the  Polarizing Beam splitter", self.name)
    """Method to receive a photon for measurement.

        Args:
            photon (Photon): photon to measure (must have polarization encoding)

        Side Effects:
            May call get method of one receiver.
        """

    assert photon.encoding_type["name"] == "polarization", "Beamsplitter should only be used with polarization."

    if self.get_generator().random() < self.fidelity: 
      #index = int((self.timeline.now() - self.start_time) * self.frequency * 1e-12)
      index=1
      if 0 > index or index >= len(self.basis_list):
          return
      else:
        res = photon.measure(photon.encoding_type["bases"][self.basis_list[index]],photon, self.get_generator())
        #print(self._receivers)
        self._receivers[res].get(photon)




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

      
###### defining counter class #####
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
#    print("CLICK", info["time"])

    if (self.name=="counter of det_zero"):
      result="0"
    elif (self.name=="counter of det_one"):
      result="1"
    elif (self.name=="counter of det_plus"):
      result="+"
    elif (self.name=="counter of det_minus"):
      result="-"

    self.reg.register_event(self.owner, self.time, result )




################ Defining quantum channels #########################

class Qchannel(optch.QuantumChannel):
  def __init__(self, name="none", timeline=None, attenuation=0, distance=1*10**(3), polarization_fidelity= 1.0, light_speed= 3*10**(-4), frequency= 8e7):
    super().__init__( name=name, timeline=timeline, attenuation=attenuation, distance=distance, polarization_fidelity= polarization_fidelity, light_speed=light_speed, frequency= frequency)
    self.sender_node=None
    self.receiver_node=None
    self.frequency=8e12


  def get(self, photon):
        source=self.sender_node
        """Method to transmit photon-encoded qubits.

        Args:
            qubit (Photon): photon to be transmitted.

        Side Effects:
            Receiver node may receive the qubit (via the `receive_qubit` method).
        """
        
#        log.logger.info("{} send qubit with state {} to {} by Channel {}".format(self.sender.name, qubit.quantum_state, self.receiver, self.name))

        assert self.delay >= 0 and self.loss <= 1, f"QuantumChannel init() function has not been run for {self.name}"
#        assert source == self.sender

        # remove lowest time bin
        if len(self.send_bins) > 0:
            time = -1
            while time < self.timeline.now():
                time_bin = hq.heappop(self.send_bins)
                time = self.timebin_to_time(time_bin, self.frequency)
            assert time == self.timeline.now(), f"qc {self.name} transmit method called at invalid time"

        # check if photon state using Fock representation
        if photon.encoding_type["name"] == "fock":
            key = photon.quantum_state  # if using Fock representation, the `quantum_state` field is the state key.
            # apply loss channel on photonic statex
            self.timeline.quantum_manager.add_loss(key, self.loss)

            # schedule receiving node to receive photon at future time determined by light speed
            future_time = self.timeline.now() + self.delay
            process = Process(self.receiver, "receive_qubit", [source.name, photon])
            event = Event(future_time, process)
            self.timeline.schedule(event)

        # if not using Fock representation, check if photon kept
        elif (self.sender.get_generator().random() > self.loss) or photon.is_null:
            if self._receiver_on_other_tl():
                self.timeline.quantum_manager.move_manage_to_server(photon.quantum_state)

            if photon.is_null:
                photon.add_loss(self.loss)

            # check if polarization encoding and apply necessary noise
            if photon.encoding_type["name"] == "polarization" and self.sender.get_generator().random() > self.polarization_fidelity:
#                print("phass")
                photon.random_noise(self.get_generator())
            
            # schedule receiving node to receive photon at future time determined by light speed
            future_time = self.timeline.now() + self.delay
            process = Process(self.receiver, "receive_qubit", [source.name, photon])
            event = Event(future_time, process)
            self.timeline.schedule(event)

        # if not using Fock representation, if photon lost, exit
        else:
            pass
 

  ######## function which connects the nodes to the channel ############
  def connect(self, sender_node=None, receiver_node=None ):
    print(sender_node.name, "Has been connected with ", receiver_node.name)
    self.sender_node=sender_node
    self.receiver_node=receiver_node
    self.set_ends(sender_node, receiver_node.name)




  ##### function to update distance and attenuation of the channels #######
  def update_channel(self, dist=-1, absor=-1):
    if dist<0 or absor<0 :
      print("Error, either distance value or absorption are invalid, parameters have not been changed")
    else:
      self.distance=dist
      self.attenuation=absor
      self.delay = round(self.distance / self.light_speed)
      self.loss = 1 - 10 ** (self.distance * self.attenuation / -10)
#      print("attenuation has been set to ", self.attenuation)





################### Defining class for a BB84 sender ####################

class BB84_laser_source(nd.Node):

  def __init__(self,name="unnamed", freq=1*10**(8), tl = Timeline(2e10), receiver=None, av_photons=0.1, phase_err=0, reg=None):
    super().__init__(name, tl)
    self.reg=reg
    self.frequency=freq
    self.name=name
    laser_name= name +" laser_source"
    polarization ={"name": "polarization","bases": [((complex(1), complex(0)), (complex(0), complex(1))),((complex(np.sqrt(1 / 2)), complex(np.sqrt(1 / 2))), (complex(-np.sqrt(1 / 2)), complex(np.sqrt(1 / 2))))]}


    self.laser_source=ls.LightSource(name=laser_name, timeline=tl,frequency=self.frequency, wavelength=1550, mean_photon_num=av_photons,phase_error=phase_err)


    ###### attaching the SPDC_source to the node ########
    self.add_component(self.laser_source)

    ######## attaching the source to the receiver ############
    self.receiver=receiver
    self.laser_source.add_receiver(self.receiver)


    ######### Informing user about creation of the BB84_laser_emitter  ##########
    print("BB84_laser_source named: ",  self.name ," has been set")


  def emit_photons(self):
    ### defining qubits to encode and to emit ####
    c = 1 / mt.sqrt(2)
    phi_plus = (c + 0j, c + 0j)
    phi_minus = (c + 0j, -c + 0j)
    zero= (1+0j,  0j)
    one=(0j, 1 + 0j)



    #### defining which qubit to encode and  calling the register#####
    randindex= random.randint(0, 3)
    if randindex==0:
      state_list=zero
      self.reg.register_event( owner_name=self.name, time=self.timeline.now(), det_result="0")
    elif randindex==1:
      state_list=one
      self.reg.register_event( owner_name=self.name, time=self.timeline.now(), det_result="1")
    elif randindex==2:
      state_list=phi_plus
      self.reg.register_event( owner_name=self.name, time=self.timeline.now(), det_result="+")
    elif randindex==3:
      state_list=phi_minus
      self.reg.register_event( owner_name=self.name, time=self.timeline.now(), det_result="-")

   ### emitting photons ####
#    print("The index of the chosen state is:  ", randindex )
    self.laser_source.emit( [state_list])


#### defining class for a BB84 receiver ######

class BB84_receiver(nd.Node):
  
  def __init__(self,name="unnamed",det_eff=0.35,dark_counts=5000,count_rate_max=3*10**(18),time_resolution= 150, tl = Timeline(2e10), reg=None, PBS_efficiency=0.99, BS_efficiency=0.99):
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
    self.bs=BeamSplitter(name_i=bs_name, owner_i=self.name,  tml=tl, eff=BS_efficiency, ph_counter=self.bs_counter, source_list=["Alice"])
    self.add_component(self.bs)

    ######### Polarization beam splitter definition, measurement on the Z basis ############
    PBSA_name= name + "Polarizing_Beam_splitter_A"
    self.PBSA=PBeamSplitter(name=PBSA_name, timeline=tl, fidelity=PBS_efficiency)
    self.add_component(self.PBSA)
    self.PBSA.set_basis_list([0,0], start_time=0, frequency=1*10**12)

    ######### Polarization beam splitter definition, measurement on the x basis ############
    PBSB_name= name + "Polarizing_Beam_splitter_B"
    self.PBSB=PBeamSplitter(name=PBSB_name, timeline=tl, fidelity=PBS_efficiency)
    self.add_component(self.PBSB)
    self.PBSB.set_basis_list([1,1], start_time=0, frequency=1*10**12)




    ######## Connecting the hardware ##################
    self.set_first_component(self.bs.name)
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


  ######redefining get function such that it can obtain photons from other objects ######################## 
  def get(self, photon):
    self.bs.get( photon)




def main():

  with open('vonfig.json', 'r', encoding='utf-8') as file:
    parameters= json.load(file)


#### Parameters settings ############

  runtime =parameters["general_parameters"]["runtime"]*(10)**(12)
  emission_frequency=parameters["emitter_parameters"]["emitter_frequency"]*10**(-12) ###in numbers of emission per picoseconds
  average_photons_emitted=parameters["emitter_parameters"]["av_ph_emitted"]
  detector_efficiency=parameters["Detector_parameters"]["det_eff"]
  dark_counts=parameters["Detector_parameters"]["dark_counts"]
  count_rate_max=parameters["Detector_parameters"]["count_rate_max"]
  time_resolution=parameters["Detector_parameters"]["time_resolution"]
  PBS_eff=parameters["Detector_parameters"]["PBS_efficiency"]
  BS_eff=parameters["Detector_parameters"]["Beam_splitter_efficiency"]
  LS=constants.speed_of_light*parameters["channel_parameters"]["light_speed"]*10**(-12)
  pol_fidelity=parameters["channel_parameters"]["Polarization_error"]


#  emission_frequency=80000*10**(-12) ###in numbers of emission per picoseconds
#  average_photons_emitted=0.01
#  detector_efficiency=0.5
#  dark_counts=500
#  count_rate_max=8*10**(33)
#  time_resolution=5
#  PBS_eff=0.99
#  BS_eff=0.99
#  runtime=2*10**(13)
  distance_AB= 2*10**(3)
  ch_attenuationAB=0.002
#  LS=constants.speed_of_light*10**(-12)
  phase_error=0
#  pol_fidelity=0.96

##### Defining timeline ################
  tl = Timeline(runtime)
  tl.show_progress = True


  ######### defining registers for Alice and Bob 
  Reg_Alice=Register("BB84_Reg_Alice")
  Reg_Bob=Register("BB84_Reg_Bob")


#### Defining nodes, Alice (the sender), Bob (the receiver), and quantum channels ################
  qc0 = Qchannel("qc0", tl, distance=distance_AB, polarization_fidelity=pol_fidelity, attenuation=ch_attenuationAB,light_speed=LS)
  Alice = BB84_laser_source(name="Alice", tl=tl, receiver=qc0, av_photons=average_photons_emitted, phase_err=phase_error, reg=Reg_Alice)
  Bob = BB84_receiver(name="Bob",  tl=tl,det_eff=detector_efficiency, dark_counts=dark_counts ,count_rate_max=count_rate_max ,time_resolution=time_resolution, reg=Reg_Bob,PBS_efficiency=PBS_eff, BS_efficiency=BS_eff)
  qc0.connect(Alice,Bob)


  
  ############## Defining the emission processes and setting events #############
  number_of_emissions=int(emission_frequency*runtime)-1
  emission_period=(int(1/emission_frequency))
  process = Process(Alice, "emit_photons", [])
  emission_event=[]
  for t in range(number_of_emissions):
    process = Process(Alice, "emit_photons", [])
    event = Event(emission_period*t, process)
    tl.schedule(event)


  ###### Running the simulation ###############
  tl.init()
  tl.run()


  #######  Getting the dictionary with the result ###########
  results_Alice=Reg_Alice.get_registered_events()
  results_Bob=Reg_Bob.get_registered_events()


  n_Alice_measurements=len(results_Alice["Time"])
  n_Bob_measurements=len(results_Bob["Time"])
  offset=(distance_AB)/LS
  coincident_counts=0 
  Alice_raw_key=[]
  Bob_raw_key=[]


  print("\n\n")
  print("Number of Alice detections is: ",n_Alice_measurements )
  print("Number of Bob detections is: ",n_Bob_measurements )


#  print(results_Alice["Basis"][:30])
#  print(results_Bob["Basis"][:30])
#  print(offset)

  print("\n","### Starting to post processing analysis ###")


  start=0
  for i in range (n_Alice_measurements): 
    for j in range(start,n_Bob_measurements,1):
#      print(results_Alice["Time"][i], results_Bob["Time"][j])
      if np.abs(results_Alice["Time"][i]-results_Bob["Time"][j]+offset)<=3*time_resolution and results_Alice["Basis"][i]==results_Bob["Basis"][j]:
        coincident_counts=coincident_counts+1
        Alice_raw_key.append(results_Alice["Det_result"][i])
        Bob_raw_key.append(results_Bob["Det_result"][j])
        start=j
        break
      elif np.abs(results_Alice["Time"][i]-results_Bob["Time"][j]-offset)> time_resolution*10**(10):
        break
      
    

  print("Total length of the raw key is: ",coincident_counts )
  key_rate=int(coincident_counts/(runtime*10**(-12)))
  print("Key transmission rate is: " ,key_rate, " Qubit/sec")

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
