
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
import pandas as pd
from datetime import datetime




################ Redefining quantum channel ################################

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
                print("phass")
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




  def reset_distance(self, dist=-1):
    if dist<0:
      print("Error, distance value is invalid, distance has not been changed")
    else:
      self.distance=dist


  def reset_absorption(self, absor=-1):
    if absor<0:
      print("Error, absorption value is invalid, absorption has not been changed")
    else:
      self.attenuation=absor




################ redefining correct Polarizing beam splitter class ###########################

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
    self.bs=BeamSplitter(name_i=bs_name, owner_i=self.name,  tml=tl, eff=BS_efficiency, ph_counter=self.bs_counter, source_list=["Charlie"])
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

############### class that represents a BBM92 source ###################################

class BBM92_SPDC_source(nd.Node):

  def __init__(self,name="unnamed", freq=1*10**(8), tl = Timeline(2e10), receiver_1=None, receiver_2=None, av_pairs=0.1):
    super().__init__(name, tl)
    self.frequency=freq
    self.name=name
    SPDC_name= name +" SPDC_source"
    polarization ={"name": "polarization","bases": [((complex(1), complex(0)), (complex(0), complex(1))),((complex(np.sqrt(1 / 2)), complex(np.sqrt(1 / 2))), (complex(-np.sqrt(1 / 2)), complex(np.sqrt(1 / 2))))]}

    self.SPDC_source = ls.SPDCSource(name=SPDC_name, timeline=tl, frequency=freq, wavelengths=[1550,1550], mean_photon_num=av_pairs,encoding_type=polarization, phase_error=0)


    ###### attaching the SPDC_source to the node ########
    self.add_component(self.SPDC_source)

    ####### attaching the SPDC source to the receiver   #######################
    self.receiver_1=receiver_1
    self.receiver_2=receiver_2
    self.SPDC_source.add_receiver(self.receiver_1)
    self.SPDC_source.add_receiver(self.receiver_2)

    ######### Informing user about creation of the BBM92_SPDC_emitter  ##########
    print("BBM92_SPDC_source named: ",  self.name ," has been set")

  ######### function that calls the SPDC source and tells it to excite the SPDC crystal ###################
  def emit_photons(self):
        c = 1 / mt.sqrt(2)
        phi_plus = [c + 0j, c + 0j]
        state_list=phi_plus
        self.SPDC_source.emit( [phi_plus])
        


##################### main program definition ############################

def main():
  ######### importing parameters from config.json #######
  with open('config.json', 'r', encoding='utf-8') as file:
    parameters= json.load(file)

  ######### defining registers for Alice and Bob 
  Reg_Alice=Register("BBM92_Reg_Alice")
  Reg_Bob=Register("BBM92_Reg_Bob")

  """
  ##### importing orbits ###########
  # Carica i CSV
  df_lijiang = pd.read_csv("fucino_satellite.csv")
  df_delingha = pd.read_csv("lario_satellite.csv")
  
  df_lijiang["Source"] = "Lijiang"
  df_delingha["Source"] = "Delingha"

  # Unisci i due DataFrame per calcolare il tempo globale
  df_all = pd.concat([df_lijiang, df_delingha], ignore_index=True)

  # Converte Date+Time in datetime assoluto
  def to_datetime(row):
      return datetime.strptime(
          row["Date"] + " " + row["Time"],
          "%d %b %Y %H:%M:%S.%f"
      )

  df_all["Datetime"] = df_all.apply(to_datetime, axis=1)

  # Trova il primo timestamp assoluto
  t0 = df_all["Datetime"].min()
  
  # Calcola i secondi dall'inizio della prima giornata
  df_all["SecondsFromStart"] = (df_all["Datetime"] - t0).dt.total_seconds()


  # Ora ri-separa i DataFrame sincronizzati
  df_lijiang = df_all[df_all["Source"] == "Lijiang"].copy()
  df_delingha = df_all[df_all["Source"] == "Delingha"].copy()


  def atmospheric_attenuation_1550(elev_deg, alpha0_db=0.1):
#    if elev_deg <= 0:
#        return np.inf  # satellite sotto l’orizzonte
    elev_rad = np.deg2rad(np.abs(elev_deg))
    airmass = 1.0 / np.sin(elev_rad)
    return alpha0_db * airmass



   ### defining normal dictionaries from Pandas Dataframes ##
  df_lijiang = df_all.merge(df_lijiang[["Date","Time"]], on=["Date","Time"])
  df_delingha = df_all.merge(df_delingha[["Date","Time"]], on=["Date","Time"])

  # Applica il modello a entrambe le stazioni
  df_lijiang["AtmosphericLoss_dB"] = df_lijiang["Elevation"].apply(atmospheric_attenuation_1550)
  df_delingha["AtmosphericLoss_dB"] = df_delingha["Elevation"].apply(atmospheric_attenuation_1550)
  


  Alice_orb_dict = df_lijiang.to_dict(orient="list")
  Bob_orb_dict = df_delingha.to_dict(orient="list")

#  print(df_delingha.head(5))
#  print(df_delingha["Time"].iloc[0])
#  print(df_delingha["Time"].iloc[1])

  ### rewriting the ranges in meters instead of kilometers ################
  Alice_orb_dict["Range"]=Alice_orb_dict["Range"]*10**(3)
  Alice_orb_dict["Range"]=Alice_orb_dict["Range"]*10**(3)
  """

  # Carica i CSV
  df_lijiang = pd.read_csv("fucino_satellite.csv")
  df_delingha = pd.read_csv("lario_satellite.csv")

  df_lijiang["Source"] = "Lijiang"
  df_delingha["Source"] = "Delingha"

  # Unisci per calcolare il tempo globale
  df_all = pd.concat([df_lijiang, df_delingha], ignore_index=True)
  print(df_all.columns)
  def to_datetime(row):
      return datetime.strptime(
          row["Date"] + " " + row["Time"],
          "%d %b %Y %H:%M:%S.%f"
      )

  df_all["Datetime"] = df_all.apply(to_datetime, axis=1)

  t0 = df_all["Datetime"].min()
  df_all["SecondsFromStart"] = (df_all["Datetime"] - t0).dt.total_seconds()

  # Ora ri-separa, SENZA altri merge
  df_lijiang = df_all[df_all["Source"] == "Lijiang"].copy()
  df_delingha = df_all[df_all["Source"] == "Delingha"].copy()

  def atmospheric_attenuation_1550(elev_deg, alpha0_db=0.1):
      elev_rad = np.deg2rad(np.abs(elev_deg))
      airmass = 1.0 / np.sin(elev_rad)
      return alpha0_db * airmass

  df_lijiang["AtmosphericLoss_dB"] = df_lijiang["Elevation"].apply(atmospheric_attenuation_1550)
  df_delingha["AtmosphericLoss_dB"] = df_delingha["Elevation"].apply(atmospheric_attenuation_1550)

  Alice_orb_dict = df_lijiang.to_dict(orient="list")
  Bob_orb_dict   = df_delingha.to_dict(orient="list")

  Alice_orb_dict["Range"] = np.array(Alice_orb_dict["Range"]) * 10**3
  Bob_orb_dict["Range"]   = np.array(Bob_orb_dict["Range"])   * 10**3












  ###rewriting atmospheric loss in db per meter ####
  att_factor= 1 #10**(0)
  for k in range(len(Bob_orb_dict["AtmosphericLoss_dB"])):
    Bob_orb_dict["AtmosphericLoss_dB"][k]=(Bob_orb_dict["AtmosphericLoss_dB"][k]/Bob_orb_dict["Range"][k])*att_factor
  for k in range(len(Alice_orb_dict["AtmosphericLoss_dB"])):
    Alice_orb_dict["AtmosphericLoss_dB"][k]=(Alice_orb_dict["AtmosphericLoss_dB"][k]/Alice_orb_dict["Range"][k])*att_factor


  ############# getting parameters from the json file ######################
  runtime =parameters["general_parameters"]["runtime"]*(10)**(12)
  emission_frequency=parameters["emitter_parameters"]["emitter_frequency"]*10**(-12) ###in numbers of emission per picoseconds
  average_pairs=parameters["emitter_parameters"]["av_ph_emitted"]
  detector_efficiency=parameters["Detector_parameters"]["det_eff"]
  dark_counts=parameters["Detector_parameters"]["dark_counts"]
  count_rate_max=parameters["Detector_parameters"]["count_rate_max"]
  time_resolution=parameters["Detector_parameters"]["time_resolution"]
  PBS_eff=parameters["Detector_parameters"]["PBS_efficiency"]
  BS_eff=parameters["Detector_parameters"]["Beam_splitter_efficiency"]
  distance_CA = 1 # dummy value
  distance_CB = 100 # dummy value
  ch_attenuationCA= 0 ## dummy value to impose that the attenuation is  high
  ch_attenuationCB= 0 ## dummy value to impose that the attenuation is  high
  LS=constants.speed_of_light*parameters["channel_parameters"]["light_speed"]*10**(-12) ##light speed expresses in meters per picosecond
  ########### setting up the parameters #################
  emission_period=(int(1/emission_frequency)) #period between two emissions
  sim_time=runtime*10**(-12)
#  print("The script will simulate events up to : ", sim_time, " Seconds")
  tl = Timeline(runtime)
  tl.show_progress = True
#### Defining nodes, Alice and Bob (the receivers) and Charlie (the sender), and quantum channels ################
  Alice = BBM92_receiver(name="Alice", tl=tl,det_eff=detector_efficiency, dark_counts=dark_counts, count_rate_max=count_rate_max, time_resolution=time_resolution, reg=Reg_Alice, PBS_efficiency=PBS_eff, BS_efficiency=BS_eff)
  
  Bob = BBM92_receiver(name="Bob",  tl=tl,det_eff=detector_efficiency, dark_counts=dark_counts ,count_rate_max=count_rate_max ,time_resolution=time_resolution, reg=Reg_Bob,PBS_efficiency=PBS_eff, BS_efficiency=BS_eff)

  qc0 = Qchannel("qc0", tl, distance=distance_CA, polarization_fidelity=1, attenuation=ch_attenuationCA,light_speed=LS)
  qc1 = Qchannel("qc1", tl, distance=distance_CB, polarization_fidelity=1, attenuation=ch_attenuationCB, light_speed=LS)
  Charlie = BBM92_SPDC_source(name="Charlie", tl=tl, receiver_1=qc0, receiver_2=qc1, av_pairs=average_pairs)
  Alice.set_seed(0)
  Bob.set_seed(1)
  Charlie.set_seed(2)
  qc0.connect(Charlie, Alice)
  qc1.connect(Charlie, Bob)

  start=[] #times in seconds at which the emitter should start emitting as the satellite can see both of the receivers
  stop=[] #times in seconds at which the emitter should stop emitting as the satellite cannot see at least one of the receivers
  start.append(max(Alice_orb_dict["SecondsFromStart"][0], Bob_orb_dict["SecondsFromStart"][0])) #instant of coincident visibility
  start_date=[] #times in date format at which the emitter should start emitting as the satellite can see both of the receivers
  stop_date=[] #times in date format at which the emitter should stop emitting as the satellite cannot see at least one of the receivers
  start_date.append(max(Alice_orb_dict["Time"][0], Bob_orb_dict["Time"][0])) #instant of coincident visibility


  #### cheching when the orbits start or stop ####
  max_silence_time=1200 ##maximum silence time of the updating of the orbits in seconds, it needs to be bigger than the orbit time 
  t=start[0]
  i=0
  j=0
  limiter="None"
  Alice_silence_time=0
  Bob_silence_time=0
  on_off_flag=True ##flag which tells us whether the source should emit or not

  if sim_time> max(Alice_orb_dict["SecondsFromStart"]):
    print("Simulation time has been reset to: ", max(Alice_orb_dict["SecondsFromStart"]), " seconds, as it exceeds the simulation time of the orbits" )
    sim_time=max(Alice_orb_dict["SecondsFromStart"])
  if sim_time> max(Bob_orb_dict["SecondsFromStart"]):
    print("Simulation time has been reset to: ", max(Bob_orb_dict["SecondsFromStart"]), " seconds as it exceeds the simulation time of the orbits" )
    sim_time=max(Bob_orb_dict["SecondsFromStart"])

  print("The script will simulate events up to : ", sim_time, " Seconds")



  while t<sim_time:

    if on_off_flag==True:
      if t>=Alice_orb_dict["SecondsFromStart"][i]:
        Alice_silence_time=0
        i=i+1
      else:
        Alice_silence_time=Alice_silence_time+1


      if t>=Bob_orb_dict["SecondsFromStart"][j]:
        Bob_silence_time=0
        j=j+1
      else:
        Bob_silence_time=Bob_silence_time+1

      if Bob_silence_time>max_silence_time :
        stop.append(Bob_orb_dict["SecondsFromStart"][j-1])
        stop_date.append(Bob_orb_dict["Time"][j-1])
        on_off_flag=False
        limiter="Bob"

      if Alice_silence_time>max_silence_time :
        stop.append(Alice_orb_dict["SecondsFromStart"][i-1])
        stop_date.append(Alice_orb_dict["Time"][i-1])
        on_off_flag=False
        limiter="Alice"



    elif on_off_flag==False:
      if t>=Alice_orb_dict["SecondsFromStart"][i] and t>=Bob_orb_dict["SecondsFromStart"][j]:
        start.append(max(Alice_orb_dict["SecondsFromStart"][i], Bob_orb_dict["SecondsFromStart"][j]))
        start_date.append(max(Alice_orb_dict["Time"][i], Bob_orb_dict["Time"][j]))
        on_off_flag=True
     
    else:
      print("Problem with window estimation")        


    t=t+1
  
  if len(start)==len(stop)+1:
    stop.append(sim_time)
    stop_date.append(sim_time)



  #### obtaining coincident windows #######
  print("\n The windows in which the distribution is possible are: \n")
  emission_windows=[]

  if len(start)==len(stop):
    for i in range(len(start)):
      emission_windows.append(stop[i]-start[i])
      print ("Window ", i+1, " from ", start_date[i] , " to ", stop_date[i], " of duration ", f"{emission_windows[i]:.1f}", "seconds \n")
  else:
    print("Error in window estimation " )




  ############## Defining the processes and events #############

  ### scheduling the events where the satellite moves, so when the distance changes for Satellite-Alice channel ####
  for i in range(len(Alice_orb_dict["SecondsFromStart"])):
    CA_channel_update_process= Process(qc0, "update_channel", [Alice_orb_dict["Range"][i] ,Alice_orb_dict["AtmosphericLoss_dB"][i]])
    CA_channel_update_event = Event(Alice_orb_dict["SecondsFromStart"][i]*10**(12), CA_channel_update_process) ## NB seconds are being converted in ps 
    tl.schedule(CA_channel_update_event)

  ### scheduling the events where the satellite moves, so when the distance changes for Satellite-Bob channel ####
  for i in range(len(Bob_orb_dict["SecondsFromStart"])):
    CB_channel_update_process= Process(qc1, "update_channel", [Bob_orb_dict["Range"][i] ,Bob_orb_dict["AtmosphericLoss_dB"][i]])
    CB_channel_update_event = Event(Bob_orb_dict["SecondsFromStart"][i]*10**(12), CB_channel_update_process) ## NB seconds are being converted in ps 
    tl.schedule(CB_channel_update_event)

  

  #### redefining start and stop in picoseconds  ##########
  for i in range(len(start)):
    start[i]=start[i]*10**(12)
  for j in range(len(stop)):
    stop[j]=stop[j]*10**(12)

  ### turning on the emission only within the windows #####
  emission_process = Process(Charlie, "emit_photons", [])

  for i in range(len(start)):
    t=start[i] 
    while  t>=start[i] and t<=stop[i]:
      emission_event= Event(t, emission_process)
      tl.schedule(emission_event)
      t=t+emission_period



  ###### Running the simulation ###############
  tl.init()
  tl.run()


  #######  Getting the dictionary with the result ###########
  results_Alice=Reg_Alice.get_registered_events()
  results_Bob=Reg_Bob.get_registered_events()


  n_Alice_measurements=len(results_Alice["Time"])
  n_Bob_measurements=len(results_Bob["Time"])
  offsets={"Update_time":[], "Distance":[]}
  coincident_counts=0
  Alice_raw_key=[]
  Bob_raw_key=[]

  print(" \nThe number of measurement performed by alice is: ", n_Alice_measurements)
  print("The number of measurement performed by bob is: ", n_Bob_measurements)


  ##### Finding all the time offsets ####

  time_for_off=0
  a=0
  b=0
#  print(len(Alice_orb_dict["SecondsFromStart"]),len(Bob_orb_dict["SecondsFromStart"]) )
  while time_for_off<=sim_time and a<len(Alice_orb_dict["SecondsFromStart"]) and b<len(Bob_orb_dict["SecondsFromStart"]):
#    print(Alice_orb_dict["SecondsFromStart"][a], "   ", Bob_orb_dict["SecondsFromStart"][b], a, "   ", b)
    if Alice_orb_dict["SecondsFromStart"][a]<=time_for_off and Bob_orb_dict["SecondsFromStart"][b]<=time_for_off:
      offsets["Distance"].append((Alice_orb_dict["Range"][a]-Bob_orb_dict["Range"][b])/LS)
#      print( Alice_orb_dict["Range"][a], "  ", Bob_orb_dict["Range"][b])
      offsets["Update_time"].append(Alice_orb_dict["SecondsFromStart"][a]*10**(12)+(Alice_orb_dict["Range"][a]/LS)) ### calculating the offsets in picosecond
      a=a+1
      b=b+1
#      print("AAAAAA both orbits changed at the same time ")
    if Alice_orb_dict["SecondsFromStart"][a]<=time_for_off:
      offsets["Distance"].append((Alice_orb_dict["Range"][a]-Bob_orb_dict["Range"][b])/LS)
      offsets["Update_time"].append(Alice_orb_dict["SecondsFromStart"][a]*10**(12)+(Alice_orb_dict["Range"][a]/LS)) ###converting only the offsets in picoseconds
      a=a+1
#      print("AAAAAA Alice orbits changed first ")
    if Bob_orb_dict["SecondsFromStart"][b]<=time_for_off:
      offsets["Distance"].append((Alice_orb_dict["Range"][a]-Bob_orb_dict["Range"][b])/LS)
      offsets["Update_time"].append(Bob_orb_dict["SecondsFromStart"][b]*10**(12)+(Alice_orb_dict["Range"][a]/LS)) ###converting only the offsets in picoseconds
      b=b+1
#      print("AAAAAA Bob orbits changed first ")
    time_for_off=time_for_off+0.01







#  for p in range(30):
#    print(results_Alice["Time"][p]*10**(-12)- results_Bob["Time"][p]*10**(-12), results_Alice["Time"][p]*10**(-12))


#  for l in range(len(offsets["Update_time"])):
#    print(offsets["Update_time"][l]*10**(-12), "   ", offsets["Distance"][l])

  ### Defining parameters and arrays for the plots and post processing ####
  research_window=10**(20)
  average_QBERS=np.zeros(int(len(start)))
  average_KTR=np.zeros(int(len(start)))
  upd_time_array = [[] for _ in range(len(start))]
  QBER_array=[[] for _ in range(len(start))]
  KR_array=[[] for _ in range(len(start))]
  init_Alice=0
  init_Bob=0

#  print(offsets)


  ##### post processing each orbit ####
  for w in range(len(start)):
    print("\n Post processing orbit ", w+1)
    #resetting every parameter for this orbit
    coincident_counts=0
    Alice_raw_key=[]
    Bob_raw_key=[]
    Alice_RK_var=[]
    Bob_RK_var=[]
    start_Alice_pp=0
    start_Bob_pp=0
    stop_Alice_pp=0
    stop_Bob_pp=0
    k=0
    offset=0
    
    ### finding the first measurement of alice which is inside the window ## 
    for a in range(init_Alice,len(results_Alice["Time"]),1):
      if results_Alice["Time"][a]>=start[w]:
        start_Alice_pp=a
        init_Alice=a
        break

    ### finding the last measurement of alice which is inside the window ## 
    for a in range(init_Alice, len(results_Alice["Time"]),1):
      if results_Alice["Time"][a]>=stop[w]:
        stop_Alice_pp=a
        init_Alice=a
        break
      else:
        stop_Alice_pp=len(results_Alice["Time"])


    ### finding the first measurement of bob which is inside the window ## 
    for b in range(init_Bob,len(results_Bob["Time"]),1):
      if results_Bob["Time"][b]>=start[w]:
        start_Bob_pp=b
        init_Bob=b
        break


    ### finding the last measurement of Bob which is inside the window ## 
    for b in range(init_Bob,len(results_Bob["Time"]),1):
      if results_Bob["Time"][b]>=stop[w]:
        stop_Bob_pp=b
        init_Bob=b
        break
      else:
        stop_Bob_pp=len(results_Bob["Time"])


    #### finding the offset at start time ###
    for h in range(len(offsets["Distance"])-1):
      if offsets["Update_time"][h]<=results_Alice["Time"][start_Alice_pp] and  offsets["Update_time"][h+1]>=results_Alice["Time"][start_Alice_pp]: 
#        print("offset updated at: ", offsets["Update_time"][h] ) 
        offset=offsets["Distance"][h]
        k=h+1
        break
    if offset==0:
      print("Problem with individuation of the offset")



    #### recovering the coincident measurements #####
    for i in range (start_Alice_pp, stop_Alice_pp, 1):

      #### updating offset and calculating QBER and key rate for that offset###
      if k<len(offsets["Update_time"]):
        if offsets["Update_time"][k]<=results_Alice["Time"][i]:
#          print(k, "   ", len(offsets["Update_time"]) )
          offset=offsets["Distance"][k]
          upd_time_array[w].append(offsets["Update_time"][k])
          if (offsets["Update_time"][k]-offsets["Update_time"][k-1])==0:
            print("Problem with kr at: ", offsets["Update_time"][k], "   " ,offsets["Update_time"][k-1])
          else:
            KR_array[w].append(len(Alice_RK_var)/((offsets["Update_time"][k]-offsets["Update_time"][k-1])*10**(-12)))
          k=k+1
          ccfr=0 #correct counts for QBER array estimation 
          cdfr=0  #correct counts for QBER array estimation

          
          for d in range(len(Alice_RK_var))  :
            if Alice_RK_var[d]==Bob_RK_var[d]:
              ccfr=ccfr+1
            else:
              cdfr=cdfr+1

          if len(Alice_RK_var)==0:
            QBER_array[w].append("None")
          else:
            QBER_array[w].append((cdfr/len(Alice_RK_var)))

          Alice_RK_var=[]
          Bob_RK_var=[]



      for j in range(start_Bob_pp,stop_Bob_pp,1):

        if np.abs(results_Alice["Time"][i]-results_Bob["Time"][j]-offset)<=2*time_resolution and results_Alice["Basis"][i]==results_Bob["Basis"][j]:
#          print(results_Alice["Time"][i], "   ", offsets["Update_time"][k])
          coincident_counts=coincident_counts+1
          Alice_raw_key.append(results_Alice["Det_result"][i])
          Bob_raw_key.append(results_Bob["Det_result"][j])
          Alice_RK_var.append(results_Alice["Det_result"][i])
          Bob_RK_var.append(results_Bob["Det_result"][j])

          start_Bob_pp=j
          break
        elif  np.abs(results_Alice["Time"][i]-results_Bob["Time"][j]-offset)> research_window: # and results_Alice["Time"][i]-results_Bob["Time"][j]-offset>0:
          break
        
      
#    print(QBER_array[w])
#    print(KR_array[w])
#    print(upd_time_array[w])
    print("Total length of the raw key for the orbit", w+1, " is: ",coincident_counts )
    key_rate=coincident_counts/((stop[w]-start[w])*10**(-12))
    average_KTR[w]=key_rate
    print("The average key transmission rate throughout the orbit is: ", f"{key_rate:.3g}" , " Qubit/sec")

    correct_counts=0
    different_counts=0
    for k in range(coincident_counts)  :
      if Bob_raw_key[k]==Alice_raw_key[k]:
        correct_counts=correct_counts+1
      else:
        different_counts=different_counts+1

    if coincident_counts==0:
      QBER="None"
#      average_QBERS[w]=QBER
      print("No coincident counts have been detected for the orbit ",w+1, ", therefore QBER cannot be estimated")
    else:
      QBER=(different_counts/coincident_counts)
      average_QBERS[w]=QBER
      print(f"The average QBER of the orbit is: {QBER:.3f}")
    
  ### rescaling time from the start of the orbit and reprocessing ###
  for w in range(len(start)):
    resc_fact=min(upd_time_array[w])
    for l in range(len(upd_time_array[w])):
      upd_time_array[w][l]=(upd_time_array[w][l]- resc_fact)*10**(-12)
      if (QBER_array[w][l]=="None"):
        QBER_array[w][l]=0.5


  ### plotting the change in QBER for each orbit #### 
  plt.figure(1)
  plt.title("QBER along each orbit")
  for w in range(len(start)):
    QBER_array[w]=QBER_array[w][1:]
    upd_time_array[w]=upd_time_array[w][1:]
    n_win=w+1
    plt.plot(upd_time_array[w],QBER_array[w],marker='o', markersize=8, label=f"Window {n_win}")
    
  plt.ylabel("QBER")
  plt.xlabel("Time since start of the orbit (s)")
  plt.legend()
  plt.savefig("Grafico_QBER.png", dpi=300)
  plt.show()
  plt.close()        


  ### plotting the change in Key_rate for each orbit ####
  plt.figure(2)
  plt.title("Key transmission rate along each orbit")
  for w in range(len(start)):
    KR_array[w]=KR_array[w][1:]  
    n_win=w+1
    plt.plot(upd_time_array[w],KR_array[w],marker='o', markersize=8,label=f"Window {n_win}")
    
  plt.ylabel("Key transmission rate (bits/s)")
  plt.xlabel("Time since start of the orbit (s)")
  plt.legend()
  plt.savefig("Grafico_Key_Rate.png", dpi=300)
  plt.show()
  plt.close()        

    

if __name__ == "__main__":
  main()
