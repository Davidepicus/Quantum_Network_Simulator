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
 

  def connect(self, sender_node=None, receiver_node=None ):
    print(sender_node.name, "Has been connected with ", receiver_node.name)
    self.sender_node=sender_node
    self.receiver_node=receiver_node
    self.set_ends(sender_node, receiver_node.name)


################### Defining class for a BB84 sender ####################

class BB84_laser_source(nd.Node):

  def __init__(self,name="unnamed", freq=1*10**(8), tl = Timeline(2e10), receiver_1=None, av_photons=0.1):
    super().__init__(name, tl)
    self.frequency=freq
    self.name=name
    laser_name= name +" laser_source"
    polarization ={"name": "polarization","bases": [((complex(1), complex(0)), (complex(0), complex(1))),((complex(np.sqrt(1 / 2)), complex(np.sqrt(1 / 2))), (complex(-np.sqrt(1 / 2)), complex(np.sqrt(1 / 2))))]}


    self.laser_source=ls.LightSource(name=laser_name, timeline=tl, wavelength=1550, mean_photon_num=av_photons)


    ###### attaching the SPDC_source to the node ########
    self.add_component(self.laser_source)

    ######## attaching the source to the receiver ############
    self.receiver_1=receiver_1
    self.laser_source.add_receiver(self.receiver_1)



    ######### Informing user about creation of the BB84_laser_emitter  ##########
    print("BB84_laser_source named: ",  self.name ," has been set")


    ######### function that calls the SPDC source and tells it to excite the SPDC crystal ###################
  def emit_photon(self):
        x = random.randint(0, 1) 
        y = random.randint(0, 1) 


        c = 1 / mt.sqrt(2)
        phi_plus = [c + 0j, c + 0j]
        state_list=phi_plus
        self.laser_source.emit( [phi_plus])
        




def main():


#### Parameters settings ############
  average_photons_emitted=0.1
  runtime=2*10**(10)
  distance_AB= 2*10**(3)
  ch_attenuationAB=0.002
  LS=1

##### Defining timeline ################
  tl = Timeline(runtime)
  tl.show_progress = True

#### Defining nodes, Alice (the sender), and quantum channels ################
  qc0 = Qchannel("qc0", tl, distance=distance_AB, polarization_fidelity=1, attenuation=ch_attenuationAB,light_speed=LS)
  Alice = BB84_laser_source(name="Alice", tl=tl, receiver_1=qc0, av_photons=average_photons_emitted)
  



if __name__ == "__main__":
  main()
