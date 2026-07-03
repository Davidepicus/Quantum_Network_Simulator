# Hello

Questa è la repository di Github per quantum network simulator. Il pacchetto usato è SeQUeNCe.

Il link alla repository di Github che contiene il codice del pacchetto SeQUeNCe è: https://github.com/sequence-toolbox/SeQUeNCe

Il link alla documentazione del pacchetto SeQUeNCe è: https://sequence-rtd-tutorial.readthedocs.io/stable/index.html

# Files

Per ora ci sono tre files nella repository. "BB84_test.py" è il test di una simulazione di un protocollo BB84. Questo file è una copia esatta del codice che si trova nella repository di Sequence. 
Il secondo file, "BBM92_simulation.py" simula un protocollo BBM92 con dei nodi di ricezione ed invio creati appositamente. 

Vi sono due classi che costituiscono i due tipi di nodi necessari per il protocollo BBM92.
"BBM92_SPDC_source" è la classe che contiene l'emettitore di fotoni, in questo caso un cristallo SPDC eccitato con una certa frequenza.

"BBM92_receiver" è invece la classe che simula un ricevitore BBM92. Il ricevitore è costituito fa un beam splitter 50/50 2X1, due polarizing beam splitter che misurano la polarizzazione dei fotoni nelle diverse (base {H,V} e base {+.-}), e 4 detector.


# Stato della Simulazione

Per ora le due classi sono funzionanti e connesse. I fotoni vengono emessi ed arrivano sui detector. Manca ancora la parte di post processing classica, ma non dovrebbe essere troppo complicata
