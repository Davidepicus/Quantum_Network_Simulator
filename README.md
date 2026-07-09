# Hello

Questa è la repository di Github per quantum network simulator. Il pacchetto usato è SeQUeNCe.

Il link alla repository di Github che contiene il codice del pacchetto SeQUeNCe è: https://github.com/sequence-toolbox/SeQUeNCe

Il link alla documentazione del pacchetto SeQUeNCe è: https://sequence-rtd-tutorial.readthedocs.io/stable/index.html



Per ora ci sono due  files nella repository. "BB84_test.py" è il test di una simulazione di un protocollo BB84. Questo file è una copia esatta del codice che si trova nella repository di Sequence. 

# BBM92 protocol simulation


Il file "BBM92_simulation.py" simula un protocollo BBM92 con dei nodi di ricezione ed invio creati appositamente. 

Vi sono due classi che costituiscono i due tipi di nodi necessari per il protocollo BBM92.
"BBM92_SPDC_source" è la classe che contiene l'emettitore di fotoni, in questo caso un cristallo SPDC.

"BBM92_receiver" è invece la classe che simula un ricevitore BBM92. Il ricevitore è costituito da un beam splitter 50/50 1X2, due polarizing beam splitter che misurano la polarizzazione dei fotoni nelle diverse basi (base {H,V} e base {+.-}), e 4 detector.

Vi è inoltre una classe chiamata "det_counter" che serve per registrare gli arrivi e le misure dei fotoni, mentre la classe "Register" permette di registrare in un dizionario I risultati delle misure di ognuno dei due nodi, il tempo di arrivo dei fotoni e la base utilizzata per la misura. Vi sono due classi beam splitter che servono per modificare alcuni bug e incongruenze del codice Sequence originario. Il tutto va ancora commentato, ottimizzato e pulito


# Stato della Simulazione

Alice, Charlie e Bob sono connessi. La simulazione funziona ma va ancora otimizzata. Vanno inseriti i conteggi oscuri, il rumore dovuto ai canali, e va settata bene la frequenza di emissione. Va inoltre introdotto un file di config che possa rendere più semplice l'inserimento dei parametri vari. Questo è tuttavia il primo prototipo funzionante.
