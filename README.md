# Benvenuti

Questa è la repository di Github per quantum network simulator. Il pacchetto usato è SeQUeNCe.

Il link alla repository di Github che contiene il codice del pacchetto SeQUeNCe è: https://github.com/sequence-toolbox/SeQUeNCe

Il link alla documentazione del pacchetto SeQUeNCe è: https://sequence-rtd-tutorial.readthedocs.io/stable/index.html



Per ora ci sono tre  files  importanti nella repository.
"BB84_simulation.py" è il test di una simulazione di un protocollo BB84. Devo ancora iniziare a lavorare a questa simulazione.
Data_organyzer.py è un codice che serve per aprire i dati delle orbite e salvarli nel formato csv per una gestione più facile 

# BBM92 protocol simulation


Il file "BBM92_simulation.py" simula un protocollo BBM92 con dei nodi di ricezione ed invio creati appositamente. 

Vi sono due classi che costituiscono i due tipi di nodi necessari per il protocollo BBM92.
"BBM92_SPDC_source" è la classe che contiene l'emettitore di fotoni, in questo caso un cristallo SPDC.

"BBM92_receiver" è invece la classe che simula un ricevitore BBM92. Il ricevitore è costituito da un beam splitter 50/50 1X2, due polarizing beam splitter che misurano la polarizzazione dei fotoni nelle diverse basi (base {H,V} e base {+.-}), e 4 detector.

Vi è inoltre una classe chiamata "det_counter" che serve per registrare gli arrivi e le misure dei fotoni, mentre la classe "Register" permette di registrare in un dizionario i risultati delle misure di ognuno dei due nodi, il tempo di arrivo dei fotoni e la base utilizzata per la misura. Vi sono due classi beam splitter che servono per modificare alcuni bug e incongruenze del codice Sequence originario. Una terza classe chiamata qChannel, è presente al fine di ridefinire il modo in cui si connettono i canali con i nodi.

# Stato della Simulazione

Alice, Charlie e Bob sono connessi. La simulazione funziona ma va ancora otimizzata, in particolare il post processing. 
I parametri che possono essere modificati sono:

I conteggi oscuri, il dead time ed il time jitter dei detector.

L'efficienza dei beam splitter 50/50 1X2, e l'efficienza dei polarizing beam splitter.

La frequenza di stimolazione del cristallo SPDC e il valore medio di coppie di fotoni emesse per stimolazione.

L'attenuazione dei canali quantistici, e la velocità della luce in essi.

Il tempo di simulazione. 

Tutti questi parametri possono essere modificati dal file di configurazione "config.json". Questo prototipo è funzionante, anche se i paramteri vanno introdotti con cura per evitare un tempo di simulazione e post processing troppo lungo. 

Manca l'inserimento degli effetti dovuti alla depolarizzazione dei fotoni all'interno dei canali quantistici, anche se per considerare questi effetti si potrebbe cambiare la fedeltà dei beam splitter. Per ora l'attenuazione dei canali è nulla in quanto il tempo di simulazione necessario sarebbe troppo lungo. 

Il satellite si muove, e i fotoni vengono emessi e ricevuti correttamente







