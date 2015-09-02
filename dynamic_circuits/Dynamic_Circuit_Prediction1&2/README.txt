Dynamic_Circuit_Prediction
takes into account mRNA concentration, txn tln rate separately


Dynamic_Circuit_Prediction2
combines txn/tln into one term when determining REU (not protein concentration)



How to run a test circuit that 
open Final.py in Spyder (python IDE).
press 'Play' to load the file.
in the prompt, type makeGraphFromNetlist(exampleNetlist, placeToSaveExample)
... the netlist is converted to JSON and saved as a file, that file is read by General.py.


It's also possible to press 'Play' on OptimalCircuit.py, and run
wrapperForNetlist("JsonFiles/SR_Latches/SRLatch.json", "JsonFiles/SR_Latches/SRLatch_Graph.json", makeBarGraph=True,makeOtherGraphs=True,useDefaultInput=False)
(add more details here...)
(could also import Final.py and press play to access the variable names)


