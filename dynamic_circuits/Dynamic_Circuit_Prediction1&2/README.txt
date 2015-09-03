Dynamic_Circuit_Prediction
Takes into account mRNA concentration, transcription and translation rates separately. The optimizing algorithm for Dynamic_Circuit_Prediction may be incomplete compared to Dynamic_Circuit_Prediction2


Dynamic_Circuit_Prediction2
Combines transcription and translation into one term when determining REU (not protein concentration).



How to run a test circuit that:
open Final.py in Spyder (python IDE).
press 'Play' to load the file.
in the prompt, type makeGraphFromNetlist(exampleNetlist,genesToUse=genesList)
... the netlist is converted to JSON and saved as a file, that file is read by General.py.


It's also possible to press 'Play' on OptimalCircuit.py, and run
wrapperForNetlist("JsonFiles/SR_Latches/SRLatch.json", "JsonFiles/SR_Latches/SRLatch_Graph.json", makeBarGraph=True,makeOtherGraphs=True,useDefaultInput=False)
(add more details here...)
(could also import Final.py and press play to access the variable names)

You can do various things with this program:

makeGraphFromNetlist and makeGraphFromCircuitString
You can run a single circuit using either default values or values from the gene library if you specify which genes to use. If you try to specify which genes to use, you must give the correct number of genes or an error will be thrown.
You can give the inputs a unique wave form or use the default by changing the 'useDefaultWaveForm' function input.
example: makeGraphFromNetlist(exampleNetlist,genesToUse=genesList)

optimizeNetlistWithLibraries, optimizeCircuitStringWithLibraries, optimizeNetlistWithLibrariesTimed, and optimizeCircuitStringWithLibrariesTimed
The last two do the same thing as the first two except you can set a time limit with the last two. Given a netlist and a library of genes they will use a Hill Climbing algorithm to determine the best arrangement of genes for the netlist. You can adjust the minimum score allowed for a gate, the number of trajectories, and the max time allowed for the last two.
example: optimizeNetlistWithLibraries(truthValueExampleFileLoc2,smallestScoreAllowed=3)

Look at examplesForUse() for more examples.

If odeint is creating problems go to General.py and comment out
soln_f = odeint(f, initc, t, args=(input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,numGates),hmax=350)
and uncomment out
soln_f = DifferentialSolver.differentialSolver(f, initc, t, args=(input_and_logic_gate_dictionaries,input_and_logic_gate_names,logic_gate_names,numGates))
