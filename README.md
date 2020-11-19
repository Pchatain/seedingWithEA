# seedingWithEA
All the code for my CS229 Final Project are located in GenericDetector and LDMX_Detector. GenericDetector contains the code I used for the 2 methods on the generic detector in ACTS, for both the muon and the ttbar samples. LDMX_Detector contains the code I used to implement both methods on the LDMX_Simulated data. The bulk of the code I wrote is in the files ParameterSweep.py and EvolutionaryStrategies.py. In addition, the LDMX_Detector contains the file I modified in order to calculate the metrics for the seeding algorithm.

testingEA.py contains an evolutionary algorithm that optimizes the seeding algorithm configuration from the acts-project.

Utilizes my local copy of acts-project to run with a modified print output and modified seeding algorithm.

testingScript.py contains a script to plot the efficiency and number of seeds generated versus a configurable parameter. Utilizes multiprocessing.
