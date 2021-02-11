# hyperpaths
## Prerequisites

The code for constructing directed hypergraph objects (required for the hyperpath heuristic, cut finding procedure, and cutting planes algorithm) is available at https://github.com/annaritz/pathway-connectivity.
Once you have converted the OWL file from Reactome or Pathway Commons, you should have a -hypernodes.txt file and a -hyperedges.txt file, which will be read in by run.py to construct the directed hyperpath object.

To run the hyperpath heuristic and cut-finding procedure, only python is required, as well as halp, which can be installed via pip.

For the cutting-planes algorithm, the LP solver CPLEX is required, which can be installed through https://www.ibm.com/support/knowledgecenter/SSSA5P_12.7.1/ilog.odms.cplex.help/CPLEX/GettingStarted/topics/set_up/Python_setup.html

## Installation

To install, you must create a folder in the directory called parsed/ and then set ROOTDIR in run.py to be the current directory. It is in this parsed/ directory where your hypernodes.txt and hyperedges.txt files should go.

## Algorithms

The hyperpath heuristic can be run using the following command:

```
 python run.py --heuristic --source <source> --target <target> --name <prefix of hypernodes.txt and hyperedges.txt file>
```

The cut finding and cutting planes procedure can be run alone (assuming the hyperpath heuristic was already run to find the taildistance list) using the following command:

```
 python run.py --cuttingplanes --source <source> --target <target> --name <prefix of hypernodes.txt and hyperedges.txt file> --taildistancelist <pkl file of taildistances from the hyperpath heuristic>
```
The heuristic, cut finding and cutting planes procedure can be run together using the following command:

```
 python run.py --heuristic --cuttingplanes --source <source> --target <target> --name <prefix of hypernodes.txt and hyperedges.txt file> 
```

Currently under construction...


