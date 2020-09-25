Expander7.0 batch mode:

executable6.5Batch.jar allows using few of Expander functionality in Batch mode.
In order to use this file, it should be placed inside a regular Expander6.5 directory.
The following functionalities are currently supported:
1) Matisse (option 3) 
2) Degas (option 4)
3) Tango (option 5)
Additional options will be added in future versions

Required parameters:
1) Matisse (option 3):
	Matisse batch mode requires the following parameters
	Gene expression data file path
	Network file path
	Organism name (as it appears in the file organisms.txt)
    Beta value (between 0 and 1 default recommanded is 0.95)
	Minimal module size
	Maximal module size;
	Path of results file;

An example comman line:
java -jar -Xms2500m -Xmx3200m executable6.5Batch.jar 3 C:/Expander6.5/sample_input_files/expressionData1.txt C:/Expander6.5/organisms/s.cerevisiae/interactions/Expander.Saccharomyces_cerevisiae.Biogrid.sif s.cerevisiae 0.95 4 100 C:/Temp/matisseBatchRes.txt

* GE data file should contain a unique gene ID column and a gene Symbol column
* Gene IDs should match the type supported by Expander for the investigated organism
 (for further details refer to the Expander user manual in the file help.html)
 
2)Degas (option 4)
Degas batch mode requires the following parameters:
Gene expression data file path
Network file path
Organism name (as it appears in the file organisms.txt
Names of all case conditions separated by ;
Names of all control conditions separated by ;
Optimization algorithm (0 for CONNECTING_GREEDY, 1 for CUSP, 2 for CUSP_STAR, 3 for EXTENDING_GREEDY, 4 for EXTENDING_GREEDY_TWICE or 5 for VANILLA_GREEDY (recommanded default=1)
Dysregulation direction (-1 for DOWN, 1 for UP, 0 for BOTH - recommanded default=0
Dysregulation significanse threshod (between 0 and 1, recommanded default=0.05)
Dysregulation ratio(recommanded default=1.3)
Max modules (recommanded default =1)
Min tested K (recommanded default =10)
Max tested K (recommanded default =50)
Tested K step (recommanded default =10)
l parameter (recommanded default =1)
Path of results file (including file name)

An example command line:
java -jar -Xms2500m -Xmx3200m executable6.5Batch.jar 4 C:/Expander6.5/sample_input_files/expressionData1.txt C:/Expander6.5/organisms/s.cerevisiae/interactions/Expander.Saccharomyces_cerevisiae.Biogrid.sif 
Cond1;Cond2;Cond3;Cond4 Cond25;Cond26;Cond27;Cond28 1 0 0.05 1.3 1 C:/Temp/DegasBatchRes.txt > degasLog.txt

* GE data file should contain a unique gene ID column and a gene Symbol column
* Gene IDs should match the type supported by Expander for the investigated organism
 (for further details refer to the Expander user manual in the file help.html)
 
 
3) Tango (option 5)
Tango batch mode operates on predefined gene clusters, and requires the following parameters:
Organism name (as it appears in the file organisms.txt)
Clustering file path - for format details see help.html file
Background file path
number of iterations
Maximal GO class size (level of generality in the tree) - 3000 is the default in the GUI
Threshold p-value
Result file full path

An example command line:
java -jar -Xms2500m -Xmx3200m executable6.5Batch.jar 3 s.cerevisiae C:/Expander6.5/sample_input_files/expressionData1Clustering.sol C:/Expander6.5/sample_input_files/scvBg.txt  1000 3000 0.05 C:\Temp\TangoRes.txt
				
* Both gene clusters file and background set file should contain gene IDs of the type supported by Expander for the investigated organism
 (for further details refer to the Expander user manual in the file help.html)

