README for "PMTAnalysis" ROOT CLASS
By: Eli Brunner-Huber for the ANNIE Collaboration
Last Edited: 07/18/24

PMTAnalysis is a ROOT class to be used for the analysis of tank PMT performance during laser runs in the ANNIE collaboration.
In order to use PMTAnalysis, one must have ROOT installed on their machine, as well as have the root file with raw laser data in the directory one calls PMTAnalysis from.

A quick list of prerequisites before using PMTAnalysis:

1.) ROOT must be installed
2.) A raw data file must be in the current directory
3.) User must have read/write permissions in the current directory
4.) Must have GhostScript installed (?)

Once these prerequisites are met, you are ready to run PMTAnalysis. To do so, enter the following command into the command prompt (NOT in root):

username@machine ~ % root PMTAnalysis.C 

This will simultaneously open ROOT and initialize PMTAnalysis. Now you have to create an instance:

root [1] PMTAnalysis instance;

This will generate what feels like an infinity of output lines in the command prompt. I'm working on suppressing these, but for the meantime, ignore them.

Now you have created an instance of PMTAnalysis, and are ready to use it. The only function from PMTAnalysis you need to call is, you guessed it, "Analyze". The parameters
for "Analyze" are:

1.) Input File
2.) Start time (of the laser run)
3.) End time (of the laser run)
4.) Cutoff voltage (to avoid high voltage noise)

Once ready, run:

root [2] instance.Analyze("inputFileName", startTime, endTime, cutoffVoltage)

This will also generate a lot of output lines. Again, ignore them for now.

Your analyzed ROOT tree is ready! In addition, a pdf has been created in that same directory that summarizes the most important information for easy sharing. In order to view,
restart ROOT and open a TBrowser in the desired directory. 

