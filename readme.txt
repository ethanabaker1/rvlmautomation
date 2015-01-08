Light Study Data Automation Tool
Ethan Baker eab109(at)pitt.edu
Yates Lab, University of Pittsburgh
####################################
Installation Notes:
This software requires MATLAB R2014b (or later, some earlier versions may also work, but MATLAB 2010 does not support this code).
This software relies on two .m files, risingEdge.m and risingEdge2.m (Latest version does not require risingEdge2.m)
All code files should be placed in the default MATLAB directory.

###################################
Running the Data Automation Tool
Open MATLAB. Open the rvlmautomationcode2.m file. Select the green run triangle on the top of the MATLAB window. 
A dialogue box will appear asking you to choose a directory. Select the directory containing SPIKE 2 files
exported as .mat files. Before exporting, it is essential to generate the HR channel. 
Output tables will be generated within MATLAB.

###################################
How the code works:
1. Creates matricies of Spike2 Data v. Time (generated from start time, end time, interval)
2. Applies a Savitzky-Golay filter to the light channel (Essentially, this eliminates noise and places a 7th order polynomial function over the data, creating local minima at the light onset and light offset)
3. Inverts all of the S-G filtered data to make local minima local maxima.
4. Gets times of all ascending crosses of a threshold (currently hard-coded, but this may present an issue in some cases)
5. Depending on number of runs, calculate light on and off times for all crosses.
6. Using a second set threshold, find all ascending crosses in the table channel. 
7. Calculates all tilt onset and offset times. 
8. Depending on what regions of interest exist in the file, calcualtes mean blood flow and HR.
9. Generates an output array
10. Converts output array to output tables. 
(11. Sends output table to Excel file) not implemented yet

####################################
Known Issues:
Due to variations in ambient light during existing runs, the light threshold is not consistent.
This may result in some instances of the program finding too many light onset times or none at all.
This can be resolved by modifying the threshold value in risingEdge2.m. I am working
on a less involved solution. UPDATE 11/21/14: Should be fixed, requires testing. risingEdge2 is no 
longer required. 

Code is very long (2700 lines) - there may be a way to shorten it for easier debugging w/loops, but this presents
a challenge given variable tilt and run counts.Probably won't provide a speed advantage (program already fast), but
would improve readability. 

1/8/15 - Modifications made to reflect new channel numbers. Analysis code now accessible as function analyzeRVLM(). 
All dependencies on risingEdge2.m removed. Threshold now calculated based on preset increase over pre-tilt light level, rather than preset value. 



