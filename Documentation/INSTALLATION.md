# Installation Instructions for Fiber-Orientation-Tools

https://github.com/charlestucker3/Fiber-Orientation-Tools

If you downloaded a ZIP file and unzipped it, you should have a folder named `Fiber-Orientation-Tools-main`.  Move this folder to a convenient long-term location for the documentation and examples files.  

The `Fiber-Orientation-Tools` sub-folder contains all the the Matlab functions.  You want to make these functions accessible from Matlab no matter which directory you are working in.  To do this:

1.  Start Matlab and type the command `userpath`.  This will show the path to a user directory that Matlab adds to its search path every time it starts up.  On my Mac the default `userpath` is `/Users/<username>/Documents/MATLAB`.  This path will look different on Windows or Linux machines.  We will call this your `userpath` directory.  
 
2.  Navigate to that folder, and look for a file named `startup.m`.  This is a set of Matlab commands that are executed every time you start Matlab.  If this file exists, open it with Matlab.  If it does not exist,  open a new script in Matlab, and save it in your `userpath` directory as `startup.m`.  

3. Edit `startup.m`, adding a line like `addpath /Users/<username>/Documents/MATLAB/Fiber-Orientation-Tools`.  Use the path to your `userpath` directory, followed by `Fiber-Orientation-Tools`.  

4.  Move the `Fiber-Orientation-Tools` directory (the one with the Matlab functions) to your `userpath` directory.  The `Documentation` folder can stay where it is.

5. Quit Matlab, then re-start it.  This should place the `Fiber-Orientation-Tools` directory in your Matlab path.  To confirm that, type `help closeA4`.  This should bring up the help text for the `closeA4` function.  If it does, your installation is complete.  If not, type `path` to see your full Matlab path, and confirm that the directory `Fiber-Orientation-Tools` is in the path.  Also check to make sure that the `Fiber-Orientation-Tools` folder is where it should be.  

6. If you plan to use the `F2A `and `A2F` functions, go to https://www.mathworks.com/matlabcentral/fileexchange/65915-elfun18 and download the `elfun19` package by Milan Batista.  Install this in the same way you installed Fiber-Orientation-Tools.

See the `FiberOrientationToolDoc.pdf` for information on the Matlab functions in the toolkit.  That PDF file also lists the Matlab live scripts in the `Documentation` folder that show how the various tools are used.  