
Instructions for running the code:

Edit the code to select the grid size and to change the initial velocity field. You should, however, be able to run the two example files without any editing.

Compile the code with

f77 ftc2db-2.f

where f77 is your fortran compiler (this could be g77, fort77, or something else). The compilation generates an executable a.out.

Run the program by

a.out<example1

This will generate a directory where your output files will be put. The program generates output to the screen. If you want to redirect the output to a file, use:

a.out<example1>outputfile

To plot the files, start MATLAB in the same directory as you a.out file is and open frd_movie. Edit the top of the program to make sure the name of the files is correct, select the fields that you want to plot, and the number of output files.
