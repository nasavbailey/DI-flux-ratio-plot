# CGI-flux-ratio-plot
This contains the python code to plot the flux ratio vs separation plot. It also contains the text files with the contrast curves and planet contrast values.

The most recent updated plot is flux_ratio_plot.pdf (no ELT curve) and flux_ratio_plot_ELT.pdf (including the ELT curve). 

Datafiles should be formatted in the following way to allow for auto import by astropy.io.ascii & auto generation of caption.txt:
#col_name1  col_name2
#Short caption: _Write a short explanation, suitable for a figure caption._
#Full explanation: _Write any/all details needed to understand this dataset._
data11  data12
data21	data22
...

