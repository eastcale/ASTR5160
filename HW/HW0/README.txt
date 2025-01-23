This code generates 10 random floats in the range 0-10 and generates corresponding line data when passed values for m and b in the expression y=m*x+b. It then randomly offsets each of the resulting y values by some magnitude determined by a Gaussian distribution with standard deviation 0.5 centered at 0. A line of best fit is found from this new data. A plot will be generated that contains the original line (with no offset, utilizing the user inputs for m and b values), the offset data points, and the line of best fit. A png titled 'line-plot.png' will be saved in your present working directory.

To run HW0.py in a python shell:

python HW0.py

You will be prompted to input a slope, as well as a y-intercept. These must be numbers (floats or integers are acceptable).



