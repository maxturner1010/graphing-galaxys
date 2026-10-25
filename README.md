# graphing-galaxy-s
this projects goal is to automatically fit good mathematical functions with an appropriate amount of paramaters to data.

Originally, this project was simply a way for me to learn the basics of scientific data analysis through learning some uses of NumPy, matplotlib, htm5, and sciPY.
I read astrophysics data from a htm5 file then sanitised the data. I then used curve fit (from sciPY) with some polynomial functions with 3,4,5 parameters and then compared them to the original data using chi2 to see which model was best before visualising it by plotting the models using matplotlib.

Since then I have added extra functionality so that the code can find a good type of function with an appropriate amount of parameters to fit multiple sets of data found in the htm5 file. I have done this by finding the reduced chi2, AIC, and BIC for multiple function types (right now polynomial, logarithmic, power) for multiple amount of parameters. I store the best models for chi2,AIC and BIC and then I can chose between them and if they pick one result I can be fairly confident it is the best model that is within my models search parameters. However, I have found that the BIC and AIC normally agree on a model while the reduced chi2 normaly fits a function that is too close to the data possibly not reflecting an acurate modle.

In plotting the models against the data they are supposed to represent using matplotlib I have found the modules look like they represent the data it is trying to fit to showing that the code does work.

Currently I am plotting astronomy data of magnitude against galaxy's per area.

Future additions:
- I would like to test how accurately my code can predict the models that data can represent by getting points from a function, adding noise and then seeing if the original model is correctly predicted. I could test for different amounts of noise and different amounts of data points to test how effective AIC, BIC, reduce chi2 are under different types of data (noise level, function type, number of data points).
- add more types of models
- possibly find more ways to access modles
- turn it into a class so different data sets can be assessed more easily
- allow for other types of data files other than htm5
