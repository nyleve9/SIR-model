# SIR-model
for spread of disease

We analyzed data from the first 40 days of an epidemic that is spreading through Kinamali island (population 6153) and constructed a biological model based on the data to ultimately predict how long the epidemic will run and how many people will ultimately die.  This emerging infectious disease is characterized by unresponsiveness to efforts to quarantine victims at an early stage, ineffectiveness in infecting foreigners, and 100% mortality of infected people.  

The first 40 days of the epidemic progressed as visualized below:  

### Methods

To find the best-fit model we used a combination of direct search and descent methods to estimate the following parameters:

	α  the probability of an infected individual infecting a susceptible person
	γ  the inverse of the mean number of days a person is infectious
	I0  the number of infected individuals at the start of the time series
  
To estimate the parameters we visually selected the minimum sum of square value corresponding to each parameter from the following sum of squares profiles; this is a descent method for finding the model from which the data deviates the least:


We then used these visually estimated parameter values in an optimization function in R that converged upon a more precise final estimate.  Using this function, we found the following parameter values for our model of the epidemic: α = 0.00007399, γ = 0.1887, and I0 =  0.9546.  Explanations for the fact that the number of infected individuals at the start of the time series is less than one include observation uncertainty caused by delayed recognition of symptoms and uncounted individuals, and density-dependence of disease dynamics.  

We input the parameters selected by R into our time series model of the disease.  The following three plots show how well our finalized biological model of susceptible, infected, and removed (deceased) Kinamalian people fit the acquired data from the first 40 days, with green lines indicating deviations of the data from the model:

### Results 

Our model predicted that the epidemic will run for 96 days before the number of infected individuals drops below one and all members of the population are either removed or no longer under threat of infection.  Ultimately, according to current data this epidemic is predicted to kill 5417 people and leave 736 survivors.  Fortunately, barring any increases in infectiousness of the disease, our model shows that the epidemic already hit its peak on day 40, when 1355 individuals were infected with the disease and 2331 people were dead.  The number of infected individuals should decrease until the epidemic runs its full course.

### Conclusions

From this model of the Kinamali epidemic we conclude that in an attempt to decrease the rate of infection, all healthy Kinamalians should hide away for the next month and prevent any contact with other people if possible.  However, since the infection rate already appears to be decreasing, this is not entirely necessary, and there is little else that may be done to contain the spread of the disease.  Though foreigners do not appear to be susceptible, foreign contact should still be limited as a preventative measure to spreading the disease elsewhere.  Research should be conducted to determine the reason why only Kinamalians seem to be affected by this deadly disease.
