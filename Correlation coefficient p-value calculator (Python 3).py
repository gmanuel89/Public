############ Load the libraries
import math
from math import *
from scipy.stats import t

################################################ INPUT
correlation_coefficient = float(input("Correlation coefficient: \n"))
correlation_type = input("Corrrelation type (default: 'pearson')\n")
number_of_samples = int(input("Number of samples:\n"))

try:
    tails = int(input("One tail or two tails? (default: 2)\n"))
except:
    tails = 2

try:
    level_of_significance = float(input("Level of significance (default: 0.05)\n"))
except:
    level_of_significance = 0.05

if correlation_type == '':
    correlation_type = "pearson"


################################### Calculation of the Student's t-distribution
degrees_of_freedom = number_of_samples-2

if correlation_type == "pearson":
    t_value = correlation_coefficient * sqrt((number_of_samples-2)/(1-correlation_coefficient**2))

# Calculate the one-tail p-value
p_value = t.sf(t_value, degrees_of_freedom)

if tails == 1:
    print ("The p-value for a %s correlation coefficient (of %s) computed on %s samples is: %s" %(correlation_type, correlation_coefficient, number_of_samples, p_value))

if tails == 2:
    p_value = p_value*2
    print ("The p-value for a %s correlation coefficient (of %s) computed on %s samples is: %s" %(correlation_type, correlation_coefficient, number_of_samples, p_value))

###################### Significance
if p_value <= level_of_significance:
    print ("The calculated correlation coefficient IS statistically significant at a level of significance of %s" %(level_of_significance))
else:
    print ("The calculated correlation coefficient is NOT statistically significant at a level of significance of %s" %(level_of_significance))
