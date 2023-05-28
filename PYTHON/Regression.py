import matplotlib.pyplot as plt
import numpy as np

# Measurement points
x = [1, 1.5, 1.3, 2, 2, 3, 3.5, 4, 4.1, 4, 5, 5.5, 6, 6, 6.5, 7, 7, 7.2, 8, 8.1]
y = [1, 2, 1.7, 3, 3.4, 2.7, 3, 3.5, 5, 4.2, 6, 7.2, 5, 7.9, 8.5, 9.5, 9, 10.5, 8.4, 10]

# Polynomial regression
degree = 2
coefficients = np.polyfit(x, y, degree)
polynomial = np.poly1d(coefficients)

# Plotting the regression line and measurement points
plt.scatter(x, y, color='red', label='Measurements')
plt.plot(x, polynomial(x), color='blue', label='Regression Line')
plt.xlabel('Redness level')
plt.ylabel('Ripening Stage')
plt.title('') #Ripening Stage in function of the redness of the sample
plt.legend()
plt.grid(False)
plt.show()
