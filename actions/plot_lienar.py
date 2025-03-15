import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score

# Load the data
data = pd.read_csv('data.csv')

# Extract the columns
X = data['De'].values.reshape(-1, 1)
y = data['Ea'].values

# Perform linear regression
model = LinearRegression()
model.fit(X, y)
y_pred = model.predict(X)

# Calculate metrics
mae = mean_absolute_error(y, y_pred)
rmse = np.sqrt(mean_squared_error(y, y_pred))
r2 = r2_score(y, y_pred)
#equation = f'Ea = {model.coef_[0]:.2f} $\Delta$E + {model.intercept_:.2f}'
equation = r'Ea = {:.2f} $\Delta$E + {:.2f}'.format(model.coef_[0], model.intercept_)

# Plot the data and regression line
plt.figure(figsize=(10, 6))
plt.scatter(X, y, color='blue', label='DFT')
plt.plot(X, y_pred, color='red', label='Reg.')

# Add text box for metrics
textstr = f'{equation}\nMAE: {mae:.2f}\nRMSE: {rmse:.2f}\nR2: {r2:.2f}'
plt.gca().text(0.65, 0.25, textstr, transform=plt.gca().transAxes,
               fontsize=16, verticalalignment='top', bbox=dict(facecolor='white', alpha=0))

# Customize the plot
plt.xlabel(r'$\Delta$E /eV', fontsize=16)
plt.ylabel('Ea /eV', fontsize = 16)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.legend(loc='upper left', frameon=False, fontsize=18)
#plt.title('Linear Regression of De vs Ea')
plt.tight_layout()
# Save the figure
plt.savefig('data.png')

plt.show()
