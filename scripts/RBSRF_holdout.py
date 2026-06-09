import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
import matplotlib.pyplot as plt
import csv

# One_hot_encode function
def one_hot_encode(seq):
    seq = seq.upper()
    nucleotides = 'AUCG'
    encoding = np.zeros((len(seq), len(nucleotides)))
    for i in range(len(seq)):
        if seq[i] in nucleotides:
            encoding[i, nucleotides.index(seq[i])] = 1
    return encoding.flatten()

# Import my sequence file
df = pd.read_csv('./E.-coli-RBS-Prediction-and-Optimization/positiveProdigalForRF.csv')
df['encodedRBS'] = df['RBS'].apply(one_hot_encode)
df['loggedExp'] = np.log10(df['V3'])
mask = np.isfinite(df['loggedExp'])
df = df.loc[mask, :]

# Split data into training and testing sets
X = np.array(list(df['encodedRBS']))
y_log = np.array(list(df['loggedExp']))
X_train, X_test, y_train, y_test = train_test_split(X, y_log, test_size=0.2, random_state=42)

# Split the test dataset into validation and hold-out test sets
X_val, X_holdout, y_val, y_holdout = train_test_split(X_test, y_test, test_size=0.5, random_state=42)

# Initialize variables to store the best model and its performance
best_rf = None
best_r = 0
coefficients = []
# Search for the best model using a for loop
for i in range(50):
    # Train a new random forest model with different random seeds
    rf = RandomForestRegressor(n_estimators=1000, random_state=i, max_depth=12, n_jobs=8)
    rf.fit(X_train, y_train)
    y_train_pred = rf.predict(X_train)
    r_train = np.corrcoef(y_train, y_train_pred)[0][1]
    coefficients.append(r_train)

    # Evaluate the model on the validation set
    y_val_pred = rf.predict(X_val)
    r_val = np.corrcoef(y_val, y_val_pred)[0, 1]
   

    # Update the best model if the current model performs better
    if r_val > best_r:
        best_r = r_val
        best_rf = rf

print(coefficients)

# Evaluate the validation set
y_val_pred = best_rf.predict(X_val)
R_val = np.corrcoef(y_val, y_val_pred)[0, 1]

# Plot actual vs. predicted (validation set)
plt.figure()
plt.scatter(y_val, y_val_pred, alpha=0.2)
plt.plot([y_val.min(), y_val.max()], [y_val.min(), y_val.max()], 'r--')
plt.xlabel('Actual expression level')
plt.ylabel('Predicted expression level')
plt.title('Actual vs. predicted (validation set)')
plt.text(0.95, 0.95, f"R={R_val:.3f}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
plt.show()

# Evaluate the best model on the hold-out test set
y_holdout_pred = best_rf.predict(X_holdout)
R_holdout = np.corrcoef(y_holdout, y_holdout_pred)[0, 1]

# Plot actual vs. predicted (hold-out test)
plt.figure()
plt.scatter(y_holdout, y_holdout_pred, alpha=0.2)
plt.plot([y_holdout.min(), y_holdout.max()], [y_holdout.min(), y_holdout.max()], 'r--')
plt.xlabel('Actual expression level')
plt.ylabel('Predicted expression level')
plt.title('Actual vs. predicted (hold-out set)')
plt.text(0.95, 0.95, f"R={R_holdout:.3f}", ha='right', va='top', transform=plt.gca().transAxes, fontsize=12)
plt.show()

# Calculate and print median and IQR of residuals and plot the residuals
y_test_flat = np.ravel(y_holdout)
y_pred_flat = np.ravel(y_holdout_pred)
residuals = y_pred_flat - y_test_flat
residuals_finite = residuals[np.isfinite(residuals)]
print("Median of residuals:", np.median(residuals_finite))
print("IQR of residuals:", np.percentile(residuals_finite, 75) - np.percentile(residuals_finite, 25))
plt.scatter(y_pred_flat, y_pred_flat - y_test_flat, c="b", s=40, alpha=0.5)
plt.hlines(y=0, xmin=y_pred_flat.min(), xmax=y_pred_flat.max())
plt.title("Residual Plot (holdout set)")
plt.xlabel("Predicted Values")
plt.ylabel("Residuals")
plt.show()