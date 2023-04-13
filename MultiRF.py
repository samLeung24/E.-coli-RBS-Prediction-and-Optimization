import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.inspection import permutation_importance
from sklearn.model_selection import cross_val_score
from sklearn.metrics import mean_squared_error
import matplotlib.pyplot as plt

# One_hot_encode function
def one_hot_encode(seq):
    seq = seq.upper()
    nucleotides = 'ATCG'
    encoding = np.zeros((len(seq), len(nucleotides)))
    for i in range(len(seq)):
        if seq[i] in nucleotides:
            encoding[i, nucleotides.index(seq[i])] = 1
    return encoding.flatten()

def encode_structure(structure):
    n = len(structure)
    encoding = np.zeros((n, n), dtype=int)
    stack = []
    for i in range(n):
        if structure[i] == '(':
            stack.append(i)
        elif structure[i] == ')':
            j = stack.pop()
            encoding[i, j] = 1
            encoding[j, i] = 1
    return encoding.flatten()

def encode_structure2(structure):
    mapping = {'.': 0, '(': 1, ')': 2}
    encoding = np.zeros((len(structure), len(mapping)))
    stack = []
    for i, s in enumerate(structure):
        encoding[i, mapping[s]] = 1
    return encoding.flatten()

# Import my sequence file
df = pd.read_csv('./E.-coli-RBS-Prediction-and-Optimization/forRF5.csv')
df['encodedRBS'] = df['RBS'].apply(one_hot_encode)
df['encodedSC'] = df['start_codon'].apply(one_hot_encode)
df['loggedExp'] = np.log10(df['explev']+0.0001)
df['encodedStructure'] = df['X30nt'].apply(encode_structure2)
df['loggedTPM'] = np.log10(df['TPM']+0.0001)

# The variables
X = np.concatenate([
    df['cai'].values.reshape(-1, 1),  
    df['gc_content'].values.reshape(-1, 1),
    df['SC_pairing'].values.reshape(-1, 1),   
    df['encodedRBS'].tolist(),
    df['encodedSC'].tolist(),
    df['encodedStructure'].tolist(),             
    df['X30nt_score'].values.reshape(-1, 1),
    df['X30nt_proportion_paired'].values.reshape(-1, 1),
    df['loggedTPM'].values.reshape(-1, 1),
    df['Complim'].values.reshape(-1, 1)
], axis=1)
y_log = np.array(list(df['loggedExp']))
X_train, X_test, y_train, y_test = train_test_split(X, y_log, test_size=0.2, random_state=42)
    
# Split the test dataset into validation and hold-out test sets
X_val, X_holdout, y_val, y_holdout = train_test_split(X_test, y_test, test_size=0.5, random_state=42)

# Initialize variables to store the best model and its performanceß
best_rf = None
best_r = 0
best_i = 0
coefficients = []

# Initialize variables to store the RMSE values
RMSE_train = []
RMSE_val = []
max_depth_range = range(1, 35)

# Search for the best model for each depth using a for loop
for i in max_depth_range:
    # Train a new random forest model with different random seeds
    rf = RandomForestRegressor(n_estimators=1000, random_state=42, max_depth=i, n_jobs=10)
    rf.fit(X_train, y_train)
    y_train_pred = rf.predict(X_train)
    r_train = np.corrcoef(y_train, y_train_pred)[0][1]
    coefficients.append(r_train)

    # Evaluate the model on the validation set
    y_val_pred = rf.predict(X_val)
    r_val = np.corrcoef(y_val, y_val_pred)[0, 1]
    
    # Record the losses
    RMSE_train.append(np.sqrt(mean_squared_error(y_train, y_train_pred)))
    RMSE_val.append(np.sqrt(mean_squared_error(y_val, y_val_pred)))

    # Update the best model if the current model performs better
    if r_val > best_r:
        best_r = r_val
        best_rf = rf
        best_i = i

print(coefficients)
print(f"The optimal depth is {best_i}")
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

# Plot the loss function of the train and validation set
# Plot RMSE values against max_depth
plt.plot(max_depth_range, RMSE_train, label='Train')
plt.plot(max_depth_range, RMSE_val, label='Validation')
plt.xlabel('max_depth')
plt.ylabel('RMSE')
plt.legend()
plt.show()
