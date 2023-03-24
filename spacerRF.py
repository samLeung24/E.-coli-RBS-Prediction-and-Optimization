import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from sklearn.inspection import permutation_importance
from sklearn.model_selection import cross_val_score
import matplotlib.pyplot as plt
from tensorflow.keras.preprocessing.sequence import pad_sequences

# One_hot_encode function
def one_hot_encode(seq):
    seq = seq.upper()
    nucleotides = 'AUCG'
    encoding = np.zeros((len(seq), len(nucleotides)))
    for i in range(len(seq)):
        if seq[i] in nucleotides:
            encoding[i, nucleotides.index(seq[i])] = 1
    return encoding.flatten()

# Padding function
def pad_nucleotide_sequences(seq_list, maxlen):
    padded_sequences = []
    for seq in seq_list:
        seq = seq.upper()
        padded_seq = seq.ljust(maxlen, 'N')
        padded_sequences.append(padded_seq)
    return padded_sequences

# Define RNA nucleotides and their corresponding indices
nucleotides = ['A', 'C', 'G', 'U']
indices = np.identity(4)

# Import my sequence file
df = pd.read_csv('./E.-coli-RBS-Prediction-and-Optimization/positiveProdigalForRF.csv')
df = df[df['spacer_seq'].apply(lambda x: isinstance(x, str))]
maxlen = df['spacer_seq'].apply(len).max()
df['paddedSD'] = pad_nucleotide_sequences(df['spacer_seq'], maxlen=maxlen)
df['encodedSD'] = df['paddedSD'].apply(one_hot_encode)
df['loggedExp'] = np.log10(df['V3'])
mask = np.isfinite(df['loggedExp'])
df = df.loc[mask, :]

# Split data into training and testing sets
X = np.array(list(df['encodedSD']))
y_log = np.array(list(df['loggedExp']))
X_train, X_test, y_train, y_test = train_test_split(X, y_log, test_size=0.2, random_state=40)

# Train random forest model
rf = RandomForestRegressor(n_estimators=5000, random_state=42, max_depth=10, n_jobs=10)
rf.fit(X_train, y_train)

# Compute cross-validation scores
scores = cross_val_score(rf, X, y_log, cv=5, scoring='neg_mean_squared_error')

# Print mean and standard deviation of scores
print("Mean squared error:", -scores.mean())
print("Standard deviation:", scores.std())
# Evaluate the model on the testing set
from sklearn.metrics import mean_squared_error
y_pred = rf.predict(list(X_test))
mse = mean_squared_error(y_test, y_pred)
print("Mean squared error:", mse)

# Calculate and print median and IQR of residuals
y_test_flat = np.ravel(y_test)
y_pred_flat = np.ravel(y_pred)
residuals = y_pred_flat - y_test_flat
residuals_finite = residuals[np.isfinite(residuals)]
print("Median of residuals:", np.median(residuals_finite))
print("IQR of residuals:", np.percentile(residuals_finite, 75) - np.percentile(residuals_finite, 25))

# Plot residuals
plt.scatter(y_pred_flat, y_pred_flat - y_test_flat, c="b", s=40, alpha=0.5)
plt.hlines(y=0, xmin=y_pred_flat.min(), xmax=y_pred_flat.max())
plt.title("Residual Plot")
plt.xlabel("Predicted Values")
plt.ylabel("Residuals")
plt.show()

# Plot actual vs. predict (train)
plt.figure()
plt.scatter(y_train, rf.predict(X_train), alpha=0.2)
plt.plot([y_train.min(), y_train.max()], [y_train.min(), y_train.max()], 'r--')
plt.xlabel('Actual expression level')
plt.ylabel('Predicted expression level')
plt.title('Actual vs. predicted (train)')
# Calculate and display correlation coefficient on plot
y_train_pred = rf.predict(X_train)
r = np.corrcoef(y_train, y_train_pred)[0][1]
plt.text(0.95, 0.95, f"R={r:.3f}", ha='right', va='top', transform=plt.gca().transAxes,fontsize=12)
plt.show()

# Plot actual vs. predicted (test)
plt.figure()
plt.scatter(y_test, y_pred, alpha=0.2)
plt.plot([y_test.min(), y_test.max()], [y_test.min(), y_test.max()], 'r--')
plt.xlabel('Actual expression level')
plt.ylabel('Predicted expression level')
plt.title('Actual vs. predicted (test)')
# Add correlation coefficient to plot
R = np.corrcoef(y_test, y_pred)[0, 1]
plt.text(0.95, 0.95, f"R={R:.3f}", ha='right', va='top', transform=plt.gca().transAxes,fontsize = 12)
plt.show()

