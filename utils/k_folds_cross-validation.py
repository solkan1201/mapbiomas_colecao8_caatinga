from sklearn.datasets import load_boston
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LinearRegression
from sklearn.model_selection import cross_val_score



"""
    K-fold Cross Validation (The Wrong Way)

        The k-fold cross-validation only works when the models are trained 
    only with data they should have access to. This rule can be violated 
    if the data is processed improperly prior to the sampling.
"""

