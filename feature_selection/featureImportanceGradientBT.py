# https://scikit-learn.org/stable/auto_examples/ensemble/plot_gradient_boosting_regression.html#plot-feature-importance
# Author: Peter Prettenhofer <peter.prettenhofer@gmail.com>
#         Maria Telenczuk <https://github.com/maikia>
#         Katrina Ni <https://github.com/nilichen>
#
# License: BSD 3 clause

import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets, ensemble
from sklearn.inspection import permutation_importance
from sklearn.metrics import mean_squared_error
from sklearn.model_selection import train_test_split

class make_featureImportanceGTB(object):
    npath = 'ROIsCSV/ROIsCol8/*.csv'
    def __init__(self, nyear= 2022, aBacia= '741'):

        self.nyear = nyear
        self.nameBacia = aBacia
        print(f"iniciar with bacia {aBacia} an year {nyear} ")

        nameCSV = f'/{aBacia}_{nyear}_c1.csv'
        # /home/superusuario/Dados/mapbiomas/col8/features/        
        lst_pathCSV = glob.glob(self.npath)
        dirCSVs = [kk for kk in lst_pathCSV[:] if nameCSV in kk]

        print("processing the file => ", dirCSVs)
        dfROIs = pd.read_csv(dirCSVs[0])
        print("the table dfROIs Shape = ", dfROIs.shape)
        # removing unimportant columns of table files
        dfROIs = dfROIs.drop(['system:index', '.geo'], axis=1)    
        dfROIs = dfROIs[dfROIs['year'] == year]
        print("----> now shape is ", dfROIs.shape)
        self.columns = [kk for kk in dfROIs.columns]
        self.columns.remove('year')
        self.columns.remove('class')
        print("=============== All columns ================ \n ", columns)

        self.X_train, self.X_test, self.y_train, self.y_test = train_test_split(
                    dfROIs[self.columns], dfROIs['class'], stratify=dfROIs['class'], random_state=13)

        self.params = {
                "n_estimators": 150,
                "max_depth": 4,
                "min_samples_split": 15,
                "learning_rate": 0.01,
                "loss": "squared_error",
            }

    def fitModel(self):
        self.reg = ensemble.GradientBoostingRegressor(**self.params)
        reg.fit(self.X_train, self.y_train)

        mse = mean_squared_error(self.y_test, self.reg.predict(self.X_test))
        print("The mean squared error (MSE) on test set: {:.4f}".format(mse))

        self.plot_training_device()

    def plot_training_device(self):
        test_score = np.zeros((self.params["n_estimators"],), dtype=np.float64)
        for i, self.y_pred in enumerate(self.reg.staged_predict(self.X_test)):
            test_score[i] = mean_squared_error(self.y_test, self.y_pred)

        fig = plt.figure(figsize=(6, 6))
        plt.subplot(1, 1, 1)
        plt.title("Deviance")
        plt.plot(
            np.arange(self.params["n_estimators"]) + 1,
            self.reg.train_score_,
            "b-",
            label="Training Set Deviance",
        )
        plt.plot(
            np.arange(params["n_estimators"]) + 1, test_score, "r-", label="Test Set Deviance"
        )
        plt.legend(loc="upper right")
        plt.xlabel("Boosting Iterations")
        plt.ylabel("Deviance")
        fig.tight_layout()
        plt.show()

    def make_features_importance(self):
        feature_importance = self.reg.feature_importances_
        sorted_idx = np.argsort(feature_importance)
        pos = np.arange(sorted_idx.shape[0]) + 0.5

        self.plot_results_FeatImportance(feature_importance)

    def plot_results_FeatImportance(self, feature_import):
        fig = plt.figure(figsize=(12, 16))
        plt.subplot(1, 2, 1)
        plt.barh(pos, feature_import[sorted_idx], align="center")
        plt.yticks(pos, np.array(self.columns)[sorted_idx])
        plt.title("Feature Importance (MDI)")

        result = permutation_importance(
            self.reg, self.X_test, self.y_test, n_repeats=10, random_state=42, n_jobs=2
        )
        sorted_idx = result.importances_mean.argsort()
        plt.subplot(1, 2, 2)
        plt.boxplot(
            result.importances[sorted_idx].T,
            vert=False,
            labels=np.array(self.columns)[sorted_idx],
        )
        plt.title("Permutation Importance (test set)")
        fig.tight_layout()
        plt.show()


if __name__ == '__main__':
    nbacias = '741'
    year = 2000
       
    makeROIs_analises = make_featureImportanceGTB(nyear= year, aBacia= nbacias)
    makeROIs_analises.fitModel()

    print("make feature feature importance ")
    makeROIs_analises.make_features_importance()
