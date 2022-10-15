# import maninetcluster here
import numpy as np
import pandas as pd
from sklearn import preprocessing

from ManinetCluster_python.alignment import CCA
from ManinetCluster_python.correspondence import Correspondence
from ManinetCluster_python.distance import SquaredL2
from ManinetCluster_python.neighborhood import neighbor_graph

class CCA_Alignment:
    def __init__(self,file1,file2,d):
        self.file1 = file1
        self.file2 = file2
        self.d = d
        
        
    def run(self):
        X = pd.read_csv(self.file1)
        #X = X.T
        #X.index = X["Unnamed: 0_x"]
        #X = X.drop(["Unnamed: 0_x"], axis=1)
        
        Y = pd.read_csv(self.file2)
        #Y = Y.T
        #Y.index = Y["Unnamed: 0"]
        #Y = Y.drop(["Unnamed: 0"], axis=1)
        
        X = X.to_numpy()
        Y = Y.to_numpy()
        
        n = len(X)
        d = self.d
        
        #actual = generate_correlation_map(X, Y)
        actual = np.eye(n)
        
        corr = Correspondence(matrix=actual)
        
        X_normalized = preprocessing.normalize(X, norm='l2')
        Y_normalized = preprocessing.normalize(Y, norm='l2')

        cca = CCA(X_normalized, Y_normalized, corr, d)
        Xnew,Ynew = cca.project(X_normalized, Y_normalized)
        Xnew_df = pd.DataFrame(Xnew)
        Ynew_df = pd.DataFrame(Ynew)
        
        Ynew_df.to_csv("cca_projected_points_mod2_"+str(Y.shape[1])+"_d"+str(d)+".csv", index=False)
        Xnew_df.to_csv("cca_projected_points_mod1_"+str(X.shape[1])+"_d"+str(d)+".csv", index=False)     
        
def generate_correlation_map(x, y):
    """Correlate each n with each m.

    Parameters
    ----------
    x : np.array
      Shape N X T.

    y : np.array
      Shape M X T.

    Returns
    -------
    np.array
      N X M array in which each element is a correlation coefficient.

    """
    mu_x = x.mean(1)
    mu_y = y.mean(1)
    n = x.shape[1]
    if n != y.shape[1]:
        raise ValueError('x and y must ' +
                         'have the same number of timepoints.')
    s_x = x.std(1, ddof=n - 1)
    s_y = y.std(1, ddof=n - 1)
 
    cov = np.dot(x,
                 y.T) - n * np.dot(mu_x[:, np.newaxis],
                                  mu_y[np.newaxis, :])
    

    deno = np.dot(s_x[:, np.newaxis], s_y[np.newaxis, :])
    final = cov / np.dot(s_x[:, np.newaxis], s_y[np.newaxis, :])
    final[np.isnan(final)] = 0    

    return final
#     return cov / np.dot(s_x[:, np.newaxis], s_y[np.newaxis, :])
        
        
