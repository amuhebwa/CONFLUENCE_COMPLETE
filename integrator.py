import numpy as np
from functools import partial
from scipy.optimize import minimize
import numpy as np
import pandas as pd
import shutil
import netCDF4 as netcdf
import os
import itertools
import code


class MeanFlow:
    """Mean flow discharge integrator."""

    def __init__(self, n, A0, H, W, S, dA, routing_table):
        self.nreaches = len(A0)
        self.A0 = np.array(A0)
        self.n = np.array(n)
        w = np.mean(W, axis=0)
        h = np.mean(H, axis=0)
        if dA is None:
            dA = np.array([(w[r] + W[np.argmin(H[:, r]), r]) / 2 * (h[r] - H[np.argmin(H[:, r]), r])
                           for r in range(self.nreaches)]).T
        self.data = (H, W, S, dA)
        self.rivs = self._riverTopology(routing_table)

    def _riverTopology(self, routing_table):
        """Create river topology from routing table. River topology is
        a dictionary with each entry being a list of upstream reaches for
        the specific reach."""
        rout = pd.read_csv(routing_table)
        rout = rout.dropna(axis=1)[['GridID', 'HydroID', 'NextDownID']]
        rivs = {k: [] for k in rout['GridID']}
        hyd = {row[1]['HydroID']: row[1]['GridID'] for row in rout.iterrows()}
        for row in rout.iterrows():
            if row[1]['NextDownID'] > 0 and row[1]['NextDownID'] in hyd:
                rivs[hyd[row[1]['NextDownID']]].append(hyd[row[1]['HydroID']])
        idx = {k: i for i, k in enumerate(rivs)}
        rivs = {i: [idx[l] for l in rivs[k]] for i, k in enumerate(rivs)}
        return rivs

    def objective(self, x):
        """Objective function for mean-flow integration."""
        n = np.array(x[1::2])
        A0 = np.array(x[0::2])
        H, W, S, dA = self.data
        w = np.mean(W, axis=0)
        h = np.mean(H, axis=0)
        da = np.mean(dA, axis=0)
        Q = 1 / n * (A0 + da) ** (5 / 3) * w ** (-2 / 3) * np.mean(S,
                                                                   axis=0) ** (1 / 2)
        error = [np.sqrt((Q[i] - np.sum(Q[j] for j in self.rivs[i])) ** 2) for i in self.rivs]
        return np.sum(error)

    def constraint(self, i, x):
        """Constrain discharge to increase downstream."""
        n = np.array(x[1::2])
        A0 = np.array(x[0::2])
        H, W, S, dA = self.data
        w = np.mean(W, axis=0)
        h = np.mean(H, axis=0)
        da = np.mean(dA, axis=0)
        Q = 1 / n * (A0 + da) ** (5 / 3) * w ** (-2 / 3) * np.mean(S,
                                                                   axis=0) ** (1 / 2)
        return Q[i] - np.sum([Q[j] for j in self.rivs[i]])

    def integrate(self):
        """Integrate discharge by forcing mean annual
        flow to increase downstream."""
        lbnds = [0.0] * (self.nreaches * 2)
        ubnds = [0.0] * (self.nreaches * 2)
        lbnds[1::2] = [0.01] * self.nreaches
        ubnds[1::2] = [0.09] * self.nreaches
        lbnds[::2] = np.abs(np.amin(self.data[-1], axis=0))
        ubnds[::2] = [1e4] * self.nreaches
        bnds = list(zip(lbnds, ubnds))
        x0 = [0.0] * (self.nreaches * 2)
        x0[1::2] = self.n
        x0[::2] = self.A0
        cons = [{
            'type': 'ineq',
            'fun': partial(self.constraint, i)
        } for i in range(self.nreaches - 1)]
        solution = minimize(self.objective,
                            x0,
                            method='trust-constr',
                            bounds=bnds,
                            tol=1e-3,
                            constraints=cons)
        return solution.x[::2], solution.x[1::2]


def write(ncfile, Ab, n, Q_int):
    """Make copy of input file and write with estimated parameters."""
    shutil.copyfile(ncfile, "integrator.nc")
    f = netcdf.Dataset("integrator.nc", 'r+')

    # 1. save the new A0 produced by the integrator
    Ab = Ab.reshape((-1, 1))
    f.createDimension('Ab_row', Ab.shape[0])
    f.createDimension('Ab_col', Ab.shape[1])
    Ab_data = f.createVariable('Ab_Intgrator', np.float32, ('Ab_row', 'Ab_col'))
    Ab_data[:] = Ab

    # 2. save the new n produced by the integrator
    nInt = n.reshape((-1, 1))
    f.createDimension('nInt_row', nInt.shape[0])
    f.createDimension('nInt_col', nInt.shape[1])
    new_n_data = f.createVariable('n_Integrator', np.float32, ('nInt_row', 'nInt_col'))
    new_n_data[:] = nInt

    # 3. save the discharge produced by integrator
    f.createDimension('QInt_row', Q_int.shape[0])
    f.createDimension('QInt_col', Q_int.shape[1])
    QInt_data = f.createVariable('Q_Integrator', np.float32, ('QInt_row', 'QInt_col'))
    QInt_data[:] = Q_int

    f.close()


def removeFlagged(H, W, S, dA):
    """Clean up data."""
    i = np.where(np.all(H > -1000, axis=1))[0]
    H = H[i, :]
    i = np.where(np.all(W > -1000, axis=1))[0]
    W = W[i, :]
    i = np.where(np.all(S > -1000, axis=1))[0]
    S = S[i, :]
    if dA is not None:
        dA = dA[i, :]
    return H, W, S, dA


def main(ncfile, routing_table):
    """Driver function for the Confluence framework."""
    f = netcdf.Dataset(ncfile)
    H = f.variables['H'][:].data
    W = f.variables['W'][:].data
    S = f.variables['S'][:].data
    n = f.variables['n'][:].data
    # Qpred = f.variables["Qpred"][:].data
    # Qstar = f.variables["Qstar"][:].data
    Ab = f.variables['Abase'][:].data

    if 'dA' in f.variables:
        dA = f.variables['dA'][:].data
    else:
        dA = None
    
    print("Original Dimensions are; ")
    print('H: {},W: {}, S: {}, dA: {}'.format(H.shape, W.shape, S.shape, dA.shape))

    tempS = S.copy()
    tempIndex = np.where(np.all(tempS > -1000, axis=1))[0]

    H, W, S, dA = removeFlagged(H, W, S, dA)
    #code.interact(local=locals())
    print("New Dimensions are; ")
    print('H: {},W: {}, S: {}, dA: {}'.format(H.shape, W.shape, S.shape, dA.shape))

    mf = MeanFlow(n, Ab, H, W, S, dA, routing_table)
    Ab_est, n_est = mf.integrate()
    Q_int = 1 / n_est * (Ab_est + dA) ** (5 / 3) * W ** (-2 / 3) * np.mean(S, axis=0) ** (1 / 2)

    #code.interact(local=locals())
    Q_int_final = tempS
    Q_int_final[tempIndex, :] = Q_int
    #code.interact(local=locals())
    #Q_int = Q_int.reshape((1, -1))

    write(ncfile, Ab_est, n_est, Q_int_final)
    #write(ncfile, Ab_est, n_est, Q_int)
    # code.interact(local=locals())


fname = os.getcwd() + "/bam_output.nc"
routing_table = os.getcwd() + "/routing_table.csv"
main(fname, routing_table)
