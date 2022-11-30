from particles import Particle
from states import State
from conservation import channel_check
from scipy.interpolate import interp1d
from scipy.special import sph_harm as Y
from sympy.physics.quantum.cg import CG
import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import numpy as np
import glob


class Interaction:
    def __init__(self, name_init: str, name_fin: str):
        par = channel_check(name_init, name_fin)
        assert 'None' not in par.items(), f'Particles in {name_init} to {name_fin} are invalid'

        self.__B_init = Particle(par['B_init'])
        self.__M_init = Particle(par['M_init'])

        self.__B_fin = Particle(par['B_fin'])
        self.__M_fin = Particle(par['M_fin'])

        assert self.B_init.MI + self.M_init.MI == self.B_fin.MI + \
            self.M_fin.MI, f'Channel: {name_init} to {name_fin} is invalid (I3 does not conserve)'
        assert self.B_init.charge + self.M_init.charge == self.B_fin.charge + \
            self.M_fin.charge, f'Channel: {name_init} to {name_fin} is invalid (Charge does not conserve)'
        assert self.B_init.strange + self.M_init.strange == self.B_fin.strange + \
            self.M_fin.strange, f'Channel: {name_init} to {name_fin} is invalid (Strangeness does not conserve)'

        self.__name_init = name_init
        self.__name_fin = name_fin
        self.__Data = {}

        path_init = par['path_init']
        path_fin = par['path_fin']

        for i in glob.glob(f'./Data/{path_init}/{path_fin}/*.csv'):
            self.__Data[i.split('/')[4].split('.')[0]] = pd.read_csv(i)
            self.__Data[i.split('/')[4].split('.')[0]
                        ]['W'] = self.__Data[i.split('/')[4].split('.')[0]]['W']/1000

    @ property
    def name_init(self):
        return self.__name_init

    @ property
    def name_fin(self):
        return self.__name_fin

    @ property
    def Data(self):
        return self.__Data

    @ property
    def B_init(self):
        return self.__B_init

    @ property
    def M_init(self):
        return self.__M_init

    @ property
    def B_fin(self):
        return self.__B_fin

    @ property
    def M_fin(self):
        return self.__M_fin

    def __repr__(self):
        return f'Interaction({self.name_init}, {self.name_fin})'

    def __str__(self):
        print(f'Channel: {self.name_init} to {self.name_fin}\n')
        print(self.B_init)
        print(self.M_init)
        print(self.B_fin)
        print(self.M_fin)
        print(self.Data.keys())
        for i in self.Data.keys():
            print(f'\nMultipole(L_2I2J): {i}')
            print(self.Data[i])
            print('For initial state: ')
            print(State(i))
            print('For final state: ')
            print(State(i))
        return ''

    def E1(self, W: float, init=False, fin=False) -> float:
        if init:
            return (W**2 + self.M_init.mass**2 - self.B_init.mass**2)/2/W
        if fin:
            return (W**2 + self.M_fin.mass**2 - self.B_fin.mass**2)/2/W
        return 0

    def E2(self, W: float, init=False, fin=False) -> float:
        if init:
            return (W**2 + self.B_init.mass**2 - self.M_init.mass**2)/2/W
        if fin:
            return (W**2 + self.B_fin.mass**2 - self.M_fin.mass**2)/2/W
        return 0

    def k(self, W: float, init=False, fin=False) -> float:
        if init:
            return np.sqrt(((W**2 + self.B_init.mass**2 - self.M_init.mass**2)/2/W)**2 - self.B_init.mass**2)
        if fin:
            return np.sqrt(((W**2 + self.B_fin.mass**2 - self.M_fin.mass**2)/2/W)**2 - self.B_fin.mass**2)
        return 0

    def rho(self, W: float, init=False, fin=False) -> float:
        return np.pi*self.k(W, init, fin)*self.E1(W, init, fin)*self.E2(W, init, fin)/W

    def t(self, W: float, state: str) -> complex:
        f_real = interp1d(self.Data[state]['W'], self.Data[state]['Re_T'], kind='cubic')
        f_imag = interp1d(self.Data[state]['W'], self.Data[state]['Im_T'], kind='cubic')
        return -(f_real(W) + f_imag(W)*1j)/np.sqrt(self.rho(W, init = True))/np.sqrt(self.rho(W, fin = True))

    def ampl2(self, phase_space, mj: dict):

        W = phase_space['W'][0]

        result = pd.DataFrame()
        result = phase_space.copy()

        result['T'] = 0 + 0*1j

        MI = self.M_init.MI + self.B_init.MI        

        for i in self.Data.keys():
            t_ampl = self.t(W, i)
            state = State(i)

            clebsh_isospin_in = float(CG(self.M_init.I, self.M_init.MI, self.B_init.I, self.B_init.MI, state.I, MI).doit())
            clebsh_isospin_fin = float(CG(self.M_fin.I, self.M_fin.MI, self.B_fin.I, self.B_fin.MI, state.I, MI).doit())

            for ML_in in state.ML:
                sph_harm_in = Y(ML_in, state.L, 0, 0)   # Y(ml, l, phi, theta)
                for MS_in in state.MS:
                    if mj['mj_M_init'] + mj['mj_B_init'] != MS_in:
                        continue
                    clebsh_part_in = float(CG(self.M_init.J, mj['mj_M_init'], self.B_init.J, mj['mj_B_init'], state.S, MS_in).doit())
                    for MJ in state.MJ:            # CG(j1, m1, j2, m2 | j3, m3) 
                        if ML_in + MS_in != MJ:
                            continue
                        clebsh_state_in = float(CG(state.L, ML_in, state.S, MS_in, state.J, MJ).doit())
                        for MS_fin in state.MS:
                            if mj['mj_M_fin'] + mj['mj_B_fin'] != MS_fin:
                                continue
                            clebsh_part_fin = float(CG(self.M_fin.J, mj['mj_M_fin'], self.B_fin.J, mj['mj_B_fin'], state.S, MS_fin).doit()) 
                            for ML_fin in state.ML: 
                                if ML_fin + MS_fin != MJ:
                                    continue                               
                                clebsh_state_fin = float(CG(state.L, ML_fin, state.S, MS_fin, state.J, MJ).doit()) 

                                theta_phi = {'Y':[]}

                                for _, row in phase_space.iterrows():
                                    theta_phi['Y'].append(Y(ML_fin, state.L, row['phi']*np.pi/180.0, row['theta']*np.pi/180.0).conj())
                                

                                result = pd.concat([result, pd.DataFrame(theta_phi)], axis=1)

                                result['T'] += t_ampl*clebsh_isospin_in*clebsh_state_in*clebsh_part_in*clebsh_isospin_fin*clebsh_state_fin*clebsh_part_fin*sph_harm_in*result['Y']
                                result.drop(columns=['Y'], inplace=True)

        result['amlitude2_inner'] = abs(result['T'])**2
        result.drop(columns=['W', 'theta', 'phi', 'T'], inplace=True)

        return result

    def diff_cs(self, phase_space):

        W = phase_space['W'][0]

        result = pd.DataFrame()
        result = phase_space.copy()

        result['amlitude2'] = 0

        GeV2_to_mbsr = 0.1973**2*10

        mj = {
            'mj_M_init' : 0,
            'mj_M_fin': 0,
            'mj_B_init': 0,
            'mj_B_fin': 0
        }

        for mj_M_init in self.M_init.MJ:
            for mj_B_init in self.B_init.MJ:
                for mj_M_fin in self.M_fin.MJ:
                    for mj_B_fin in self.B_fin.MJ:
                        if mj_M_init + mj_B_init != mj_M_fin + mj_B_fin:
                            continue

                        mj = {
                            'mj_M_init' : mj_M_init,
                            'mj_M_fin': mj_M_fin,
                            'mj_B_init': mj_B_init,
                            'mj_B_fin': mj_B_fin
                        }

                        result = pd.concat([result, self.ampl2(phase_space, mj)], axis=1)
                        result['amlitude2'] = result['amlitude2'] + result['amlitude2_inner']
                        result.drop(columns=['amlitude2_inner'], inplace=True)

        result['cs'] = GeV2_to_mbsr*(4*np.pi)**2/self.k(W, init=True)**2*self.rho(W, init=True)*self.rho(W, fin=True)/(2*self.M_init.J + 1)/(2*self.B_init.J + 1)*result['amlitude2']
        result.drop(columns=['amlitude2'], inplace=True)

        return result

    def check_interpolation(self, state: str):
        interp_x, interp_y = [], []

        W_min = self.Data[state]['W'].min()
        W_max = self.Data[state]['W'].max()

        for W in np.arange(W_min, W_max, 0.0025):
            interp_x.append(self.t(W, state).real)
            interp_y.append(self.t(W, state).imag)


        fig = px.scatter(x=self.Data[state]['W'].tolist(), y=self.Data[state]['Re_T'].tolist(), title=f'Real part of T for {state}')
        fig.add_trace(
            go.Scatter(
                x=np.arange(W_min, W_max, 0.0025),
                y=interp_x,
                mode="lines",
                line=go.scatter.Line(color="gray"),
                showlegend=False)
        )

        fig.update_layout(
            title=f'Real part of T for {state}',
            xaxis_title="W [GeV]",
            yaxis_title="Re(T)",
        )

        fig.show()

        fig = px.scatter(x=self.Data[state]['W'].tolist(), y=self.Data[state]['Im_T'].tolist(), title=f'Imaginary part of T for {state}')
        fig.add_trace(
            go.Scatter(
                x=np.arange(W_min, W_max, 0.0025),
                y=interp_y,
                mode="lines",
                line=go.scatter.Line(color="gray"),
                showlegend=False)
        )

        fig.update_layout(
            title=f'Imaginary part of T for {state}',
            xaxis_title="W [GeV]",
            yaxis_title="Im(T)",
        )

        fig.show()

    def check_diff_cs(self, W):

        X = np.arange(0, 180, 1)

        df = pd.DataFrame(X, columns =['theta'])
        df['W'] = W
        df['phi'] = 0

        Y = self.diff_cs(df)['cs'].tolist()

        fig = px.scatter(x=X, y=Y)
        fig.update_layout(
            title=f'Differential cross section for {self.M_init.name}{self.B_init.name} to {self.M_fin.name}{self.B_fin.name}; W = {W} GeV',
            xaxis_title=r'$\theta\;[degree]$',
            yaxis_title=r'$\dfrac{\text{d}\sigma_{\gamma^{*}}}{\text{d}\Omega} \text{ [mb]/sr}$',
        )

        fig.show()
