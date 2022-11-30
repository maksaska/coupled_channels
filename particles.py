from fractions import Fraction
from conservation import projection, relevant_particle
import numpy as np


class Particle:
    def __init__(self, name: str):
        par = relevant_particle(name)

        assert par['check'] == True, f'Particle name: {name} is invalid'

        self.__name = name

        self.__J = par['J']
        self.__I = par['I']
        self.__MI = par['I3']
        self.__MJ = projection(self.J)

        self.__mass = par['mass']
        self.__charge = par['charge']
        self.__strange = par['strange']
        self.__P = par['parity']
        

    @property
    def name(self):
        return self.__name

    @property
    def J(self):
        return self.__J

    @property
    def I(self):
        return self.__I

    @property
    def MJ(self):
        return self.__MJ

    @property
    def MI(self):
        return self.__MI

    @property
    def mass(self):
        return self.__mass

    @property
    def charge(self):
        return self.__charge

    @property
    def strange(self):
        return self.__strange

    @property
    def P(self):
        return self.__P

    def __repr__(self):
        return f'Particle({self.name})'

    def __str__(self):
        return f'''Particle: {self.name} 
        
        Charge: {self.charge}
        Strange: {self.strange}

        Mass = {self.mass} GeV 
        J = {self.J}   MJ = {self.MJ}
        I = {self.I}   MI = {self.MI}
        P = {self.P}
                        '''
