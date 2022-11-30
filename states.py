from fractions import Fraction
from conservation import projection


class State:
    def __init__(self, state: str):

        self.__state = state
        self.__J = Fraction(int(list(state)[2]), 2)
        self.__I = Fraction(int(list(state)[1]), 2)
        self.__L = ['S', 'P', 'D', 'F', 'G', 'H', 'I'].index(list(state)[0])
        self.__S = Fraction(1, 2)

        self.__MJ = projection(self.J)
        self.__MI = projection(self.I)
        self.__ML = projection(self.L)
        self.__MS = projection(self.S)

    @property
    def state(self):
        return self.__state

    @property
    def J(self):
        return self.__J

    @property
    def I(self):
        return self.__I

    @property
    def L(self):
        return self.__L

    @property
    def S(self):
        return self.__S

    @property
    def MJ(self):
        return self.__MJ

    @property
    def MI(self):
        return self.__MI

    @property
    def ML(self):
        return self.__ML

    @property
    def MS(self):
        return self.__MS

    def __repr__(self):
        return f'State({self.state})'

    def __str__(self):
        return f'''State: {self.state}

        S = {self.S}      MS = {self.MS}
        L = {self.L}      ML = {self.ML}
        J = {self.J}      MJ = {self.MJ}
        I = {self.I}      MI = {self.MI}
                        '''
