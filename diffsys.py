import scipy as sp

from abc import ABC, abstractmethod


class DiffSystem(ABC):
    """
    System of ordinary differential equations.
    """
    @abstractmethod
    def __init__(self, *args):
        pass

    @abstractmethod
    def ode(self, ys, time):
        pass

    def solve(self):
        return sp.integrate.odeint(
            self.ode, self.initial_conditions, self.time
        ).T

    @abstractmethod
    def animate(self, trace_length=20):
        pass
