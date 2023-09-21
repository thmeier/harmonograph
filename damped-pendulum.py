import numpy as np
import scipy as sp
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.animation as animation

from collections import deque
from diffsys import DiffSystem
from params import Params

# TODO: move into helper.py
def history2segments(hx, hy):
    history2points = lambda hx, hy: np.array([hx, hy]).T.reshape(-1, 1, 2)
    points2segments = lambda pts: np.concatenate([pts[:-1], pts[1:]], axis=1)
    return points2segments(history2points(hx, hy))

class DampedPendulum(DiffSystem):
    """
    System of ordinary differential equations, i.e. a damped pendulum.
    """
    def __init__(self, length, mass, damping, gravity, stop_time, initial_conditions, origin=(0, 0), debug=False):
        """
        length :: float
            length of bar of pendulum, measured in [m].
        mass :: float
            mass of bob of pendulum, measured in [kg].
        damping :: float
            frictional coefficient, unitless, i.e. [1].
        gravity :: float
            gravitational force acting on the pendulum's bob, measured in [m/s^2].
            Note: unintuitively we should use the value 9.81 instead of -9.81.
        stop_time :: float
            how many seconds to simulate the system, measured in [s].
        initial_conditions :: (float, float)
            tuple holding the initial angle of the pendulum, measured in [rad] and the initial angular velocity, measured in [rad/s]
        origin :: ([1], [1])
        """
        self.debug = debug

        self.params = Params()
        self.params.l = length
        self.params.m = mass
        self.params.d = damping
        self.params.g = gravity

        self.pen_conditions = (0.5 * self.params.l, 0.3 * self.params.l)
        self.pen_length = np.sqrt(np.sum(np.power(self.pen_conditions, 2)))

        self.origin = origin
        self.initial_conditions = initial_conditions
        self.stop_time = stop_time
        self.time = np.linspace(0.0,  self.stop_time, self.stop_time * 100)
        self.dt = self.time[1] - self.time[0]

        # solve ODE with given initial value
        self.solution = self.solve()

    def ode(self, ys, time):
        """
        # reduce second-order ODE system to one of first-order
        theta'' = -d * theta' - (g/l) * sin(theta)

        # let new variables be as follows
        y1 = theta1'
        y2 = theta1

        # their derivative yields
        y1' = theta1''
        y2' = theta1'

        # perform substitution
        y1' = -d * y1 - (g/l) * sin(y2)
        y2' = y1
        """
        y1, y2 = ys
        return np.array([
            -self.params.d * y1 - (self.params.g/self.params.l) * np.sin(y2),
            y1
        ])

    def animate(self, trace_length=20):
        if trace_length < 0:
            trace_length = self.time.shape[-1]

        bob_angle = self.solution[0, :]
        pen_angle = np.array(bob_angle) + np.arctan2(*self.pen_conditions)

        # (x/y) position of the mass
        bob_pos_x = self.origin[0] + self.params.l * np.sin(bob_angle)
        bob_pos_y = self.origin[1] - self.params.l * np.cos(bob_angle)

        # (x/y) position of the pen
        pen_pos_x = self.origin[0] + self.pen_length * np.sin(pen_angle)
        pen_pos_y = self.origin[1] - self.pen_length * np.cos(pen_angle)

        fig = plt.figure(figsize=(10, 10))
        axlim = (np.fmax(*self.origin) + np.fmax(self.params.l, self.pen_length)) * 1.25
        ax = fig.add_subplot(
            autoscale_on=False, xlim=(-axlim, axlim), ylim=(-axlim, axlim)
        )
        ax.set_aspect('equal')

        if not self.debug:
            ax.axis('off')
        else:
            ax.grid()
            ax.tick_params(
                # affect x-axis and y-axis
                axis='both',
                # affect major and minor ticks
                which='both',
                # disable ticks
                top=False, left=False, bottom=False, right=False,
                # disable labels
                labeltop=False, labelleft=False, labelbottom=False, labelright=False,
            )

        pendulum_arm, = ax.plot([], [], 'o-w', lw=2)
        pen_arm, = ax.plot([], [], '-r', lw=2)
        # TODO: REMOVE
        pen_trace = mpl.collections.LineCollection(
            [],
            cmap=plt.get_cmap('YlOrRd'),
            norm=plt.Normalize(0, self.time.shape[-1]),
            lw=1, ls=':',
        )
        # TODO: seems not to be needed, maybe remove?
        #pen_trace.set_array(self.time)
        ax.add_collection(pen_trace)
        time_template = 'time = %.1fs'
        time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
        # TODO: REMOVE
        pen_history_x = deque(maxlen=trace_length)
        pen_history_y = deque(maxlen=trace_length)
        # TODO: UNTIL HERE

        def animate(i):
            if i == 0:
                pen_history_x.clear()
                pen_history_y.clear()
                pen_trace.set_segments([])

            pendulum_arm_curr_x = [self.origin[0], bob_pos_x[i]]
            pendulum_arm_curr_y = [self.origin[1], bob_pos_y[i]]

            projected_pen_length = self.pen_length * np.cos(np.arctan2(*self.pen_conditions))
            gamma = projected_pen_length / self.params.l

            pen_arm_curr_x = [self.origin[0] + gamma * np.diff(pendulum_arm_curr_x)[0], pen_pos_x[i]]
            pen_arm_curr_y = [self.origin[1] + gamma * np.diff(pendulum_arm_curr_y)[0], pen_pos_y[i]]

            pen_history_x.appendleft(pen_pos_x[i])
            pen_history_y.appendleft(pen_pos_y[i])

            pen_trace.set_segments(history2segments(
                pen_history_x, pen_history_y
            ))
            pen_trace.set_array(
                list(range(i, np.max(i-trace_length, 0), -1))
            )

            pen_arm.set_data(pen_arm_curr_x, pen_arm_curr_y)
            pendulum_arm.set_data(pendulum_arm_curr_x, pendulum_arm_curr_y)
            time_text.set_text(time_template % (i * self.dt))

            return pendulum_arm, pen_arm, pen_trace, time_text

        _ = animation.FuncAnimation(
            fig, animate, self.solution.shape[-1], interval=self.dt*1000, blit=True
        )
        plt.show()

if __name__ == '__main__':
    plt.style.use(['dark_background'])

    pendulum = DampedPendulum(
        length=2.0,
        mass=2.0,
        #damping=20.,             # over damped
        damping=0.4,              # under damped
        #damping=2 * np.sqrt(20), # critically damped, damping^2 = 4 * mass
        gravity=9.81,
        stop_time=10,
        initial_conditions=(
            np.deg2rad(-90),   # initial angle
            np.deg2rad(0),  # initial angular velocity
        ),
        origin=(1, 1),
        debug=True #False
    )

    pendulum.animate(trace_length=120)
