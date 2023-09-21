#!/usr/bin/python3

import numpy as np
import scipy as sp
import argparse as ap
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

# TODO: move into helper.py
def float_range(arg, val_min=0.0, val_max=1.0):
    """ Type function for argparse - a float within some predefined bounds """
    try:
        f = float(arg)
    except ValueError:
        raise argparse.ArgumentTypeError(f'Float expected, got {type(arg)}')
    if f < val_min or f > val_max:
        raise argparse.ArgumentTypeError(f'Argument must be < {val_max} and > {val_min}')
    return f


class Harmonograph(DiffSystem):
    """
    System of ordinary differential equations, i.e. coupled damped pendula.
    """
    def __init__(self,
                 lengths,
                 masses,
                 spring,
                 damping,
                 gravity,
                 stop_time,
                 initial_conditions,
                 pen_properties,
                 canvas_properties,
                 origin,
                 pendula_distance,
                 debug=False
                 ):
        self.debug = debug

        self.params = Params()
        self.params.ls = lengths
        self.params.ms = masses
        self.params.k = spring
        self.params.d = damping
        self.params.g = gravity

        self.pendula_distance = pendula_distance
        # TODO: make pen_properties a dict with 'offset', 'armlength'
        #       then make the field in the correct reihenfolge such that arctan can be used
        self.pen_properties = pen_properties
        self.pen_length = np.linalg.norm(self.pen_properties)

        # TODO: make canvas_properties a dict with 'offset', 'sidelength'
        self.canvas_properties = canvas_properties
        self.canvas_offset = canvas_properties[0]
        self.canvas_sidelength = canvas_properties[1]

        self.origin1 = origin
        self.origin2 = (origin[0] + self.pendula_distance, origin[1])
        self.initial_conditions = initial_conditions
        self.stop_time = stop_time
        self.time = np.linspace(0.0,  self.stop_time, self.stop_time * 100)
        self.dt = self.time[1] - self.time[0]

        # solve ODE with given initial value
        self.solution = self.solve()

    def ode(self, ys, time):
        """
        Damped Coupled Pendula
        """
        y00, y01, y10, y11 = ys
        m0, m1 = self.params.ms
        l0, l1 = self.params.ls
        g, k, d = self.params.g, self.params.k, self.params.d
        # NOTE: uses small angle approximation (sin phi = phi)
        # return np.array([
        #     ((m0 * (l0 * y00 ** 2 - g) - k * l0) * np.sin(y01) + k * l1 * np.sin(y11)) / (m0 * l0 * np.cos(y01)) - d / m0 * y00,
        #     y00,
        #     ((m1 * (l1 * y10 ** 2 - g) - k * l1) * np.sin(y11) + k * l0 * np.sin(y01)) / (m1 * l1 * np.cos(y11)) - d / m1 * y10,
        #     y10
        # ])
        # NOTE: does not use small angle approximation (sin phi = phi)
        return np.array([
            (y00 ** 2 + k / (m0 * l0) * (l1 * np.sin(y11) - l0 * np.sin(y01)) - g / l1) * np.tan(y01) - d / m0 * y00,
            y00,
            (y10 ** 2 + k / (m1 * l1) * (l0 * np.sin(y01) - l1 * np.sin(y11)) - g / l0) * np.tan(y11) - d / m1 * y10,
            y10
        ])

    def show_angle(self):
        angle_pendulum1 = self.solution[0, :]
        angle_pendulum2 = self.solution[2, :]

        fig, ax = plt.subplots(figsize=(5, 7))
        ax.plot(self.time, angle_pendulum1)
        ax.plot(self.time, angle_pendulum2)
        plt.show()

    def show_harmonogram(self):

        # CALCULATIONS

        # angle (w.r.t. time) of first pendulum
        angle_pendulum1 = self.solution[0, :]
        # angle (w.r.t. time) of second pendulum
        angle_pendulum2 = self.solution[2, :]
        # angle (w.r.t. time) of pencil
        angle_pen = np.array(angle_pendulum1) + np.arctan2(*self.pen_properties)

        # position of the pencil
        pen_pos_x = self.origin1[0] + self.pen_length * np.sin(angle_pen)
        pen_pos_y = self.origin1[1] - self.pen_length * np.cos(angle_pen)

        # position of the canvas
        canvas_center_pos_x = self.origin2[0] + self.canvas_offset * np.sin(angle_pendulum2),
        canvas_center_pos_y = self.origin2[1] - self.canvas_offset * np.cos(angle_pendulum2),

        # x/y and ys of curve
        harmonogram_x = (pen_pos_x - canvas_center_pos_x)
        harmonogram_y = (pen_pos_y - canvas_center_pos_y)

        # calculate nice limits for matplotlib
        scale = 1.1
        xlim = scale * np.array([np.min(harmonogram_x), np.max(harmonogram_x)])
        ylim = scale * np.array([np.min(harmonogram_y), np.max(harmonogram_y)])

        # PLOTTING

        fig, ax = plt.subplots(figsize=(10, 10))
        ax.set(
            autoscale_on=False, aspect='equal',
            xlim=xlim, ylim=ylim,
        )

        if not self.debug:
            ax.axis('off')
        else:
            ax.grid(color='darkgray', ls=':')
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

        harmonogram = mpl.collections.LineCollection(
            history2segments(harmonogram_x, harmonogram_y),
            array=np.arange(0, self.time.shape[-1]),
            # c.f. https://matplotlib.org/stable/tutorials/colors/colormaps.html
            # YlOrRd, inferno, bone, coolwarm, twilight, Pastel1
            cmap=plt.get_cmap('plasma'), alpha=0.85, lw=1, ls=':',
        )
        ax.add_collection(harmonogram)
        plt.show()

    def animate(self, trace_length=20):
        if trace_length < 0:
            trace_length = self.time.shape[-1]

        # PENDULA

        angle1 = self.solution[0, :]
        angle2 = self.solution[2, :]

        # position of the bob of pendulum1
        bob1_pos_x = self.origin1[0] + self.params.ls[0] * np.sin(angle1)
        bob1_pos_y = self.origin1[1] - self.params.ls[0] * np.cos(angle1)

        # position of the bob of pendulum2
        bob2_pos_x = self.origin2[0] + self.params.ls[1] * np.sin(angle2)
        bob2_pos_y = self.origin2[1] - self.params.ls[1] * np.cos(angle2)

        # PEN (TRACE)

        pen_angle = np.array(angle1) + np.arctan2(*self.pen_properties)

        # position of the pen
        pen_pos_x = self.origin1[0] + self.pen_length * np.sin(pen_angle)
        pen_pos_y = self.origin1[1] - self.pen_length * np.cos(pen_angle)

        # CANVAS

        _canvas_bottom_corner = (
            0.5 * self.canvas_sidelength, self.canvas_offset + 0.5 * self.canvas_sidelength
        )
        _canvas_top_corner = (
            0.5 * self.canvas_sidelength, self.canvas_offset - 0.5 * self.canvas_sidelength
        )

        canvas_angle = {
            'bottom_left' :
                np.array(angle2) + np.arctan2(*_canvas_bottom_corner),
            'bottom_right' :
                np.array(angle2) - np.arctan2(*_canvas_bottom_corner),
            'top_left' :
                np.array(angle2) + np.arctan2(*_canvas_top_corner),
            'top_right' :
                np.array(angle2) - np.arctan2(*_canvas_top_corner),
        }

        canvas_pos_x = {
            'center' :
                self.origin2[0] + self.canvas_offset * np.sin(angle2),
            'bottom_left' :
                self.origin2[0] + np.linalg.norm(_canvas_bottom_corner) * np.sin(canvas_angle['bottom_left']),
            'bottom_right' :
                self.origin2[0] + np.linalg.norm(_canvas_bottom_corner) * np.sin(canvas_angle['bottom_right']),
            'top_left' :
                self.origin2[0] + np.linalg.norm(_canvas_top_corner) * np.sin(canvas_angle['top_left']),
            'top_right' :
                self.origin2[0] + np.linalg.norm(_canvas_top_corner) * np.sin(canvas_angle['top_right']),
        }

        canvas_pos_y = {
            'center' :
                self.origin2[1] - self.canvas_offset * np.cos(angle2),
            'bottom_left' :
                self.origin2[1] - np.linalg.norm(_canvas_bottom_corner) * np.cos(canvas_angle['bottom_left']),
            'bottom_right' :
                self.origin2[1] - np.linalg.norm(_canvas_bottom_corner) * np.cos(canvas_angle['bottom_right']),
            'top_left' :
                self.origin2[1] - np.linalg.norm(_canvas_top_corner) * np.cos(canvas_angle['top_left']),
            'top_right' :
                self.origin2[1] - np.linalg.norm(_canvas_top_corner) * np.cos(canvas_angle['top_right']),
        }

        # MPL SETUP

        #fig = plt.figure(figsize=(10, 10))
        fig, axd = plt.subplot_mosaic([
            ['upper left', 'right'],
            ['lower left', 'right']
        ], gridspec_kw=dict(
            width_ratios=[1, 3], height_ratios=[2, 1]
        ), figsize=(10, 10), layout='tight')

        ax = axd['upper left']
        # TODO improve limits
        axlim = np.fmax(
            self.origin1[0] + self.pendula_distance + np.fmax(*self.params.ls),
            self.origin1[1] + np.fmax(*self.params.ls)
        ) * 1.25
        ax.set(
            autoscale_on=False, aspect='equal',
            xlim=(-axlim, axlim), ylim=(-axlim, axlim)
        )
        ax_canvas = axd['right']
        ax_canvas.set(
            autoscale_on=False, aspect='equal',
            xlim=(-self.canvas_sidelength, self.canvas_sidelength),
            ylim=(-self.canvas_sidelength, self.canvas_sidelength),
        )

        ax_angle = axd['lower left']
        ax_angle.set(
            autoscale_on=False,
            xlim=(0, self.stop_time),
            ylim=(-2 * np.pi, 2 * np.pi),
        )

        for k in axd:
            if not self.debug:
                axd[k].axis('off')
            else:
                axd[k].grid(color='darkgray', ls=':')
                axd[k].tick_params(
                    # affect x-axis and y-axis
                    axis='both',
                    # affect major and minor ticks
                    which='both',
                    # disable ticks
                    top=False, left=False, bottom=False, right=False,
                    # disable labels
                    labeltop=False, labelleft=False, labelbottom=False, labelright=False,
                )

        _, = ax_angle.plot(self.time, angle1)
        color1 = mpl.artist.getp(_, 'color')
        _, = ax_angle.plot(self.time, angle2)
        color2 = mpl.artist.getp(_, 'color')

        angle1_trace, = ax_angle.plot([], [], '.', color=color1)
        angle2_trace, = ax_angle.plot([], [], '.', color=color2)

        pendulum1_arm, = ax.plot([], [], 'o-', color=color1, lw=2)
        pendulum2_arm, = ax.plot([], [], 'o-', color=color2, lw=2)

        pen_arm, = ax.plot([], [], '-', color=color1, lw=2)

        harmonogram = mpl.collections.LineCollection(
            [],
            cmap=plt.get_cmap('YlOrRd'), alpha=0.6,
            norm=plt.Normalize(0, self.time.shape[-1]),
            lw=1, ls=':',
        )
        ax_canvas.add_collection(harmonogram)

        canvas_bottom, = ax.plot([], [], '-', color=color2, lw=2)
        canvas_top,    = ax.plot([], [], '-', color=color2, lw=2)
        canvas_left,   = ax.plot([], [], '-', color=color2, lw=2)
        canvas_right,  = ax.plot([], [], '-', color=color2, lw=2)

        time_template = 'time = %.1fs'
        time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)

        harmonogram_history_x = deque(maxlen=trace_length)
        harmonogram_history_y = deque(maxlen=trace_length)

        def animate(i):

            # RESET

            if i == 0:
                harmonogram_history_x.clear()
                harmonogram_history_y.clear()
                harmonogram.set_segments([])

            # CALCULATE

            angle_x = [i * self.dt]
            angle1_y = [angle1[i]]
            angle2_y = [angle2[i]]

            pendulum1_arm_curr_x = [self.origin1[0], bob1_pos_x[i]]
            pendulum1_arm_curr_y = [self.origin1[1], bob1_pos_y[i]]

            pendulum2_arm_curr_x = [self.origin2[0], bob2_pos_x[i]]
            pendulum2_arm_curr_y = [self.origin2[1], bob2_pos_y[i]]

            projected_pen_length = self.pen_length * np.cos(np.arctan2(*self.pen_properties))
            gamma = projected_pen_length / self.params.ls[0]

            pen_arm_curr_x = [self.origin1[0] + gamma * np.diff(pendulum1_arm_curr_x)[0], pen_pos_x[i]]
            pen_arm_curr_y = [self.origin1[1] + gamma * np.diff(pendulum1_arm_curr_y)[0], pen_pos_y[i]]

            canvas_bottom_curr_x = [
                canvas_pos_x['bottom_left'][i],
                canvas_pos_x['bottom_right'][i]
            ]
            canvas_top_curr_x = [
                canvas_pos_x['top_left'][i],
                canvas_pos_x['top_right'][i]
            ]
            canvas_left_curr_x = [
                canvas_pos_x['bottom_left'][i],
                canvas_pos_x['top_left'][i]
            ]
            canvas_right_curr_x = [
                canvas_pos_x['bottom_right'][i],
                canvas_pos_x['top_right'][i]
            ]

            canvas_bottom_curr_y = [
                canvas_pos_y['bottom_left'][i],
                canvas_pos_y['bottom_right'][i]
            ]
            canvas_top_curr_y = [
                canvas_pos_y['top_left'][i],
                canvas_pos_y['top_right'][i]
            ]
            canvas_left_curr_y = [
                canvas_pos_y['bottom_left'][i],
                canvas_pos_y['top_left'][i]
            ]
            canvas_right_curr_y = [
                canvas_pos_y['bottom_right'][i],
                canvas_pos_y['top_right'][i]
            ]

            # UPDATE

            angle1_trace.set_data(angle_x, angle1_y)
            angle2_trace.set_data(angle_x, angle2_y)

            pendulum1_arm.set_data(pendulum1_arm_curr_x, pendulum1_arm_curr_y)
            pendulum2_arm.set_data(pendulum2_arm_curr_x, pendulum2_arm_curr_y)

            pen_arm.set_data(pen_arm_curr_x, pen_arm_curr_y)

            canvas_bottom.set_data(canvas_bottom_curr_x, canvas_bottom_curr_y)
            canvas_top.set_data(canvas_top_curr_x, canvas_top_curr_y)
            canvas_left.set_data(canvas_left_curr_x, canvas_left_curr_y)
            canvas_right.set_data(canvas_right_curr_x, canvas_right_curr_y)

            time_text.set_text(time_template % (i * self.dt))

            harmonogram_x = pen_pos_x[i] - canvas_pos_x['center'][i]
            harmonogram_y = pen_pos_y[i] - canvas_pos_y['center'][i]

            if np.abs(harmonogram_x) <= self.canvas_sidelength:
                harmonogram_history_x.appendleft(harmonogram_x)
            if np.abs(harmonogram_y) <= self.canvas_sidelength:
                harmonogram_history_y.appendleft(harmonogram_y)

            harmonogram.set_segments(history2segments(
                harmonogram_history_x, harmonogram_history_y
            ))
            harmonogram.set_array(
                list(range(i, np.max(i-trace_length, 0), -1))
            )

            # DRAW

            artists = [
                angle1_trace, angle2_trace,
                harmonogram,
                pendulum1_arm, pendulum2_arm,
                pen_arm,
                canvas_bottom, canvas_top, canvas_left, canvas_right,
                time_text
            ]

            return artists

        _ = animation.FuncAnimation(
            fig, animate, self.solution.shape[-1], interval=1, blit=True
        )
        plt.show()

if __name__ == '__main__':
    parser = ap.ArgumentParser(
        prog='harmonograph',
        description='Generating harmonic, yet chaotic images using damped coupled pendula',
    )
    parser.add_argument(
        '-l1', '--length_pendulum1',
        help='length of the left pendulum\'s arm. SI unit: [m]',
        action='store',
        default=3.0,
    )
    parser.add_argument(
        '-l2', '--length_pendulum2',
        help='length of the right pendulum\'s arm. SI unit: [m]',
        action='store',
        default=3.0,
    )
    parser.add_argument(
        '-m1', '--mass_pendulum1',
        help='mass of the left pendulum\'s bob. SI unit: [kg]',
        action='store',
        default=200.0,
    )
    parser.add_argument(
        '-m2', '--mass_pendulum2',
        help='mass of the right pendulum\'s bob. SI unit: [kg]',
        action='store',
        default=400.0,
    )
    parser.add_argument(
        '-s', '--spring_stiffness',
        help='stiffness of the spring connecting the pendula. SI unit: [kg/s^2]',
        action='store',
        default=1.0,
    )
    parser.add_argument(
        '-d', '--damping_coefficient',
        help='stiffness of the spring connecting the pendula. SI unit: [kg/s]',
        action='store',
        default=1.0,
    )
    parser.add_argument(
        '-g', '--gravity',
        help='gravitational acceleration. SI unit: [m/s^2]',
        action='store',
        default=9.81,
    )
    parser.add_argument(
        '-t', '--time',
        help='how many seconds to simulate the system. SI unit: [s]',
        action='store',
        default=2000,
    )
    parser.add_argument(
        '-c', '--cmap',
        help='defines which colormap the pen uses for drawing (c.f. https://matplotlib.org/stable/users/explain/colors/colormaps.html)',
        choices=[
            # Perceptually Uniform Sequential
            'viridis', 'plasma', 'inferno', 'magma', 'cividis',
            # Sequential
            'Greys', 'Purples', 'Blues', 'Greens', 'Oranges', 'Reds',
            'YlOrBr', 'YlOrRd', 'OrRd', 'PuRd', 'RdPu', 'BuPu',
            'GnBu', 'PuBu', 'YlGnBu', 'PuBuGn', 'BuGn', 'YlGn',
            # Sequential 2
            'binary', 'gist_yarg', 'gist_gray', 'gray', 'bone',
            'pink', 'spring', 'summer', 'autumn', 'winter', 'cool',
            'Wistia', 'hot', 'afmhot', 'gist_heat', 'copper',
            # Diverging
            'PiYG', 'PRGn', 'BrBG', 'PuOr', 'RdGy', 'RdBu', 'RdYlBu',
            'RdYlGn', 'Spectral', 'coolwarm', 'bwr', 'seismic',
            # Cyclic
            'twilight', 'twilight_shifted', 'hsv',
        ]
    )
    parser.add_argument(
        '-a', '--alpha',
        help='defines the alpha/transparency of the pen stroke',
        type=float_range,
    )
    parser.add_argument(
        '--dark',
        help='whether or not to use dark background',
        action='store_true',
    )

    args = parser.parse_args()

    if args.dark:
        plt.style.use(['dark_background'])

    harmonograph = Harmonograph(
        lengths=(
            3.0, # pendulum1 : length of arm
            3.0  # pendulum2 : length of arm
        ),
        masses=(
            20.0, # pendulum1 : mass of bob
            30.0  # pendulum2 : mass of bob
        ),
        spring=10.0, # spring stiffness in [kg/s]
        damping=1.0, # (under/critical/over) damping (</=/>) 2 * np.sqrt(mass * spring)
        gravity=9.81,
        stop_time=200, # number of seconds to simulate system
        initial_conditions=(
            # pendulum1 : initial angle
            np.deg2rad(30),
            # pendulum1 : initial angular velocity
            np.deg2rad(-8),
            # pendulum2 : initial angle
            np.deg2rad(-10),
            # pendulum2 : initial angular velocity
            np.deg2rad(-12),
        ),
        # TODO: internalize into Harmonograph()
        origin=(0, 0), # origin of pendulum1
        pendula_distance = 1.0,
        pen_properties=(
            1.15, # pen arm length
            2.0, # pen offset
        ),
        canvas_properties=(
            2.0, # canvas offset
            2.0, # canvas side length
        ),
        debug=False
    )

    #harmonograph.animate_full(trace_length=-120)
    #harmonograph.show_angle()
    harmonograph.show_harmonogram()

