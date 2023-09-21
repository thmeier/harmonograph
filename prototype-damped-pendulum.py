import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import collections

def ode(ys, ts, B, g, l, A, omega):
    """
    # reduce second-order ODE system to one of first-order
    theta'' = -B * theta' - (g/l) * sin(theta) + A * cos(omega * t)

    # let new variables be as follows
    y1 = theta1'
    y2 = theta1

    # their derivative yields
    y1' = theta1''
    y2' = theta1'

    # perform substitution
    y1' = -B * y1 - (g/l) * sin(y2) + A * cos(omega * t)
    y2' = y1
    """
    y1, y2 = ys
    return np.array([
        -B * y1 - (g/l) * np.sin(y2) + A * np.cos(omega),
        y1
    ])

width, height = 1000, 1000
# length of pendula [m]
l = 2.0
# masses of pendula [kg]
m = 20.0
# gravity [m/s^2]
g = -9.81
# NOTE: gravity has enorm influence on pendulum overswinging

A = 0.0     # forcing factor
B = 0.6     # damping factor
omega = np.pi * 0.5 # angular frequency of driving oscillations

y0 = np.array([  # initial conditions
    0.0,         # angular acceleration
    0.1,         # angular velocity
])
t  = np.linspace(
    0.0,        # start time
    60.0,       # stop time
    6000        # time steps
)
dt = np.diff(t)[0]


history_len = 20  # how many trajectory points to display

y = sp.integrate.odeint(ode, y0, t, args=(B, g, l, A, omega))

ang_acc = y[:, 0]
#ang_vel = y[:, 1]

# x(t+1) = x(t) + v(t) * dt
# v(t+1) = v(t) + a(t) * dt

pos_x = l * np.sin(ang_acc)
pos_y = - l * np.cos(ang_acc)

#
# PLOT
#

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(
    autoscale_on=False, 
    xlim=(-l * 1.25, l * 1.25), 
    ylim=(-l * 1.25, l * 1.25)
)
ax.set_aspect('equal')
ax.grid()

line, = ax.plot([], [], 'o-', lw=2)
trace, = ax.plot([], [], '.-', lw=1, ms=2)
time_template = 'time = %.1fs'
time_text = ax.text(0.05, 0.9, '', transform=ax.transAxes)
history_x = collections.deque(maxlen=history_len) 
history_y = collections.deque(maxlen=history_len)

def animate(i):
    thisx = [0, pos_x[i]]
    thisy = [0, pos_y[i]]

    if i == 0:
        history_x.clear()
        history_y.clear()

    history_x.appendleft(thisx[1])
    history_y.appendleft(thisy[1])

    line.set_data(thisx, thisy)
    trace.set_data(history_x, history_y)
    time_text.set_text(time_template % (i*dt))
    return line, trace, time_text


ani = animation.FuncAnimation(
    fig, animate, len(y), interval=dt*1000, blit=True
)
plt.show()
