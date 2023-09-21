---
title: 'Harmonograph'
subtitle: 'Producing Harmonic Yet Chaotic Images'
author: 'Thierry Meier'
date: '2nd of September, 2023'
header-includes:
  - \usepackage{mathtools}
  - \usepackage{caption}
---

# Experiment Setup

\begin{figure}[h!]
    \centering
    \includegraphics[width=0.75\textwidth]{harmonograph.jpg}
    %\caption{goes here}
    \label{fig:experiment_setup}
\end{figure}

Here we have a depiction of the Harmonograph as found in the Swiss Science Center [Technorama](www.technorama.ch), located in the municipality of Winterthur, Switzerland. It consists of two identical pendula coupled via a spring with one arm holding a pen whilst the other holds the canvas. This contraption then produces chaotic, yet harmonic images known as Harmonograms by dislocating each of the pendula.
Since the nature of this system is chaotic, tiny changes in the initial conditions produce vastly different images as can be seen in the background.

\pagebreak

# Derivation of Equation of Motion of Damped Coupled Pendula

The acting forces on pendulum $(i)$, where $i \in \{0,1\}$, are as follows:

\begin{align}
    F_g^{(i)} &= -g m_i \sin(\theta_i)                    \quad &&\text{(gravitational force)}\\
    F_s^{(i)} &= k (x_{i+1 \bmod 2} - x_i) \sin(\theta_i) \quad &&\text{(spring force)}\\
    F_d^{(i)} &= -d v_i                                   \quad &&\text{(damping force)}
\end{align}

For notational convenience let $j \coloneqq i + 1 \bmod 2$.

Thus, the net force acting on pendulum $(i)$ is

\begin{align*}
F_\text{net}^{(i)} &= (k (x_j - x_i) - g m_i) \sin(\theta_i) - d v_i
\end{align*}

By Newton's second law of motion, we have $F=ma$, thus

\begin{align*}
F_\text{net}^{(i)} = m_i a_i &\stackrel{!}{=} (k (x_j - x_i) - g m_i) \sin(\theta_i) - d v_i\\
                   \iff a_i &= (\frac{k}{m_i} (x_j - x_i) - g) \sin(\theta_i) - \frac{d}{m_i} v_i\\
\end{align*}

Because of the relationship between position, velocity and acceleration, i.e. $v = \dot{x}$ and $a = \ddot{x}$, we have

\begin{align}
\label{eq:sys}
\ddot{x} &= (\frac{k}{m_i} (x_j - x_i) - g) \sin(\theta_i) - \frac{d}{m_i} \dot{x}
\end{align}

By further inspection of the experiment setup we observe that

\begin{align}
\label{eq:pos}
x_i &= l_i \sin(\theta_i)\\
\label{eq:vel}
\implies \dot{x}_i &= l_i \dot{\theta}_i \cos{\theta_i}\\
\label{eq:acc}
\implies \ddot{x}_i &= l_i \ddot{\theta}_i \cos{\theta_i} - l_i \dot{\theta}_i^2 \sin{\theta_i}
\end{align}

\pagebreak

By substituting equation (\ref{eq:pos}), (\ref{eq:vel}) and (\ref{eq:acc}) into (\ref{eq:sys}), we get

\begin{align*}
l_i \ddot{\theta}_i \cos{\theta_i} - l_i \dot{\theta}_i^2 \sin{\theta_i} &= 
\Big (\frac{k}{m_i} (l_j \sin{\theta_j} - l_i \sin{\theta_i}) - g \Big) \sin{\theta_i} - \frac{d}{m_i} l_i \dot{\theta}_i \cos{\theta_i}
\end{align*}

We rearrange and simplify as follows, yielding

\begin{align*}
l_i \ddot{\theta}_i \cos{\theta_i} &= 
%\Big (\frac{k}{m_i} (l_j \sin{\theta_j} - l_i \sin{\theta_i}) - g \Big) \sin{\theta_i} 
%+ l_i \dot{\theta}_i^2 \sin{\theta_i} 
%- \frac{d}{m_i} l_i \dot{\theta}_i \cos{\theta_i}\\
%&=
\Big (l_i \dot{\theta}_i^2 + \frac{k}{m_i} (l_j \sin{\theta_j} - l_i \sin{\theta_i}) - g \Big) \sin{\theta_i}
- \frac{d}{m_i} l_i \dot{\theta}_i \cos{\theta_i}\\
\iff \ddot{\theta}_i &=
\Big (l_i \dot{\theta}_i^2 + \frac{k}{m_i} (l_j \sin{\theta_j} - l_i \sin{\theta_i}) - g \Big) \frac{\sin{\theta_i}}{l_i \cos{\theta_i}}
%+ \frac{k (l_j \sin{\theta_j} - l_i \sin{\theta_i})}{m_i l_i \cos{\theta_i}}
- \frac{d l_i \dot{\theta}_i \cos{\theta_i}}{m_i l_i \cos{\theta_i}}\\
&=
\Big (\dot{\theta}_i^2 + \frac{k}{m_i l_i} (l_j \sin{\theta_j} - l_i \sin{\theta_i}) - \frac{g}{l_i} \Big) \frac{\sin{\theta_i}}{\cos{\theta_i}}
%\frac{m_i (l_i \dot{\theta}_i^2 - g) \sin{\theta_i} + k (l_j \sin{\theta_j} - l_i \sin{\theta_i})}
%{m_i l_i \cos{\theta_i}}
- \frac{d}{m_i} \dot{\theta}_i\\
&=
\Big (\dot{\theta}_i^2 + \frac{k}{m_i l_i} (l_j \sin{\theta_j} - l_i \sin{\theta_i}) - \frac{g}{l_i} \Big) \tan{\theta_i}
- \frac{d}{m_i} \dot{\theta}_i
\end{align*}

This results in the following coupled system of second order ordinary differential equations

\begin{equation*}
\begin{cases}
    \quad \ddot{\theta}_0 &= \quad

\displaystyle \Big (\dot{\theta}_0^2 + \frac{k}{m_0 l_0} (l_1 \sin{\theta_1} - l_0 \sin{\theta_0}) - \frac{g}{l_0} \Big) \tan{\theta_0}
- \frac{d}{m_0} \dot{\theta}_0\\[15pt]

        %\displaystyle \frac{(m_0 (l_0 \dot{\theta}_0^2 - g) - k l_0) \sin{\theta_0} + k l_1 \sin{\theta_1}}{m_0 l_0 \cos{\theta_0}}
        %- \frac{d}{m_0} \dot{\theta}_0\\[15pt]

    \quad \ddot{\theta}_1 &= \quad

\displaystyle \Big (\dot{\theta}_1^2 + \frac{k}{m_1 l_1} (l_0 \sin{\theta_0} - l_1 \sin{\theta_1}) - \frac{g}{l_1} \Big) \tan{\theta_1}
- \frac{d}{m_1} \dot{\theta}_1

        %\displaystyle \frac{(m_1 (l_1 \dot{\theta}_1^2 - g) - k l_1) \sin{\theta_1} + k l_0 \sin{\theta_0}}{m_1 l_1 \cos{\theta_1}}
        %- \frac{d}{m_1} \dot{\theta}_1
\end{cases}
\end{equation*}

\pagebreak

# Numerical integration with `scipy.integrate.odeint`

In order for our system of differential equations to be numerically solvable by the `scipy.integrate` module, we need to reduce the system from second to first order.

One can achieve this by introducing fresh variables as follows

\begin{align*}
    \gamma_{00} & \coloneqq \dot{\theta}_0\\
    \gamma_{01} & \coloneqq \theta_0\\
    \gamma_{10} & \coloneqq \dot{\theta}_1\\
    \gamma_{11} & \coloneqq \theta_1
\end{align*}

Trivially, their derivatives then are

\begin{align*}
    \dot{\gamma}_{00} &= \ddot{\theta}_0\\
    \dot{\gamma}_{01} &= \dot{\theta}_0\\
    \dot{\gamma}_{10} &= \ddot{\theta}_1\\
    \dot{\gamma}_{11} &= \dot{\theta}_1
\end{align*}

Now we can simply put everything together via substitution, since we derived the formulas for $\ddot{\theta}_i$ and by definition $\dot{\theta}_i = \gamma_{i0}$ and $\theta_i = \gamma_{i1}$ for $i \in \{0,1\}$. Thus our final system, ready for implementation, becomes

\begin{equation*}
\begin{cases}
    \quad \dot{\gamma}_{00} 
    &= 
    %\quad \displaystyle \frac{(m_0 (l_0 \gamma_{00}^2 - g) - k l_0) \sin{\gamma_{01}} + k l_1 \sin{\gamma_{11}}}{m_0 l_0 \cos{\gamma_{01}}}
    %- \frac{d}{m_0} \gamma_{00}\\[15pt]
    \quad \displaystyle \Big (\gamma_{00}^2 + \frac{k}{m_0 l_0} (l_1 \sin{\gamma_{11}} - l_0 \sin{\gamma_{01}}) - \frac{g}{l_0} \Big) \tan{\gamma_{01}}
    - \frac{d}{m_0} \gamma_{00}\\[15pt]
    %
    \quad \dot{\gamma}_{01} 
    &= 
    \quad \gamma_{00}\\[15pt]
    %
    \quad \dot{\gamma}_{10} 
    &= 
    %\quad \displaystyle \frac{(m_1 (l_1 \gamma_{10}^2 - g) - k l_1) \sin{\gamma_{11}} + k l_0 \sin{\gamma_{01}}}{m_1 l_1 \cos{\gamma_{11}}}
    %- \frac{d}{m_1} \gamma_{10}\\[15pt]
    \quad \displaystyle \Big (\gamma_{10}^2 + \frac{k}{m_1 l_1} (l_0 \sin{\gamma_{01}} - l_1 \sin{\gamma_{11}}) - \frac{g}{l_1} \Big) \tan{\gamma_{11}}
    - \frac{d}{m_1} \gamma_{10}\\[15pt]
    %
    \quad \dot{\gamma}_{11} 
    &= 
    \quad \gamma_{10}\\[10pt]
\end{cases}
\end{equation*}
