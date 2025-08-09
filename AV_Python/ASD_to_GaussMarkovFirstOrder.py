"""
This file is the companion of the article:
[1] Jay A. Farrell, Felipe O. Silva, Farzana Rahman, and J. Wendel. 
    "IMU Error Modeling Tutorial: INS state estimation and real-time sensor calibration
    IEEE Control Systems Magazine.
The article's main point-of-contact is J. Farrell (j.a.farrell@ieee.org)
A supplementary reference is:
[2] J. A. Farrell, "Aided Navigation: GPS with High Rate Sensors", McGraw-Hill, 2008.

This software is distributed for academic purposes under the MIT License

Copyright (c) 2025 JAY A FARRELL

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import expm

class ASD_to_GaussMarkovFirstOrder:
    """
    Class to convert AV (Allan Variance) to Gauss-Markov First Order error model.
    y(t) = b(t) + n(t), where n(t) is white noise with Power Spectral Density (PSD) Sn
    db(t)/dt = -mu_b * b(t) + w(t), mu_b >=0 and w(t) is white noise with Power Spectral Density (PSD) Sw
    """
    def __init__(self, N, B, Tp, Fs):
        """
        Initialize with Allan Standard Deviation (not Allan Variance) Parameters.
        :param N: AV Random walk parameter. Assumed input units are:
                  Accelerometer units: m/s^2/sqrt(Hz) = m/s/sqrt(s). 
                  Gyro units: deg/sqrt(hr)
        :param B: AV Bias instability parameter. Assumed input units are:
                  Accelerometer units: m/s^2. 
                  Gyro units: deg/s
        :param Tb: Delay in seconds at which GM AV is desired
        :param Fs: Sampling frequency, Hertz (Hz)
        """
        self.N = N  # See section of [1] containing eqns. (21-23)
        self.B = B
        self.Tp = Tp
        self.Fs = Fs     # Sampling frequency in Hz
        self.convert_ASD_to_continuous_time_FirstOrder()
        self.print_continuous_time_model()
        self.convert_continuous_to_discrete_time()
        self.print_discrete_time_model()

    def convert_ASD_to_continuous_time_FirstOrder(self):    
        # Convert to other parameters
        self.Ts = 1 / self.Fs # Sampling period in seconds
        self.Sn = self.N * self.N # Power Spectral Density of white measurement error, eqn. 22 in [1]
        self.Tb = self.Tp / 1.89 # Eqn. (34) in [1]  
        self.mu_b = 1 / self.Tb  # Gauss-Markov process rate parameter, eqn. (34) in [1]
        self.Sb = self.B * self.B / self.Tb / (0.4365**2)  # eqn. (37) in [1]
        self.Pb_ss_c = self.Sb / (2 * self.mu_b)  # Steady-state covariance of the bias process

    def print_continuous_time_model(self):
        """
        Print the continuous-time model of the Gauss-Markov First Order process.
        """
        print(f"Continuous-time model is (units depend on context): ")
        print(f"\t    z(t) = b(t) + n(t),")
        print(f"\tdb(t)/dt = -{self.mu_b:.4e} * b(t) + w(t),")
        print(f"where n(t) is white noise with PSD Sn = {self.Sn:.4e}")
        print(f"and w(t) is the bias process noise with PSD Sw = {self.Sb:.4e}")
        print(f"The steady-state covariance of the bias process is Pb_ss_c = {self.Pb_ss_c:.4e}.\n")

    def convert_continuous_to_discrete_time(self):
        """
        Convert the continuous-time Gauss-Markov First-order process to discrete-time.
        Returns the discrete-time model parameters.
        """
        # Continuous to discrete transformation using eqn. (4.114) of [2]
        A  = -self.mu_b
        Ac = np.array([[-A, self.Sb],
                       [0, A]]) * self.Ts
        B_c = expm(Ac)  
        # Discrete-time model parameters
        self.Phi = B_c[1, 1]              # From eqn. (4.115) of [2]
        self.Qb  = self.Phi * B_c[0, 1]   # Process noise covariance, from eqn. (4.116) of [3]
        self.Qn = self.Sn * self.Fs       # Discrete time random walk process noise. See (62) in [1]
        self.Pb_ss_d = self.Qb / (1 - self.Phi**2)  # Steady-state covariance of the bias process

    def print_discrete_time_model(self):
        """
        Print the discrete-time model of the Gauss-Markov First Order process.
        """
        print(f"Discrete-time model is (units depend on context): ")
        print(f"\t   z[k] = b[k] + n[k],")
        print(f"\t b[k+1] = {self.Phi:.4f} * b[k] + w[k],")
        print(f"where n[k] is white noise with covariance Qn = {self.Qn:.4e}")
        print(f"and w[k] is the bias process with PSD Qw = {self.Qb:.4e}.")
        print(f"The steady-state covariance of the bias process is Pb_ss_d = {self.Pb_ss_d:.4e}.\n")

    def simulate_discrete_time_model(self, duration=10, show_plots=True):
        """
        Simulate the discrete-time Gauss-Markov First Order process.
        """
        print("Simulating discrete-time Gauss-Markov First Order process for duration (seconds):", duration)

        # Generate time vector
        t = np.arange(0, duration, self.Ts)

        # Initialize state vectors
        b = np.zeros_like(t)  # Bias
        n = np.random.normal(0, np.sqrt(self.Qn), size=t.shape) # Measurement noise
        # Process noise
        w = np.random.normal(0, np.sqrt(self.Qb), size=t.shape) # Process noise

        # Simulate the process
        for i in range(1, len(t)):
            b[i] = self.Phi * b[i-1] + w[i]  # Update bias

        z = b + n  # Measurement    

        # Plot the results
        if show_plots:
            plt.figure(figsize=(10, 6))
            plt.scatter(t, z, label='Measurement (z[k])', color='blue')
            plt.scatter(t, b, label='Bias (b[k])', color='orange')
            plt.title('Simulated Discrete-Time Gauss-Markov First Order Process')
            plt.xlabel('Time (seconds)')
            plt.ylabel('Value')
            plt.legend()
            plt.grid()
            plt.show()

            plt.figure(figsize=(10, 6))
            plt.scatter(t, b, color='orange')
            plt.xlabel('Time (seconds)')
            plt.ylabel('Bias, b[k]')
            plt.legend()
            plt.grid()
            plt.show()

        return z, b, n 
