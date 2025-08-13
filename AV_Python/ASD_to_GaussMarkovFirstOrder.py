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
        self.phi = B_c[1, 1]              # From eqn. (4.115) of [2]
        self.Qb  = self.phi * B_c[0, 1]   # Process noise covariance, from eqn. (4.116) of [3]
        self.Qn = self.Sn * self.Fs       # Discrete time random walk process noise. See (62) in [1]
        self.Pb_ss_d = self.Qb / (1 - self.phi**2)  # Steady-state covariance of the bias process

    def get_Pb_ss(self):
        if (self.Pb_ss_c - self.Pb_ss_d)/self.Pb_ss_c > 1e-3:
            print(f"Warning: Continuous-time steady-state covariance Pb_ss_c ({self.Pb_ss_c:.4e}.) "
                  f"differs from discrete-time Pb_ss_d ({self.Pb_ss_d:.4e}).")   
        return self.Pb_ss_d

    def print_discrete_time_model(self):
        """
        Print the discrete-time model of the Gauss-Markov First Order process.
        """
        print(f"Discrete-time model is (units depend on context): ")
        print(f"\t   z[k] = b[k] + n[k],")
        print(f"\t b[k+1] = {self.phi:.4f} * b[k] + w[k],")
        print(f"where n[k] is white noise with covariance Qn = {self.Qn:.4e} (std {np.sqrt(self.Qn):.4e})")
        print(f"and w[k] is the bias process with covariance Qw = {self.Qb:.4e} (std {np.sqrt(self.Qb):.4e}).")
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
            b[i] = self.phi * b[i-1] + w[i]  # Update bias

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
    
    def robust_standard_deviation(self, data):
        """
        Compute the Median Absolute Deviation (MAD) of the data.
        :param data: Input data for which the MAD is to be computed.
        :return: The computed MAD value.
        """
        median = np.median(data)
        mad = np.median(np.abs(data - median))
        std = mad * 1.4826  # Scale factor for normal distribution
        return std 

    def test_cov(self, data, T_avg = 1.0, T_int = 1.0, T_shft = 1.0, MaxRep = 1000):
        """
        Plot sample error trajectories against the model prediction of their standard deviation.
               average                 integrate 
        --|---------------------|----------------------|--------
         t0                   t0+T_avg            t0+T_avg+T_int  
        t1 = t0 + T_shft
                           average                 integrate 
        ----------|---------------------|----------------------|--------
                 t1                   t1+T_avg            t1+T_avg+T_int  
        :param data: Input data for which the covariance is to be computed.
        :param T_avg: Duration over which the data is averaged to estimate the bias.
        :param T_int: Duration over which the data and covariance are integrated.
        :param T_shft: Interval over which the data is shifted before repeating the process
        """  
        if T_shft > (T_avg + T_int):  # min value for no overlap
            T_shft = T_avg + T_int
            print(f"Warning: T_shft ({T_shft}) changed to T_avg + T_int = {T_avg + T_int}.")
        N_avg_smpls = int(self.Fs * T_avg)  # Number of samples for averaging
        N_int_smpls = int(self.Fs * T_int)  # Number of samples for integration
        N_shft_smpls = int(self.Fs * T_shft)  # Number of samples for shifting

        print(f"T_avg = {T_avg} seconds, T_int = {T_int} seconds, T_shft = {T_shft} seconds.")
        print(f"Number of samples for averaging: {N_avg_smpls}, "
              f"Number of samples for integration: {N_int_smpls}, "
              f"Number of samples for shifting: {N_shft_smpls}.")
        print(f"Fs = {self.Fs} Hz, Ts = {self.Ts} seconds.")

        N_data = len(data)  # Total number of samples in the data
        N_rep = min(int(N_data / (N_avg_smpls+N_int_smpls) ), MaxRep)
        # simulate the covariance
        P_theory = np.zeros(N_int_smpls)  # Preallocate the covariance array
        P_b = np.zeros(N_int_smpls)  # Preallocate the bias covariance array
        P_z = np.zeros(N_int_smpls)  # Preallocate the measurement covariance array
        self.print_discrete_time_model()
        
        # setup the discrete-time state-space model 
        # with state vector x = [phase, bias]
        # x(k+1) = Phi * x(k) + Gam * w(k)
        Phi = np.zeros((2, 2))
        Phi[0, 0] = 1         # Phase integration
        Phi[0, 1] = self.Ts   # Phase integration of frequency error
        Phi[1, 1] = self.phi  # Bias state transition
        Gam = np.zeros((2, 2))
        Gam[0, 0] = self.Ts  # integrate phase RW noise
        Gam[1, 1] = self.Ts  # integrate bias process noise
        Q = np.zeros((2, 2))  # Preallocate the process noise covariance
        Q[0, 0] = self.Qn
        Q[1, 1] = self.Qb
        P = np.zeros((2, 2))  # Preallocate the error state covariance
        P[1, 1] = self.get_Pb_ss()
        P_b      = np.zeros(N_int_smpls)  # Preallocate the bias covariance array
        P_theory = np.zeros(N_int_smpls)  # Preallocate the theoretical covariance array
        P_theory[0] = P[0,0]
        P_b[0] = P[1,1]
        for i in range(1, N_int_smpls):
            P = Phi @ P @ Phi.T + Gam @ Q @ Gam.T  # Update the error state covariance
            P_theory[i] = P[0, 0]
            P_b[i] = P[1, 1]
        #std_b = np.sqrt(P_b)            # Standard deviation of the frequency bias
        #std_z = np.sqrt(P_z)            # Standard deviation of the frequency error
        std_theory = np.sqrt(P_theory)  # Standard deviation of the phase error

        # compute sample error trajectories
        bias = np.zeros(N_rep)  # Preallocate the bias array
        phserr = np.zeros((N_rep, N_int_smpls))  # Preallocate the error array
        for i in range(N_rep):
            N_start = i * N_shft_smpls  # Starting index for the current segment
            N_middle = N_start + N_avg_smpls  # Middle index for the current segment
            N_end = N_middle + N_int_smpls  # End index for the current segment
            # Average the data to estimate the bias 
            idx4bias    = range(N_start, N_middle)
            bias[i]     = np.mean(data[idx4bias])  # Average the data to estimate the bias
            idx4int     = range(N_middle, N_middle + N_int_smpls)
            phserr[i,:] = np.cumsum(data[idx4int] - bias[i]) * self.Ts  # Integrate the phase error trajectory

        phsErrStd = np.zeros(N_int_smpls)  # Preallocate the phase error standard deviation array
        # Calculate the standard deviation of the phase error at each sample
        for j in range(0,N_int_smpls):
            phsErrStd[j] = self.robust_standard_deviation(phserr[:,j])  # Standard deviation of the phase error at each sample

        # Plot the results
        tm = np.arange(0, N_int_smpls) * self.Ts
        plt.figure(figsize=(10, 6))
        for i in range(N_rep):
            plt.plot(tm, phserr[i,:],  linestyle='--', linewidth=0.5, alpha = 0.5)
        #plt.plot(tm,  std_b,  linestyle='--', color='black', linewidth=1)
        #plt.plot(tm, -std_b,  linestyle='--', color='black', linewidth=1)
        plt.plot(tm,  std_theory, linestyle='-', color='black', linewidth=1)
        plt.plot(tm, -std_theory, linestyle='-', color='black', linewidth=1)
        plt.plot(tm,  phsErrStd, linestyle='-', color='blue', linewidth=1)
        plt.plot(tm, -phsErrStd, linestyle='-', color='blue', linewidth=1)
        plt.title('Sample Error Trajectories and Standard Deviation')
        plt.xlabel('Time (samples)')
        plt.ylabel('Error Trajectories and Standard Deviation')
        plt.ylim([-5*std_theory.max(), 5*std_theory.max()])
        plt.grid()
        plt.show()

        #plt.figure(figsize=(10, 6))
        #plt.scatter(range(N_rep), bias, label='Bias Estimate', color='orange')
        #plt.title('Bias Estimate from Sample Data')
        #plt.xlabel('Repetition Index')
        #plt.ylabel('Bias Estimate')
        #plt.grid()
        #plt.show()