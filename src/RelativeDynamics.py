import numpy as np

class RelativeDynamics:
    def __init__(self, mc, Jc, initial_state):
        """
        Initializes the RelativeDynamics class.

        Args:
            mc (float): Mass of the chaser spacecraft.
            Jc (np.ndarray): Inertia matrix of the chaser spacecraft (3x3).
        """
        self.mc = mc
        self.Jc = Jc
        self.mu = 3.986004418e14  # Standard gravitational parameter for Earth (m^3/s^2)
        self.sigma  = initial_state.sigma   # Modified Rodrigues Parameters describing relative attitude between the chaser and the target
        self.omega  = initial_state.omega   # The relative angular velocity between the chaser and the target, expressed in the chaser frame
        self.rho    = initial_state.rho     # The relative position between the chaser and the target, expressed in the chaser frame
        self.v      = initial_state.v       # The relative velocity between the chaser and the target, expressed in the chaser frame

    def ode_model(self, input_force, input_torque, target_state):
        """
        Computes the time derivative of the state vector.
        
        Args:
            input_force (np.ndarray): Control force applied to the chaser (3D vector).
            input_torque (np.ndarray): Control torque applied to the chaser (3D vector).
            target_state (object): An object containing the target's state information, including:
                - omega_it_t (np.ndarray): Target's inertial angular velocity expressed in the target frame (3D vector).
                - dot_omega_it_t (np.ndarray): Time derivative of target's inertial angular velocity (3D vector).
                - v_t (np.ndarray): Target's inertial velocity expressed in the target frame (3D vector).
                - dot_v_t (np.ndarray): Time derivative of target's inertial velocity (3D vector).
                - r_t (np.ndarray): Target's inertial position expressed in the target frame (3D vector).
                - dot_r_t (np.ndarray): Time derivative of target's inertial position (3D vector).
                - sigma (np.ndarray): Target's modified Rodrigues parameters (3D vector).
                
        Returns:
            np.ndarray: Time derivative of the state vector (12D vector).
        """
        G_mrp = self.G_mrp()
        sigma_dot = G_mrp @ self.omega
        omega_dot = np.linalg.inv(self.Jc) @ (input_torque - self.compute_C1_matrix(self.omega, self.sigma, target_state.omega_it_t) @ self.omega - self.compute_D1_vector(self.sigma, target_state.omega_it_t, target_state.dot_omega_it_t))
        rho_dot = self.v - self.compute_C2_matrix(self.omega, self.sigma, target_state.omega_it_t) @ self.rho
        v_dot = (1/self.mc) * input_force + self.compute_D2_vector(self.rho, self.sigma, target_state.dot_v_t)
        x_dot = np.vstack((sigma_dot, omega_dot, rho_dot, v_dot))
        return x_dot
        
        
    def skew_symmetric(self, vec):
        """
        Generates a skew-symmetric matrix from a 3D vector.

        Args:
            vec (np.ndarray): A 3D vector.

        Returns:
            np.ndarray: The skew-symmetric matrix.
        """
        return np.array([
            [0, -vec[2], vec[1]],
            [vec[2], 0, -vec[0]],
            [-vec[1], vec[0], 0]
        ])

    def G_mrp(self):
        """
        Calculates the G matrix for modified Rodrigues parameters (MRPs) kinematics.

        Args:
            sigma (np.ndarray): MRPs (3D vector).

        Returns:
            np.ndarray: The G matrix (3x3).
        """
        sigma_norm_sq = np.dot(self.sigma, self.sigma)
        I3 = np.eye(3)
        Omega_sigma = self.skew_symmetric(self.sigma)
        return 0.25 * ((1 - sigma_norm_sq) * I3 + 2 * Omega_sigma + 2 * np.outer(self.sigma, self.sigma))

    def R_t_c(self, sigma):
        """
        Calculates the rotation matrix from target frame to chaser frame using MRPs.

        Args:
            sigma (np.ndarray): Modified Rodrigues Parameters (3D vector).

        Returns:
            np.ndarray: The rotation matrix R_t^c (3x3).
        """
        sigma_norm_sq = np.dot(sigma, sigma)
        I3 = np.eye(3)
        Omega_sigma = self.skew_symmetric(sigma)
        den = (1 + sigma_norm_sq)**2
        R = I3 - (4 * (1 - sigma_norm_sq) / den) * Omega_sigma + (8 / den) * (Omega_sigma @ Omega_sigma)
        return R

    def compute_C1_matrix(self, omega, sigma, omega_it_t):
        """
        Calculates the C1 matrix for attitude dynamics.
        C1(omega) = -Jc Omega(R_t^c omega_it_t) - Omega(R_t^c omega_it_t) Jc + Omega(Jc(omega + R_t^c omega_it_t))
        """
        R_tc = self.R_t_c(sigma)
        omega_it_t_c = R_tc @ omega_it_t
        Omega_omega_it_t_c = self.skew_symmetric(omega_it_t_c)
        
        C1 = -self.Jc @ Omega_omega_it_t_c - Omega_omega_it_t_c @ self.Jc + self.skew_symmetric(self.Jc @ (omega + omega_it_t_c))
        return C1

    def compute_D1_vector(self, sigma, omega_it_t, dot_omega_it_t):
        """
        Calculates the D1 vector for attitude dynamics.
        D1(omega) = -Omega(R_t^c omega_it_t) Jc R_t^c omega_it_t - Jc R_t^c dot_omega_it_t
        """
        R_tc = self.R_t_c(sigma)
        omega_it_t_c = R_tc @ omega_it_t
        dot_omega_it_t_c = R_tc @ dot_omega_it_t
        Omega_omega_it_t_c = self.skew_symmetric(omega_it_t_c)
        
        D1 = -Omega_omega_it_t_c @ self.Jc @ omega_it_t_c - self.Jc @ dot_omega_it_t_c
        return D1

    def compute_C2_matrix(self, omega, sigma, omega_it_t):
        """
        Calculates the C2 matrix for translational dynamics.
        C2 = -Omega(omega_ic^c)
        omega_ic^c = omega + R_t^c omega_it^t
        """
        R_tc = self.R_t_c(sigma)
        omega_it_t_c = R_tc @ omega_it_t
        omega_ic_c = omega + omega_it_t_c
        C2 = -self.skew_symmetric(omega_ic_c)
        return C2

    def compute_D2_vector(self, r_c, sigma, dot_v_t):
        """
        Calculates the D2 vector for translational dynamics.
        D2 = -mu/||r_c||^3 * r_c - R_t^c dot_v_t
        """
        R_tc = self.R_t_c(sigma)
        D2 = -self.mu / (np.linalg.norm(r_c)**3) * r_c - R_tc @ dot_v_t
        return D2
