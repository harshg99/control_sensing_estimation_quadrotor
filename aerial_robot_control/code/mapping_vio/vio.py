#%% Imports

import numpy as np
from numpy.linalg import inv
from numpy.linalg import norm
from scipy.spatial.transform import Rotation

from pdb import set_trace as T
#%% Functions

def nominal_state_update(nominal_state, w_m, a_m, dt):
    """
    function to perform the nominal state update

    :param nominal_state: State tuple (p, v, q, a_b, w_b, g)
                    all elements are 3x1 vectors except for q which is a Rotation object
    :param w_m: 3x1 vector - measured angular velocity in radians per second
    :param a_m: 3x1 vector - measured linear acceleration in meters per second squared
    :param dt: duration of time interval since last update in seconds
    :return: new tuple containing the updated state
    """
    # Unpack nominal_state tuple
    p, v, q, a_b, w_b, g = nominal_state

    # YOUR CODE HERE
    #T()
    new_p = p + v*dt + 0.5*(q.as_matrix()@(a_m-a_b)+g)*dt*dt
    new_v = v + (q.as_matrix()@(a_m-a_b)+g)*dt
    new_q = Rotation.from_matrix(q.as_matrix()@Rotation.from_rotvec((w_m-w_b).flatten()*dt).as_matrix())

    return new_p, new_v, new_q, a_b, w_b, g

def skew(vec):
    return np.array([[0 ,  -vec[2]  ,  vec[1]],
        [  vec[2] ,    0  , -vec[0]],
        [ -vec[1] , vec[0],   0]])

def error_covariance_update(nominal_state, error_state_covariance, w_m, a_m, dt,
                            accelerometer_noise_density, gyroscope_noise_density,
                            accelerometer_random_walk, gyroscope_random_walk):
    """
    Function to update the error state covariance matrix

    :param nominal_state: State tuple (p, v, q, a_b, w_b, g)
                        all elements are 3x1 vectors except for q which is a Rotation object
    :param error_state_covariance: 18x18 initial error state covariance matrix
    :param w_m: 3x1 vector - measured angular velocity in radians per second
    :param a_m: 3x1 vector - measured linear acceleration in meters per second squared
    :param dt: duration of time interval since last update in seconds
    :param accelerometer_noise_density: standard deviation of accelerometer noise
    :param gyroscope_noise_density: standard deviation of gyro noise
    :param accelerometer_random_walk: accelerometer random walk rate
    :param gyroscope_random_walk: gyro random walk rate
    :return:
    """

    # Unpack nominal_state tuple
    p, v, q, a_b, w_b, g = nominal_state
    #T()
    # YOUR CODE HERE
    R = q.as_matrix()
    Fx = np.identity(18)
    #error_state_covariance = np.diag([0, 0, 0, 0.15, 0.15, 0.15, 0.008, 0.008, 0.008, 0.06, 0.06, 0.06, 0.001, 0.001, 0.001, 0, 0, 0])
    Fx[0:3,3:6]   = np.identity(3)*dt
    Fx[3:6,6:9]   = -R@skew((a_m-a_b))*dt
    Fx[3:6,9:12]  = -R*dt 
    Fx[3:6,15:18] = np.identity(3)*dt 
    Fx[6:9,6:9]  = Rotation.from_rotvec((w_m-w_b).flatten()*dt).as_matrix().T
    Fx[6:9,12:15] = -np.identity(3)*dt

    Fi = np.zeros((18,12))
    Fi[3:15,0:12] = np.identity(12)
    Qi =  np.identity(12)
    Qi[6:9,6:9] = np.power(accelerometer_random_walk,2)*dt*dt*np.identity(3)
    Qi[9:12,9:12] = np.power(gyroscope_random_walk,2)*dt*dt*np.identity(3)
    Qi[0:3,0:3] = np.power(accelerometer_noise_density,2)*np.identity(3)
    Qi[3:6,3:6] = np.power(gyroscope_noise_density,2)*np.identity(3)

    error_cov_new = Fx@error_state_covariance@Fx.T + Fi@Qi@Fi.T
    # return an 18x18 covariance matrix
    return error_cov_new


def measurement_update_step(nominal_state, error_state_covariance, uv, Pw, error_threshold, Q):
    """
    Function to update the nominal state and the error state covariance matrix based on a single
    observed image measurement uv, which is a projection of Pw.

    :param nominal_state: State tuple (p, v, q, a_b, w_b, g)
                        all elements are 3x1 vectors except for q which is a Rotation object
    :param error_state_covariance: 18x18 initial error state covariance matrix
    :param uv: 2x1 vector of image measurements
    :param Pw: 3x1 vector world coordinate
    :param error_threshold: inlier threshold
    :param Q: 2x2 image covariance matrix
    :return: new_state_tuple, new error state covariance matrix
    """
    
    # Unpack nominal_state tuple
    p, v, q, a_b, w_b, g = nominal_state
    #T()
    # YOUR CODE HERE - compute the innovation next state, next error_state covariance
    R = q.as_matrix()
    Pc = R.T@(Pw-p)
    dzdp = np.array([[1, 0, -Pc[0]/Pc[2]],[0, 1, -Pc[1]/Pc[2]]])/Pc[2]
    dpdthe = skew(Pc)
    dpdP = -R.T
    H = np.zeros((2,18))
    H[:,0:3] = dzdp@dpdP
    H[:,6:9] = dzdp@dpdthe

    Kt = error_state_covariance@H.T@np.linalg.inv(H@error_state_covariance@H.T + Q)
    factor = np.identity(18) - Kt@H
    error_state_covariance = factor@error_state_covariance@factor.T + Kt@Q@Kt.T
    
    x = np.stack([p.flatten(),v.flatten(),q.as_rotvec(),a_b.flatten(),w_b.flatten(),g.flatten()]).reshape(-1)
    innovation = uv - np.array([Pc[0]/Pc[2],Pc[1]/Pc[2]])
    deltax = Kt@innovation
    if np.linalg.norm(innovation)<error_threshold:
        p +=deltax[0:3]
        v +=deltax[3:6]
        q = Rotation.from_matrix((R@Rotation.from_rotvec(deltax[6:9].flatten()).as_matrix()))
        a_b += deltax[9:12]
        w_b += deltax[12:15]
        g += deltax[15:18]

    return (p, v, q, a_b, w_b, g), error_state_covariance, innovation
