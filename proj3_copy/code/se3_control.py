import numpy as np
from scipy.spatial.transform import Rotation
import math
class SE3Control(object):
    """

    """
    def __init__(self, quad_params):
        """
        This is the constructor for the SE3Control object. You may instead
        initialize any parameters, control gain values, or private state here.

        For grading purposes the controller is always initialized with one input
        argument: the quadrotor's physical parameters. If you add any additional
        input arguments for testing purposes, you must provide good default
        values!

        Parameters:
            quad_params, dict with keys specified by crazyflie_params.py

        """

        # Quadrotor physical parameters.
        self.mass            = quad_params['mass'] # kg
        self.Ixx             = quad_params['Ixx']  # kg*m^2
        self.Iyy             = quad_params['Iyy']  # kg*m^2
        self.Izz             = quad_params['Izz']  # kg*m^2
        self.arm_length      = quad_params['arm_length'] # meters
        self.rotor_speed_min = quad_params['rotor_speed_min'] # rad/s
        self.rotor_speed_max = quad_params['rotor_speed_max'] # rad/s
        self.k_thrust        = quad_params['k_thrust'] # N/(rad/s)**2
        self.k_drag          = quad_params['k_drag']   # Nm/(rad/s)**2

        # You may define any additional constants you like including control gains.
        self.inertia = np.diag(np.array([self.Ixx, self.Iyy, self.Izz])) # kg*m^2
        self.g = 9.81 # m/s^2

        # STUDENT CODE HERE
        # Gains 

        #Linear Controller Gains
        self.position_damp = 0.75
        self.k_p1 =  5.0
        self.k_d1 = 2*self.position_damp*np.sqrt(self.k_p1)
        
        self.k_p2 =  5.0
        self.k_d2 = 2*self.position_damp*np.sqrt(self.k_p2)
        
        self.k_p3 =  5.0
        self.k_d3 = 2*self.position_damp*np.sqrt(self.k_p3)


        self.attitude_damp = 0.75
        self.k_pphi =  350.0
        self.k_dphi = 2*self.attitude_damp*np.sqrt(self.k_pphi)
        
        self.k_ptheta =  400.0
        self.k_dtheta = 2*self.attitude_damp*np.sqrt(self.k_ptheta)
        
        self.k_ppsi =  150.0
        self.k_dpsi = 2*self.attitude_damp*np.sqrt(self.k_ppsi)
        self.gamma = self.k_drag/self.k_thrust

        # Non Linear Controller Gains
        # self.position_damp = 0.80
        # self.nlk_p1 = 8.0
        # self.nlk_d1 = 2*self.position_damp*np.sqrt(self.nlk_p1)
    
        # self.nlk_p2 =  8.0
        # self.nlk_d2 = 2*self.position_damp*np.sqrt(self.nlk_p2)
        
        # self.nlk_p3 =  10.0
        # self.nlk_d3 = 2*self.position_damp*np.sqrt(self.nlk_p3)


        # self.attitude_damp = 0.75
        # self.nlk_pphi =  190.0
        # self.nlk_dphi = 2*self.attitude_damp*np.sqrt(self.nlk_pphi)
        
        # self.nlk_ptheta =  190.0
        # self.nlk_dtheta = 2*self.attitude_damp*np.sqrt(self.nlk_ptheta)
        
        # self.nlk_ppsi =  190.0
        # self.nlk_dpsi = 2*self.attitude_damp*np.sqrt(self.nlk_ppsi)

        # Non Linear Controller Gains
        self.position_damp = 0.85
        self.nlk_p1 = 6.0
        self.nlk_d1 = 2*self.position_damp*np.sqrt(self.nlk_p1)
    
        self.nlk_p2 =  6.0
        self.nlk_d2 = 2*self.position_damp*np.sqrt(self.nlk_p2)
        
        self.position_damp = 0.80
        self.nlk_p3 =  7.5
        self.nlk_d3 = 2*self.position_damp*np.sqrt(self.nlk_p3)


        self.attitude_damp = 0.65
        self.nlk_pphi =  350
        self.nlk_dphi = 2*self.attitude_damp*np.sqrt(self.nlk_pphi)
        
        self.nlk_ptheta =  350
        self.nlk_dtheta = 2*self.attitude_damp*np.sqrt(self.nlk_ptheta)
        
        self.attitude_damp = 0.65
        self.nlk_ppsi =  150
        self.nlk_dpsi = 2*self.attitude_damp*np.sqrt(self.nlk_ppsi)

        self.gamma = self.k_drag/self.k_thrust

        self.nonLinear = True

        self.Kp = np.diag([self.nlk_p1,self.nlk_p2,self.nlk_p3])
        self.Kd = np.diag([self.nlk_d1,self.nlk_d2,self.nlk_d3])

        self.Kr = np.diag([self.nlk_pphi,self.nlk_ptheta,self.nlk_ppsi])
        self.Kw = np.diag([self.nlk_dphi,self.nlk_dtheta,self.nlk_dpsi])


        self.cmd_matrix = np.array([[1,1,1,1],[0,self.arm_length,0,-self.arm_length],\
            [-self.arm_length,0 ,self.arm_length,0],[self.gamma,-self.gamma,self.gamma,-self.gamma]])

    def update(self, t, state, flat_output):
        """
        This function receives the current time, true state, and desired flat
        outputs. It returns the command inputs.

        Inputs:
            t, present time in seconds
            state, a dict describing the present state with keys
                x, position, m
                v, linear velocity, m/s
                q, quaternion [i,j,k,w]
                w, angular velocity, rad/s
            flat_output, a dict describing the present desired flat outputs with keys
                x,        position, m
                x_dot,    velocity, m/s
                x_ddot,   acceleration, m/s**2
                x_dddot,  jerk, m/s**3
                x_ddddot, snap, m/s**4
                yaw,      yaw angle, rad
                yaw_dot,  yaw rate, rad/s

        Outputs:
            control_input, a dict describing the present computed control inputs with keys
                cmd_motor_speeds, rad/sya
                cmd_thrust, N (for debugging and laboratory; not used by simulator)
                cmd_moment, N*m (for debugging; not used by simulator)
                cmd_q, quaternion [i,j,k,w] (for laboratory; not used by simulator)
        """
        cmd_motor_speeds = np.zeros((4,))
        cmd_thrust = 0
        cmd_moment = np.zeros((3,))
        cmd_q = np.zeros((4,))

        # STUDENT CODE HERE

        x_cmd = flat_output['x']
        xd_cmd = flat_output['x_dot']
        xdd_cmd = flat_output['x_dddot']


        x = state['x']
        xd = state['v']
        R = Rotation.from_quat(state['q'])
        roll,pitch,yaw = R.as_euler('xyz',degrees = False)

        R_mat = R.as_matrix()
        #self.euler_from_quaternion(state['q'][0],state['q'][1],state['q'][2],state['q'][3])
        w = state['w']
        w_des = np.array([0,0,flat_output['yaw_dot']])
        yaw_des = flat_output['yaw']

        if self.nonLinear:
            rdd_des =  xdd_cmd - self.Kd@(xd-xd_cmd)-self.Kp@(x-x_cmd)
            #print(rdd_des)            
            # print(xdd_cmd)
            Fdes = self.mass*(rdd_des + np.array([0,0,self.g]))
            cmd_thrust = np.array([R_mat[:,2].T@Fdes])
            Rdes_mat = np.zeros((3,3))
            Rdes_mat[:,2] = Fdes/np.linalg.norm(Fdes)
            
            yawvec = np.array([np.cos(yaw_des),np.sin(yaw_des),0])
            Rdes_mat[:,1] = np.cross(Rdes_mat[:,2],yawvec)/np.linalg.norm(np.cross(Rdes_mat[:,2],yawvec))
            Rdes_mat[:,0] = np.cross(Rdes_mat[:,1],Rdes_mat[:,2])
            # print(Rdes_mat)
            R_error = 0.5*(Rdes_mat.T@R_mat - R_mat.T@Rdes_mat)
            R_errorvec = np.array([-R_error[1,2],R_error[0,2],-R_error[0,1]])
            cmd_q=Rotation.from_matrix(Rdes_mat).as_quat()
            #print(R_errorvec)
            u2 = (-self.Kr@R_errorvec - self.Kw@(np.array([w[0]-w_des[0],w[1]-w_des[1],w[2]-w_des[2]]))) #+ np.cross(w,self.inertia@w)
            cmd_moment = self.inertia@u2
        else:
            r1dd_des =  xdd_cmd[0] - self.k_d1*(xd[0]-xd_cmd[0])-self.k_p1*(x[0]-x_cmd[0])
            r2dd_des = xdd_cmd[1] - self.k_d2*(xd[1]-xd_cmd[1])-self.k_p2*(x[1]-x_cmd[1])
            r3dd_des = xdd_cmd[2] - self.k_d3*(xd[2]-xd_cmd[2])-self.k_p3*(x[2]-x_cmd[2])
            cmd_thrust = np.array([self.mass*(r3dd_des + self.g)])
            pitch_des = (r1dd_des*np.cos(yaw_des) + r2dd_des*np.sin(yaw_des))/self.g
            roll_des = (r1dd_des*np.sin(yaw_des) - r2dd_des*np.cos(yaw_des))/self.g
            u2 = np.zeros((3,))
            u2[0] = -self.k_pphi*(roll-roll_des) - self.k_dphi*(w[0] - w_des[0])
            u2[1] = -self.k_ptheta*(pitch-pitch_des) - self.k_dtheta*(w[1] - w_des[1])
            u2[2] = -self.k_ppsi*(yaw-yaw_des) - self.k_dpsi*(w[2] - w_des[2])
            
            cmd_moment = self.inertia@u2
           

        cmd_control = np.concatenate((cmd_thrust,cmd_moment),axis=0)
        cmd_forces = np.linalg.inv(self.cmd_matrix)@cmd_control
        cmd_forces[cmd_forces<0] = 0
        cmd_motor_speeds = np.sqrt(cmd_forces/self.k_thrust)

        #print(cmd_control)
        #print(x_cmd)

        control_input = {'cmd_motor_speeds':cmd_motor_speeds,
                         'cmd_thrust':cmd_thrust,
                         'cmd_moment':cmd_moment,
                         'cmd_q':cmd_q}
        return control_input
