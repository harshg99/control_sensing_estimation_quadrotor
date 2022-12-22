import numpy as np

from .graph_search import graph_search
from pdb import set_trace as T
import scipy.optimize as opt


class WorldTraj(object):
    """

    """
    def __init__(self, world, start, goal):
        """
        This is the constructor for the trajectory object. A fresh trajectory
        object will be constructed before each mission. For a world trajectory,
        the input arguments are start and end positions and a world object. You
        are free to choose the path taken in any way you like.

        You should initialize parameters and pre-compute values such as
        polynomial coefficients here.

        Parameters:
            world, World object representing the environment obstacles
            start, xyz position in meters, shape=(3,)
            goal,  xyz position in meters, shape=(3,)

        """

        # You must choose resolution and margin parameters to use for path
        # planning. In the previous project these were provided to you; now you
        # must chose them for yourself. Your may try these default values, but
        # you should experiment with them!
        self.resolution = np.array([0.1, 0.1, 0.2])
        self.margin = 0.60

        # You must store the dense path returned from your Dijkstra or AStar
        # graph search algorithm as an object member. You will need it for
        # debugging, it will be used when plotting results.
        self.path, _ = graph_search(world, self.resolution, self.margin, start, goal, astar=True)

        # You must generate a sparse set of waypoints to fly between. Your
        # original Dijkstra or AStar path probably has too many points that are
        # too close together. Store these waypoints as a class member; you will
        # need it for debugging and it will be used when plotting results.
        #T()
        #print(self.path)
        self.Dmax = 6.5
        self.points = self.simplify_path(self.path)

        #print(self.points.shape)
        # adding first points for better localisation
        #self.points = np.array([self.points.tolist()[0]]+self.points.tolist())
        self.N = self.points.shape[0]
        print(self.points)

        # Finally, you must compute a trajectory through the waypoints similar
        # to your task in the first project. One possibility is to use the
        # WaypointTraj object you already wrote in the first project. However,
        # you probably need to improve it using techniques we have learned this
        # semester.

        self.trajectory_points = self.points
        self.distances = np.zeros(self.points.shape[0]-1)
        self.polynomials = 1
        #self.max_velocity = 0.1
        # Expected Arrival times
        self.expected_arrivals = np.zeros(self.points.shape[0])



        if self.points.shape[0]!=1:
                self.unit_vectors = np.zeros((self.points.shape[0]-1,3))

        for j in range(self.points.shape[0]-1):
            self.distances[j] = self.euclidean_distance(self.trajectory_points[j],self.trajectory_points[j+1])
            if self.points.shape[0]!=1:
                self.unit_vectors[j] = (self.trajectory_points[j+1,:] - self.trajectory_points[j,:])/self.distances[j]         

        
        self.total_distance = self.distances.sum()
        self.vMax = 3.2
        self.tMax = self.total_distance/self.vMax
        self.velocity_min = self.total_distance/self.tMax
        #T()
        for j in range(self.points.shape[0]-1):
            velocity_min = self.velocity_min
            if j==0:
                self.expected_arrivals[j+1] = 0.3
            else:    
                self.expected_arrivals[j+1] = self.expected_arrivals[j] + self.distances[j]/velocity_min

        self.min_accel = 3.6
        self.max_vel = 6.0
        self.order = 5
        self.min_jerk = 5.0
        #self.spline_times = self.distances/self.velocity_min
        self.spline_times = np.sqrt(np.clip(self.distances/self.min_accel,0,None))
        #self.spline_times = np.power((self.distances/self.min_jerk),0.33)
        #self.spline_times[self.distances>2] = np.sqrt(np.clip(self.distances[self.distances>2]/self.min_accel,0,None))
        #self.spline_times[self.distances>=self.Dmax] = 1.1*self.spline_times[self.distances>=self.Dmax]
        #self.spline_times[self.distances<=1] = 0.95*self.spline_times[self.distances<=1]
        #T()

        self.T_delay = 0.25
        #T()
        #self.spline_times[0] = self.spline_times[0]*1.1
        print(self.spline_times)
        #print(self.spline_times2)
        print(self.spline_times.sum())
        #self.spline_times[1]*=1.1
        for j in range(self.unit_vectors.shape[0]-1):
            if np.dot(self.unit_vectors[j+1],self.unit_vectors[j])<=0.2:
                # if self.distances[j+1]<1.5:
                #     self.spline_times[j] = 1.1*self.spline_times[j]
                self.spline_times[j+1] = 1.1*self.spline_times[j+1]

        if self.distances[-1]<2:
            self.spline_times[-1] = self.spline_times[-1]*1.1
        else:
            self.spline_times[-1] = self.spline_times[-1]*1.1

        # parameters

        print("Finding spline params")
        self.coeffsx = self.compute_splines(self.points[0:,0],self.spline_times[0:])
        self.coeffsx = self.coeffsx.reshape((self.N - 1, -1))
        self.coeffsy= self.compute_splines(self.points[0:,1],self.spline_times[0:])
        self.coeffsy = self.coeffsy.reshape((self.N - 1, -1))
        self.coeffsz = self.compute_splines(self.points[0:,2],self.spline_times[0:])
        self.coeffsz = self.coeffsz.reshape((self.N - 1, -1))


        # self.polynomials_coefficients = np.zeros((self.trajectory_points.shape[0],3,self.polynomials+1))
        # self.polynomials_coefficients[self.polynomials_coefficients.shape[0]-1,:,0]=self.trajectory_points[-1]
        # for j in range(self.polynomials_coefficients.shape[0]-1):
        #     idx = self.polynomials_coefficients.shape[0]-2-j
        #     if self.polynomials==1:
        #         self.polynomials_coefficients[idx,:,0] = self.trajectory_points[idx]
        #         self.polynomials_coefficients[idx,:,1] = (self.polynomials_coefficients[idx+1,:,0] - self.polynomials_coefficients[idx,:,0])\
        #         /(self.expected_arrivals[j+1]-self.expected_arrivals[j])
        #     elif self.polynomials==2:
        #         self.polynomials_coefficients[idx,:,0] = self.trajectory_points[idx]
        #         self.polynomials_coefficients[idx,:,1] = 2*(self.polynomials_coefficients[idx+1,:,0] - self.polynomials_coefficients[idx,:,0])\
        #             -self.polynomials_coefficients[idx+1,:,1]
        #         self.polynomials_coefficients[idx,:,2] = (self.polynomials_coefficients[idx+1,:,1] - self.polynomials_coefficients[idx,:,1])/2
        #     elif self.polynomials==3:
        #         self.polynomials_coefficients[idx,:,0] = self.trajectory_points[idx]
        #         self.polynomials_coefficients[idx,:,2] = 3*(self.polynomials_coefficients[idx+1,:,1] - (self.polynomials_coefficients[idx+1,:,0]-self.polynomials_coefficients[idx,:,0]))\
        #          -2*self.polynomials_coefficients[idx+1,:,2]
        #         self.polynomials_coefficients[idx,:,1] = self.polynomials_coefficients[idx+1,:,1] - self.polynomials_coefficients[idx+1,:,2]\
        #         -self.polynomials_coefficients[idx,:,2]
        #         self.polynomials_coefficients[idx,:,3] = (self.polynomials_coefficients[idx+1,:,2] - self.polynomials_coefficients[idx,:,2])/3

    def simplify_path(self,path):
        #T()
        direction = path[1]-path[0]
        # waypoints = [path[0]]
        # for j in range(path.shape[0]-2):
        #     new_dir = path[j+2] - path[j+1]
        #     if np.linalg.norm(new_dir-direction)>0.01:
        #         waypoints.append(path[j+1])
        #         direction = new_dir
        # waypoints.append(path[-1])

        waypoints = np.array(self.rdp(np.array(path),0.3)) 

        simplified_waypoints = [waypoints[0].tolist()]
        for idx in range(len(waypoints)-1):
            while np.linalg.norm(simplified_waypoints[-1]-waypoints[idx+1])>=self.Dmax:
                frac = self.Dmax/np.linalg.norm(simplified_waypoints[-1]-waypoints[idx+1])
                if frac>0.75:
                    frac = 1.0
                    simplified_waypoints.append((((1-frac)*np.array(simplified_waypoints[-1])+frac*waypoints[idx+1])).tolist())
                elif frac>0.5:
                    frac = 0.5;
                    simplified_waypoints.append((((1-frac)*np.array(simplified_waypoints[-1])+frac*waypoints[idx+1])).tolist())
                else:
                    simplified_waypoints.append((((1-frac)*np.array(simplified_waypoints[-1])+frac*waypoints[idx+1])).tolist())
            if np.linalg.norm(simplified_waypoints[-1]-waypoints[idx+1])>=0.2:
                #frac=0.5
                #simplified_waypoints.append((((1-frac)*np.array(simplified_waypoints[-1])+frac*waypoints[idx+1])).tolist())
                simplified_waypoints.append(waypoints[idx+1].tolist())

        
        if path[-1].tolist() not in simplified_waypoints:
            simplified_waypoints.append(path[-1])


        return np.array(simplified_waypoints)

    def euclidean_distance(self,p1,p2):
        return np.sqrt(np.power(p1[0]-p2[0],2)+np.power(p1[1]-p2[1],2)+np.power(p1[2]-p2[2],2))   

    def rdp(self,points, epsilon):
        # get the start and end points
        #T()
        start = np.tile(np.expand_dims(points[0], axis=0), (points.shape[0], 1))
        end = np.tile(np.expand_dims(points[-1], axis=0), (points.shape[0], 1))

        # find distance from other_points to line formed by start and end
        dist_point_to_line = np.linalg.norm(np.cross(end - start, points - start, axis=-1),axis=-1) / np.linalg.norm(end - start, axis=-1)
        # get the index of the points with the largest distance
        max_idx = np.argmax(dist_point_to_line)
        max_value = dist_point_to_line[max_idx]

        result = []
        if max_value > epsilon:
            partial_results_left = self.rdp(points[:max_idx+1], epsilon)
            result += [list(i) for i in partial_results_left if list(i) not in result]
            partial_results_right = self.rdp(points[max_idx:], epsilon)
            result += [list(i) for i in partial_results_right if list(i) not in result]
        else:
            result += [points[0], points[-1]]

        return result

    def update(self, t):
        """
        Given the present time, return the desired flat output and derivatives.

        Inputs
            t, time, s
        Outputs
            flat_output, a dict describing the present desired flat outputs with keys
                x,        position, m
                x_dot,    velocity, m/s
                x_ddot,   acceleration, m/s**2
                x_dddot,  jerk, m/s**3
                x_ddddot, snap, m/s**4
                yaw,      yaw angle, rad
                yaw_dot,  yaw rate, rad/s
        """
        x        = np.zeros((3,))
        x_dot    = np.zeros((3,))
        x_ddot   = np.zeros((3,))
        x_dddot  = np.zeros((3,))
        x_ddddot = np.zeros((3,))
        yaw = 0
        yaw_dot = 0

        # STUDENT CODE HERE

        # wp_idx = self.expected_arrivals.shape[0]-1
        # for j in range(self.expected_arrivals.shape[0]-1):
        #     if self.expected_arrivals[j]<=t and self.expected_arrivals[j+1]>t:
        #         wp_idx = j


        # if wp_idx ==self.expected_arrivals.shape[0]-1:
        #     x_dot = np.zeros((3,))
        #     x = self.trajectory_points[-1]
        # else:
      
        #     if self.trajectory_points.shape[0]!=1:      
        #         x_dot_ = self.unit_vectors[wp_idx]*self.velocity_min
        #     x_dot = x_dot_*0.2
        #     x = self.trajectory_points[wp_idx] + x_dot_*(t-self.expected_arrivals[wp_idx])
        #     # if(t>self.tMax):
        #     #     x_dot = np.zeros((3,))
        #     #     x = self.trajectory_points[-1]

        #     # for j in range(self.polynomials_coefficients.shape[2]):
        #     #     if self.polynomials>=2:     
        #     #         x += self.polynomials_coefficients[wp_idx,:,j]*\
        #     #         np.power((t-self.expected_arrivals[wp_idx])/(self.expected_arrivals[wp_idx+1] - self.expected_arrivals[wp_idx]),j)
        #     #     else:
        #     #         x += self.polynomials_coefficients[wp_idx,:,j]*np.power(t-self.expected_arrivals[wp_idx],j)
        #     #     if j>=1:
        #     #         if self.polynomials>=2:
        #     #             x_dot += j*self.polynomials_coefficients[wp_idx,:,j]\
        #     #             *np.power((t-self.expected_arrivals[wp_idx])/(self.expected_arrivals[wp_idx+1] - self.expected_arrivals[wp_idx]),j-1)
        #     #         else:
        #     #             x_dot += j*self.polynomials_coefficients[wp_idx,:,j]*np.power(t-self.expected_arrivals[wp_idx],j-1)
        #     #     elif j>=2:
        #     #         if self.polynomials>=2:
        #     #             x_ddot += j*(j-1)*self.polynomials_coefficients[wp_idx,:,j]\
        #     #             *np.power((t-self.expected_arrivals[wp_idx])/(self.expected_arrivals[wp_idx+1] - self.expected_arrivals[wp_idx]),j-2)
        #     #         else:
        #     #             x_ddot += j*(j-1)*self.polynomials_coefficients[wp_idx,:,j]*np.power(t-self.expected_arrivals[wp_idx],j-2)
        #     #     elif j>=3:
        #     #         if self.polynomials>=2:
        #     #             x_dddot += j*(j-1)*(j-2)*self.polynomials_coefficients[wp_idx,:,j]\
        #     #             *np.power((t-self.expected_arrivals[wp_idx])/(self.expected_arrivals[wp_idx+1] - self.expected_arrivals[wp_idx]),j-3)
        #     #         else:
        #     #             x_dddot+=j*(j-1)*(j-2)*self.polynomials_coefficients[wp_idx,:,j]*np.power(t-self.expected_arrivals[wp_idx],j-3)
            
        # #print(t)
        # #print(x_dot.shape)
        # #print(x.shape)
        #  # STUDENT CODE HERE
        # # calculate spline segment index & handle weird cases i.e. t == np.inf
        # seg_index = np.searchsorted(np.cumsum(self.list_reach_time), t)

        t = t - self.T_delay

        if t<0:
            x = self.points[0,:]
            flat_output = {
            'x': x,
            'x_dot': x_dot,
            'x_ddot': x_ddot,
            'x_dddot': x_dddot,
            'x_ddddot': x_ddddot,
            'yaw': yaw,
            'yaw_dot': yaw_dot
            }
            flat_output = { 'x':x, 'x_dot':x_dot, 'x_ddot':x_ddot, 'x_dddot':x_dddot, 'x_ddddot':x_ddddot,
                        'yaw':yaw, 'yaw_dot':yaw_dot}
            return flat_output

        index = np.searchsorted(np.cumsum(self.spline_times), t)
        
        #for hover
        if index >= self.N - 1:
            x = self.points[-1, :]
        else:
            duration = t - np.insert(np.cumsum(self.spline_times), 0, 0)[index]
            if duration == 0:
                pos_t = np.zeros(self.order + 1)
                vel_t = np.zeros(self.order + 1)
                acc_t = np.zeros(self.order + 1)
                # the last element should be a constant, does not contain t
                #   so should never be 0
                pos_t[-1] = 1
                vel_t[-2] = 1
                acc_t[-3] = 2
            else:
                
                # if self.distances[-1]<2 and index==self.N-2:
                #     duration = duration/1.1
                pos_t = np.array(
                    [duration**(self.order - i) for i in range(self.order + 1)])
                vel_t = self.get_vel_matrix(duration,self.order) 
                
                acc_t = self.get_acc_matrix(duration,self.order)

            x[0] = self.coeffsx[index, :].T @ pos_t
            x[1] = self.coeffsy[index, :].T @ pos_t
            x[2] = self.coeffsz[index, :].T @ pos_t


            # x_dot[0] = 0.10*self.coeffsx[index, :].T @ vel_t
            # x_dot[1] = 0.10*self.coeffsy[index, :].T @ vel_t
            #x_dot[2] = 0.10*self.coeffsz[index, :].T @ vel_t
              
            # x_ddot[0] = 0.2*self.coeffsx[index, :].T @ acc_t
            # x_ddot[1] = 0.2*self.coeffsy[index, :].T @ acc_t
            #x_ddot[2] = 0.10*self.coeffsz[index, :].T @ acc_t

        flat_output = {
            'x': x,
            'x_dot': x_dot,
            'x_ddot': x_ddot,
            'x_dddot': x_dddot,
            'x_ddddot': x_ddddot,
            'yaw': yaw,
            'yaw_dot': yaw_dot
        }
        flat_output = { 'x':x, 'x_dot':x_dot, 'x_ddot':x_ddot, 'x_dddot':x_dddot, 'x_ddddot':x_ddddot,
                        'yaw':yaw, 'yaw_dot':yaw_dot}
        return flat_output

    


    # def objective(self,const, expected_arrivals,order = 5):
    #     W = order + 1
    #     obj = 0
    #     for i, Ti in enumerate(expected_arrivals):
    #         H = np.array([[100800 * Ti**7, 50400 * Ti**6, 20160 * Ti**5, 5040 * Ti**4],
    #               [50400 * Ti**6, 25920 * Ti**5, 10800 * Ti**4, 2880 * Ti**3],
    #               [20160 * Ti**5, 10800 * Ti**4, 4800 * Ti**3, 1440 * Ti**2],
    #               [5040 * Ti**4, 2880 * Ti**3, 1440 * Ti**2, 576 * Ti]])
    #         x = const[i * W:i * W + 4]
    #         obj += x.T @ H @ x
    #     return obj

    def objective(self,const, expected_arrivals,order = 5):
        W = order + 1
        obj = 0
        for i, Ti in enumerate(expected_arrivals):
            H = np.array([[720 * Ti**5, 360 * Ti**4, 120*Ti**3],
                      [360 * Ti**4, 192 * Ti**3, 72 * Ti**2],
                      [120*Ti**3, 72 * Ti**2, 36 * Ti]])
            x = const[i * W:i * W + 3]
            obj += x.T @ H @ x
        return obj


    def constraints_row(self,M,idx,num_segments):
        if not M:
            return None

        H, W = M[0].shape
        constraints = np.zeros((H, W * num_segments))
        for m, id in zip(M, idx):
            constraints[:, id * W:(id + 1) * W] = m
        return constraints

    def get_acc_matrix(self,time,order):
        arr = np.array([(order - j-1)*(order - j)*time**(order - j - 2) for j in range(order + 1)])
        #arr [-2:] = 0
        return arr

    def get_vel_matrix(self,time,order):
        arr = np.array([(order - j)*time**(order - j - 1) for j in range(order + 1)])
        #arr[-1] = 0
        return arr

    def get_jerk_matrix(self,time,order):
        arr = np.array([((order-j-1) *(order-j)*(order-j-2))*time**(order - j - 3) for j in range(order + 1)])
        #arr[-3:] = 0
        return arr

    def get_snap_matrix(self,time,order):
        arr = np.array([(order - j)*(order-j-1)*(order-j-2)*(order-j-3)*time**(order - j - 4) for j in range(order + 1)])
        #arr[-4:] = 0
        return arr

    def get_crackle_matrix(self,time,order):
        arr = np.array([(order - j)*(order-j-1)*(order-j-2)*(order-j-3)*(order-j-4)*time**(order - j - 5) for j in range(order + 1)])
        #arr[-4:] = 0
        return arr


    
    def compute_splines(self,points, times):
        
        N = len(points)
        segs  = N - 1
        order = self.order
        
        

        # define equality constraint matrix
        M_mat = []
        b_mat = []

        # define inequality constrainsts
        M_ineq_mat = []
        b_ineq_mat = []
        # initialize list of 
        constraints_mat = []

        #T()
        #coefficients
        coeffs = np.random.rand((segs) * (order + 1))
        # 1. Initial boundary condition

        loc_cons_t = np.zeros((order+1,))
        loc_cons_t[-1] = 1 
        #(order * [0] + [1])
        loc_cons_tp = np.array([times[0]**(order - i) for i in range(order + 1)])
        vel_cons_t = np.zeros((order+1,))
        vel_cons_t[-2] = 1
        acc_cons_t = np.zeros((order+1,))
        acc_cons_t[-3] = 2

        jerk_cons_t = np.zeros((order+1,))
        jerk_cons_t[-4] = 3*2
        snap_cons_t = np.zeros((order+1,))
        snap_cons_t[-5] = 4*3*2
        crack_cons_t = np.zeros((order+1,))
        crack_cons_t[-6] = 5*4*3*2
        
        M1 = np.vstack((loc_cons_t, loc_cons_tp, vel_cons_t))
        
        M_mat.append(self.constraints_row([M1], [0], segs))
        b_mat.append(np.array([points[0], points[1], 0]))

        # intermediate points
        for i in range(1, segs-1):
          # boundary conditions
            
            # location constraint
            loc_cons_tp = np.array([times[i]**(order - j) for j in range(order + 1)])
            loc_cons_mat = np.vstack((loc_cons_t, loc_cons_tp))

            M_mat.append(self.constraints_row([loc_cons_mat], [i], segs))
            b_mat.append(np.array([points[i], points[i + 1]]))

            # velocity constraitns at previous time step
            vel_cons_tp = self.get_vel_matrix(times[i-1],order)


            #acceleration constraints at previous timestep
            acc_cons_tp = self.get_acc_matrix(times[i-1],order)

            #jerk constraints at previous timestep
            jerk_cons_tp = self.get_jerk_matrix(times[i-1],order)

            #snap constraints at previous timestep
            snap_cons_tp = self.get_snap_matrix(times[i-1],order)
            
            
            vel_cons_mat = np.expand_dims(vel_cons_t,0)
            vel_final_mat = np.expand_dims(vel_cons_tp,0)

            M_mat.append(self.constraints_row([vel_cons_mat, -vel_final_mat], [i, i - 1], segs))
            b_mat.append(np.zeros(1))

            acc_cons_mat = np.expand_dims(acc_cons_t,0)
            acc_final_mat = np.expand_dims(acc_cons_tp,0)

            M_mat.append(self.constraints_row([acc_cons_mat, -acc_final_mat], [i, i - 1], segs))
            b_mat.append(np.zeros(1))



            #Higher order derivatives
            jerk_cons_mat = np.expand_dims(jerk_cons_t,0)
            jerk_final_mat = np.expand_dims(jerk_cons_tp,0)

            M_mat.append(self.constraints_row([jerk_cons_mat, -jerk_final_mat], [i, i - 1], segs))
            b_mat.append(np.zeros(1))

            snap_cons_mat = np.expand_dims(snap_cons_t,0)
            snap_final_mat = np.expand_dims(snap_cons_tp,0)

            M_mat.append(self.constraints_row([snap_cons_mat, -snap_final_mat], [i, i - 1], segs))
            b_mat.append(np.zeros(1))

            # if order==7:
            #     crack_cons_tp = self.get_crackle_matrix(times[i-1],order)
            #     crack_cons_mat = np.expand_dims(crack_cons_t,0)
            #     crack_final_mat = np.expand_dims(crack_cons_tp,0)

            #     M_mat.append(self.constraints_row([crack_cons_mat, -crack_final_mat], [i, i - 1], segs))
            #     b_mat.append(np.zeros(1))

            # M_ineq_mat.append(self.constraints_row([vel_cons_mat],[i],segs))
            # b_ineq_mat.append(np.array([self.max_vel]))
            # M_ineq_mat.append(self.constraints_row([acc_cons_mat],[i],segs))
            # b_ineq_mat.append(np.array([self.min_accel]))
            # M_ineq_mat.append(self.constraints_row([-vel_cons_mat],[i],segs))
            # b_ineq_mat.append(np.array([-self.max_vel]))
            # M_ineq_mat.append(self.constraints_row([-acc_cons_mat],[i],segs))
            # b_ineq_mat.append(np.array([-self.min_accel]))



        # 3. Final point of the entire path
        loc_cons_tp = np.array([times[-1]**(order - i) for i in range(order + 1)])

        vel_cons_tp = self.get_vel_matrix(times[-1],order)
        acc_cons_tp = self.get_vel_matrix(times[-1],order)
        jerk_cons_tp = self.get_jerk_matrix(times[-1],order)
        #snap_cons_tp = self.get_snap_matrix(times[-1],order)
        matrix = np.vstack((loc_cons_t, loc_cons_tp, vel_cons_tp))
        #matrix = np.vstack((loc_cons_t, loc_cons_tp, vel_cons_tp,jerk_cons_tp,snap_cons_tp))
        #matrix = np.vstack((loc_cons_t, loc_cons_tp, vel_cons_tp,acc_cons_tp))
        M_mat.append(self.constraints_row([matrix], [segs - 1], segs))

        b_mat.append(np.array([points[segs-1], points[segs], 0]))
        #b_mat.append(np.array([points[segs-1], points[segs], 0,0]))

        vel_cons_tp = np.expand_dims(self.get_vel_matrix(times[-2],order),0)
        acc_cons_tp = np.expand_dims(self.get_acc_matrix(times[-2],order),0)

        jerk_cons_tp = np.expand_dims(self.get_jerk_matrix(times[-2],order),0)
        snap_cons_tp = np.expand_dims(self.get_snap_matrix(times[-2],order),0)
        # Adding more self.constraints_row to make system fully determinable
        M_mat.append(self.constraints_row([np.expand_dims(vel_cons_t,0),-vel_cons_tp], [segs - 1, segs - 2], segs))

        b_mat.append(np.zeros(1))

        M_mat.append(self.constraints_row([np.expand_dims(acc_cons_t,0),-acc_cons_tp], [segs - 1, segs - 2], segs))
        
        b_mat.append(np.zeros(1))

        M_mat.append(self.constraints_row([np.expand_dims(jerk_cons_t,0),-jerk_cons_tp], [segs - 1, segs - 2], segs))
        
        b_mat.append(np.zeros(1))

        M_mat.append(self.constraints_row([np.expand_dims(snap_cons_t,0),-snap_cons_tp], [segs - 1, segs - 2], segs))
        
        b_mat.append(np.zeros(1))


        # if order==7:
        #     crack_cons_tp = self.get_crackle_matrix(times[-2],order)
        #     crack_cons_mat = np.expand_dims(crack_cons_t,0)
        #     crack_final_mat = np.expand_dims(crack_cons_tp,0)

        #     M_mat.append(self.constraints_row([crack_cons_mat, -crack_final_mat], [segs-1, segs - 2], segs))
        #     b_mat.append(np.zeros(1))

        # finalize list of self.constraints_row
        M = np.concatenate(M_mat, axis=0)
        b = np.concatenate(b_mat)
        

        constraints_mat.append({'type': 'eq', 'fun': lambda c: b - M @ c})

        # M_ineq = np.concatenate(M_ineq_mat,axis=0)
        # b_ineq = np.concatenate(b_ineq_mat)

        # constraints_mat.append({'type': 'ineq', 'fun': lambda c: b_ineq - M_ineq @ c})

        # Run optimization problem
        print("Optimising")
        #before = time.time()
        solution = opt.minimize(
            self.objective,
            coeffs,
            args=(times, order),
            constraints=constraints_mat,
            options={
                'maxiter': 10000,
                'disp': True
            })

        #after = time.time()
        #print("----- Optimization finished: {} secs.".format(after - before))
        if not solution['success']:
            print("Optimization failed: {}".format(solution['message']))
            return None
        return solution['x']