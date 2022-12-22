# control_sensing_estimation_quadrotor

## Introduction 

The goal of the project was to develop a complete
planning and control system with state estimation using
visual inertial odometry (VIO) for a quadrotor, enabling
it to fly from an initial point to a goal point efficiently
and without collision with obstacles or crashing. 
With perfect information of the state of
the system aggressive trajectories and controllers were
feasible.
With state estimation, however, aggressive trajecto-
ries designed in prior projects are not feasible due to
tracking errors in position and attitude of the quadrotor
introduced by the Visual Inertial Odometry(VIO) that is
sensitive to themotion of the quad-rotor. Thus, this project aims to
outline the changes made in the trajectory generation
and controls of the system, in order for the quadrotor to
perform well under the constraints introduced by state
estimation.

## Controller performance

![image](https://user-images.githubusercontent.com/28558013/209144857-b078cebc-545c-4e33-a51e-07f2fab31827.png)

## Path Planning

![image](https://user-images.githubusercontent.com/28558013/209145117-ddffb1d5-f71f-47d8-8d40-033a0aa61a0d.png)

## Visual Inertial Odometry Estimation Performance

![image](https://user-images.githubusercontent.com/28558013/209145499-18e42b12-4ee2-4789-a896-2856b6b3199d.png)

## Trajectory optimisation vs Path taken by quadrotor

![image](https://user-images.githubusercontent.com/28558013/209145635-a0969293-786a-43e2-883d-bc018fd4b858.png)


