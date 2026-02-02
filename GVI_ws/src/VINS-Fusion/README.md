# GVI-switch
## A Robust Multi-Sensor Loosely-Coupled Pose Estimation Based on Error-State Kalman Filter​

## 1. Prerequisites
The compiler should at least support c++14 standards.
### 1.1 **Ubuntu** and **ROS**
GVI-switch is developed under [ROS-Noetic](http://wiki.ros.org/noetic/Installation/Ubuntu) with its default OpenCV4 library. However, GVI-switch should be working on ROS-Melodic with OpenCV3.

### 1.2. **Ceres Solver**
Follow [Ceres Installation](http://ceres-solver.org/installation.html).

### 1.3. **Eigen 3.3.7**
Eigen is a fantastic matrix computation library. GVI-switch is developed under [Eigen3.3.7](https://eigen.tuxfamily.org/index.php?title=Main_Page). Other Eigen 3 versions should be OK for GVI-switch.


## 2. Build GVI-switch
Clone the repository and catkin_make:
```
    cd ~/GVI_ws/src
    git clone https://github.com/zhenke1719/GVI-Switch.git
    cd ../
    catkin_make
    source ~/GVI_ws/devel/setup.bash
```
(if you fail in this step, try to find another computer with clean system or reinstall Ubuntu and ROS)

## 3. Run with your devices Stereo cameras + IMU + GNSS or Run bags
cd ~/GVI_ws
source ~/GVI_ws/devel/setup.bash
roslaunch vins vins_rviz.launch

open a new terminal
source ~/GVI_ws/devel/setup.bash
rosrun vins vins_node /workspace/GVI_ws/src/VINS-Fusion/config/realsense_d435i/realsense_stereo_imu_configqian.yaml

## 4. License
This software is open-sourced under the MIT License. If you use our code or datasets, please cite our paper.


