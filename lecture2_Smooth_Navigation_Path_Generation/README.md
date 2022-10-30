sudo apt update  
sudo apt install libompl-dev  
git clone https://github.com/xinya-9527/Numerical-Optimization-in-Robotics-ShenLan.git  
cd Numerical-Optimization-in-Robotics-ShenLan/lecture2_Smooth_Navigation_Path_Generation  
catkin_make  
source devel/setup.bash  
roslaunch gcopter curve_gen.launch  

