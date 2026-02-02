/*******************************************************
 * Copyright (C) 2025, IUSL
 * Author: zhenke(zhenke@mail.nwpu.edu.cn)
 *  This file is based on VINS.
*******************************************************/

#include <stdio.h>
#include <queue>
#include <map>
#include <thread>
#include <mutex>
#include <ros/ros.h>
#include <cv_bridge/cv_bridge.h>
#include <opencv2/opencv.hpp>
#include "estimator/estimator.h"
#include "estimator/parameters.h"
#include "utility/visualization.h"
#include <sensor_msgs/NavSatFix.h>
#include <sensor_msgs/NavSatStatus.h>
#include <gnss_comm/gnss_ros.hpp>
#include <gnss_comm/gnss_utility.hpp>

using namespace gnss_comm;

Estimator estimator;

queue<sensor_msgs::ImuConstPtr> imu_buf;
queue<sensor_msgs::PointCloudConstPtr> feature_buf;
queue<sensor_msgs::ImageConstPtr> img0_buf;
queue<sensor_msgs::ImageConstPtr> img1_buf;
std::mutex m_buf;
bool init_imu = false;

void img0_callback(const sensor_msgs::ImageConstPtr &img_msg)
{
    m_buf.lock();
    img0_buf.push(img_msg);
    m_buf.unlock();
}

void img1_callback(const sensor_msgs::ImageConstPtr &img_msg)
{
    m_buf.lock();
    img1_buf.push(img_msg);
    m_buf.unlock();
}

// void img0_callback(const sensor_msgs::CompressedImageConstPtr &compressed_img_msg)
// {
//     try {
//         // 将压缩图像转换为cv::Mat格式
//         cv_bridge::CvImagePtr cv_ptr = cv_bridge::toCvCopy(compressed_img_msg, sensor_msgs::image_encodings::BGR8);
        
//         // 将cv::Mat转换回sensor_msgs::Image
//         sensor_msgs::ImagePtr img_msg = cv_ptr->toImageMsg();
//         img_msg->header = compressed_img_msg->header; // 保持原有的header信息
        
//         m_buf.lock();
//         img0_buf.push(img_msg);
//         m_buf.unlock();
//     }
//     catch (cv_bridge::Exception& e) {
//         ROS_ERROR("cv_bridge exception in img0_callback: %s", e.what());
//         return;
//     }
// }

// void img1_callback(const sensor_msgs::CompressedImageConstPtr &compressed_img_msg)
// {
//     try {
//         // 将压缩图像转换为cv::Mat格式
//         cv_bridge::CvImagePtr cv_ptr = cv_bridge::toCvCopy(compressed_img_msg, sensor_msgs::image_encodings::BGR8);
        
//         // 将cv::Mat转换回sensor_msgs::Image
//         sensor_msgs::ImagePtr img_msg = cv_ptr->toImageMsg();
//         img_msg->header = compressed_img_msg->header; // 保持原有的header信息
        
//         m_buf.lock();
//         img1_buf.push(img_msg);
//         m_buf.unlock();
//     }
//     catch (cv_bridge::Exception& e) {
//         ROS_ERROR("cv_bridge exception in img1_callback: %s", e.what());
//         return;
//     }
// }

void gpsCallback(const GnssPVTSolnMsgConstPtr &GPS_msg)
{
    // double t = GPS_msg->header.stamp.toSec();
    double t = ros::Time::now().toSec();
    // estimator.t_we_need = GPS_msg->header.stamp.toSec();
    double lon = GPS_msg->longitude;//jing
    double lat = GPS_msg->latitude;//wei
    double alt = GPS_msg->altitude;
    Vector3d lon_lat_alt(lon, lat, alt);
    double temp = double(GPS_msg->h_acc);
    Eigen::MatrixXd covariance = temp*Eigen::MatrixXd::Identity(3, 3);
    // covariance(0,0) = temp;
    // cout<<"11111111111111111111111 = "<<covariance<<endl;
    estimator.inputGPS(t,lon_lat_alt,covariance);
}

// void gpsCallback(const sensor_msgs::NavSatFixConstPtr &GPS_msg)
// {
//     // double t = GPS_msg->header.stamp.toSec();
//     double t = ros::Time::now().toSec();
//     // estimator.t_we_need = GPS_msg->header.stamp.toSec();
//     double lon = GPS_msg->longitude;//jing
//     double lat = GPS_msg->latitude;//wei
//     double alt = GPS_msg->altitude;
//     Vector3d lon_lat_alt(lon, lat, alt);
//     Eigen::MatrixXd covariance = Eigen::MatrixXd::Identity(3, 3);
//     covariance << GPS_msg->position_covariance[0],GPS_msg->position_covariance[1],GPS_msg->position_covariance[2],
//         GPS_msg->position_covariance[3],GPS_msg->position_covariance[4],GPS_msg->position_covariance[5],
//         GPS_msg->position_covariance[6],GPS_msg->position_covariance[7],GPS_msg->position_covariance[8];
//     estimator.inputGPS(t,lon_lat_alt,covariance);
// }

void optCallback(const ros::TimerEvent &e)
{
    // cout<<"gps_"<<estimator.global_odom_buf.size()<<"vio_"<<estimator.local_odom_buf.size()<<endl;
    if(estimator.local_odom_buf.size()>0&&estimator.global_odom_buf.size()>0)
    {
        estimator.make_pair_to_get_R();
    }
    else
    {
        cout<<"not_ready"<<estimator.local_odom_buf.size()<<"   "<<estimator.global_odom_buf.size()<<endl;
    }
}

cv::Mat getImageFromMsg(const sensor_msgs::ImageConstPtr &img_msg)
{
    cv_bridge::CvImageConstPtr ptr;
    if (img_msg->encoding == "8UC1")
    {
        sensor_msgs::Image img;
        img.header = img_msg->header;
        img.height = img_msg->height;
        img.width = img_msg->width;
        img.is_bigendian = img_msg->is_bigendian;
        img.step = img_msg->step;
        img.data = img_msg->data;
        img.encoding = "mono8";
        ptr = cv_bridge::toCvCopy(img, sensor_msgs::image_encodings::MONO8);
    }
    else
        ptr = cv_bridge::toCvCopy(img_msg, sensor_msgs::image_encodings::MONO8);

    cv::Mat img = ptr->image.clone();
    return img;
}

// extract images with same timestamp from two topics
void sync_process()
{
    while(1)
    {
        if(STEREO)//STEREO==1
        {
            cv::Mat image0, image1;
            std_msgs::Header header;
            double time = 0;
            m_buf.lock();
            if (!img0_buf.empty() && !img1_buf.empty())//
            {
                double time0 = img0_buf.front()->header.stamp.toSec();
                double time1 = img1_buf.front()->header.stamp.toSec();
                // 0.003s sync tolerance
                if(time0 < time1 - 0.007)
                {
                    img0_buf.pop();
                    printf("throw img0\n");
                }
                else if(time0 > time1 + 0.007)
                {
                    img1_buf.pop();
                    printf("throw img1\n");
                }
                else//左右两个相机得到图像的时间小于0.003s才对其进行处理
                {
                    time = img0_buf.front()->header.stamp.toSec();
                    header = img0_buf.front()->header;
                    image0 = getImageFromMsg(img0_buf.front());
                    img0_buf.pop();
                    image1 = getImageFromMsg(img1_buf.front());
                    img1_buf.pop();
                    //printf("find img0 and img1\n");
                }
            }
            m_buf.unlock();
            if(!image0.empty())
                estimator.inputImage(time, image0, image1);//
        }//
        else
        {
            cv::Mat image;
            std_msgs::Header header;
            double time = 0;
            m_buf.lock();
            if(!img0_buf.empty())
            {
                time = img0_buf.front()->header.stamp.toSec();
                header = img0_buf.front()->header;
                image = getImageFromMsg(img0_buf.front());
                img0_buf.pop();
            }
            m_buf.unlock();
            if(!image.empty())
                estimator.inputImage(time, image);
        }

        std::chrono::milliseconds dura(2);
        std::this_thread::sleep_for(dura);//
    }
}


void imu_callback(const sensor_msgs::ImuConstPtr &imu_msg)
{
    double t = imu_msg->header.stamp.toSec();
    if(!init_imu)
    {
        estimator.t_init_tag = imu_msg->header.stamp.toSec();
        init_imu = true;
    }
    else
    {
        estimator.t_we_need = imu_msg->header.stamp.toSec() - estimator.t_init_tag;
        // cout<<estimator.t_we_need<<endl;
    }


    double dx = imu_msg->linear_acceleration.x;
    double dy = imu_msg->linear_acceleration.y;
    double dz = imu_msg->linear_acceleration.z;
    double rx = imu_msg->angular_velocity.x;
    double ry = imu_msg->angular_velocity.y;
    double rz = imu_msg->angular_velocity.z;
    Vector3d acc(dx, dy, dz);//
    Vector3d gyr(rx, ry, rz);//
    estimator.inputIMU(t, acc, gyr);
    return;
}


void feature_callback(const sensor_msgs::PointCloudConstPtr &feature_msg)
{
    map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> featureFrame;
    for (unsigned int i = 0; i < feature_msg->points.size(); i++)
    {
        int feature_id = feature_msg->channels[0].values[i];
        int camera_id = feature_msg->channels[1].values[i];
        double x = feature_msg->points[i].x;
        double y = feature_msg->points[i].y;
        double z = feature_msg->points[i].z;
        double p_u = feature_msg->channels[2].values[i];
        double p_v = feature_msg->channels[3].values[i];
        double velocity_x = feature_msg->channels[4].values[i];
        double velocity_y = feature_msg->channels[5].values[i];
        if(feature_msg->channels.size() > 5)
        {
            double gx = feature_msg->channels[6].values[i];
            double gy = feature_msg->channels[7].values[i];
            double gz = feature_msg->channels[8].values[i];
            pts_gt[feature_id] = Eigen::Vector3d(gx, gy, gz);
            //printf("receive pts gt %d %f %f %f\n", feature_id, gx, gy, gz);
        }
        ROS_ASSERT(z == 1);//
        Eigen::Matrix<double, 7, 1> xyz_uv_velocity;
        xyz_uv_velocity << x, y, z, p_u, p_v, velocity_x, velocity_y;
        featureFrame[feature_id].emplace_back(camera_id,  xyz_uv_velocity);//
    }
    double t = feature_msg->header.stamp.toSec();
    estimator.inputFeature(t, featureFrame);//
    return;
}

void restart_callback(const std_msgs::BoolConstPtr &restart_msg)
{
    if (restart_msg->data == true)
    {
        ROS_WARN("restart the estimator!");
        estimator.clearState();
        estimator.setParameter();
    }
    return;
}

void imu_switch_callback(const std_msgs::BoolConstPtr &switch_msg)
{
    if (switch_msg->data == true)
    {
        //ROS_WARN("use IMU!");
        estimator.changeSensorType(1, STEREO);
    }
    else
    {
        //ROS_WARN("disable IMU!");
        estimator.changeSensorType(0, STEREO);
    }
    return;
}

void cam_switch_callback(const std_msgs::BoolConstPtr &switch_msg)
{
    if (switch_msg->data == true)
    {
        //ROS_WARN("use stereo!");
        estimator.changeSensorType(USE_IMU, 1);
    }
    else
    {
        //ROS_WARN("use mono camera (left)!");
        estimator.changeSensorType(USE_IMU, 0);
    }
    return;
}

int main(int argc, char **argv)
{
    ros::init(argc, argv, "vins_estimator");
    ros::NodeHandle n("~");
    ros::console::set_logger_level(ROSCONSOLE_DEFAULT_NAME, ros::console::levels::Info);//

    if(argc != 2)
    {
        printf("please intput: rosrun vins vins_node [config file] \n"
               "for example: rosrun vins vins_node "
               "~/catkin_ws/src/VINS-Fusion/config/euroc/euroc_stereo_imu_config.yaml \n");
        return 1;
    }

    string config_file = argv[1];
    printf("config_file: %s\n", argv[1]);

    readParameters(config_file);//
    estimator.setParameter();//

#ifdef EIGEN_DONT_PARALLELIZE
    ROS_DEBUG("EIGEN_DONT_PARALLELIZE");
#endif

    ROS_WARN("waiting for image and imu...");

    registerPub(n);

    ros::Subscriber sub_imu = n.subscribe(IMU_TOPIC, 2000, imu_callback, ros::TransportHints().tcpNoDelay());//
    ros::Subscriber sub_feature = n.subscribe("/feature_tracker/feature", 2000, feature_callback);//
    ros::Subscriber sub_img0 = n.subscribe(IMAGE0_TOPIC, 100, img0_callback);//
    ros::Subscriber sub_img1 = n.subscribe(IMAGE1_TOPIC, 100, img1_callback);
    // ros::Subscriber sub_gnss = n.subscribe("/GPS", 20,gpsCallback);
    // ros::Subscriber sub_gnss = n.subscribe("/mavros/global_position/raw/fix", 20,gpsCallback);
    // ros::Subscriber sub_gnss = n.subscribe("/x36d/gnss", 20,gpsCallback);
    ros::Subscriber sub_gnss = n.subscribe("/ublox_driver/receiver_pvt", 20,gpsCallback);
    ros::Subscriber sub_restart = n.subscribe("/vins_restart", 100, restart_callback);//
    ros::Subscriber sub_imu_switch = n.subscribe("/vins_imu_switch", 100, imu_switch_callback);//
    ros::Subscriber sub_cam_switch = n.subscribe("/vins_cam_switch", 100, cam_switch_callback);//
    ros::Timer opt_timer = n.createTimer(ros::Duration(0.03), &optCallback);

    
    std::thread sync_thread{sync_process};//对图片信息进行处理
    ros::spin();
    return 0;
}
