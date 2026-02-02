/*******************************************************
 * Copyright (C) 2025, IUSL
 * Author: zhenke(zhenke@mail.nwpu.edu.cn)
 *  This file is based on VINS.
*******************************************************/

#include "visualization.h"

using namespace ros;
using namespace Eigen;
ros::Publisher pub_odometry, pub_latest_odometry , pub_final_odometry , global_pub_odometry;
ros::Publisher pub_path , global_pub_path , gps_pub_path;
ros::Publisher pub_point_cloud, pub_margin_cloud;
ros::Publisher pub_key_poses;
ros::Publisher pub_camera_pose;
ros::Publisher pub_camera_pose_visual;
nav_msgs::Path path , global_path , gps_path;

ros::Publisher pub_keyframe_pose;
ros::Publisher pub_keyframe_point;
ros::Publisher pub_extrinsic;

ros::Publisher pub_image_track;

CameraPoseVisualization cameraposevisual(1, 0, 0, 1);
static double sum_of_path = 0;
static Vector3d last_path(0.0, 0.0, 0.0);

size_t pub_counter = 0;

void registerPub(ros::NodeHandle &n)
{
    pub_latest_odometry = n.advertise<nav_msgs::Odometry>("imu_propagate", 1000);
    pub_final_odometry = n.advertise<nav_msgs::Odometry>("final_odom", 1000);
    pub_path = n.advertise<nav_msgs::Path>("path", 1000);
    global_pub_path = n.advertise<nav_msgs::Path>("fixed_path", 1000);
    gps_pub_path = n.advertise<nav_msgs::Path>("gps_raw_path", 1000);
    pub_odometry = n.advertise<nav_msgs::Odometry>("odometry", 1000);
    global_pub_odometry = n.advertise<nav_msgs::Odometry>("fixed_odom", 1000);
    pub_point_cloud = n.advertise<sensor_msgs::PointCloud>("point_cloud", 1000);
    pub_margin_cloud = n.advertise<sensor_msgs::PointCloud>("margin_cloud", 1000);
    pub_key_poses = n.advertise<visualization_msgs::Marker>("key_poses", 1000);
    pub_camera_pose = n.advertise<nav_msgs::Odometry>("camera_pose", 1000);
    pub_camera_pose_visual = n.advertise<visualization_msgs::MarkerArray>("camera_pose_visual", 1000);
    pub_keyframe_pose = n.advertise<nav_msgs::Odometry>("keyframe_pose", 1000);
    pub_keyframe_point = n.advertise<sensor_msgs::PointCloud>("keyframe_point", 1000);
    pub_extrinsic = n.advertise<nav_msgs::Odometry>("extrinsic", 1000);
    pub_image_track = n.advertise<sensor_msgs::Image>("image_track", 1000);

    cameraposevisual.setScale(0.1);
    cameraposevisual.setLineWidth(0.01);
}


void pubLatestOdometry(const Eigen::Vector3d &P, const Eigen::Quaterniond &Q, const Eigen::Vector3d &V, double t)
{
    nav_msgs::Odometry odometry;
    odometry.header.stamp = ros::Time(t);
    odometry.header.frame_id = "world";
    odometry.pose.pose.position.x = P.x();
    odometry.pose.pose.position.y = P.y();
    odometry.pose.pose.position.z = P.z();
    odometry.pose.pose.orientation.x = Q.x();
    odometry.pose.pose.orientation.y = Q.y();
    odometry.pose.pose.orientation.z = Q.z();
    odometry.pose.pose.orientation.w = Q.w();
    odometry.twist.twist.linear.x = V.x();
    odometry.twist.twist.linear.y = V.y();
    odometry.twist.twist.linear.z = V.z();
    pub_latest_odometry.publish(odometry);
}

void global_pubLatestOdometry(const Eigen::Vector3d &P, const Eigen::Quaterniond &Q, const Eigen::Vector3d &V, double t)
{
    nav_msgs::Odometry odometry;
    odometry.header.stamp = ros::Time(t);
    odometry.header.frame_id = "world";
    odometry.pose.pose.position.x = P.x();
    odometry.pose.pose.position.y = P.y();
    odometry.pose.pose.position.z = P.z();
    odometry.pose.pose.orientation.x = Q.x();
    odometry.pose.pose.orientation.y = Q.y();
    odometry.pose.pose.orientation.z = Q.z();
    odometry.pose.pose.orientation.w = Q.w();
    odometry.twist.twist.linear.x = V.x();
    odometry.twist.twist.linear.y = V.y();
    odometry.twist.twist.linear.z = V.z();
    pub_final_odometry.publish(odometry);
}

void pubTrackImage(const cv::Mat &imgTrack, const double t)
{
    std_msgs::Header header;
    header.frame_id = "world";
    header.stamp = ros::Time(t);
    sensor_msgs::ImagePtr imgTrackMsg = cv_bridge::CvImage(header, "bgr8", imgTrack).toImageMsg();
    pub_image_track.publish(imgTrackMsg);
}


void printStatistics(const Estimator &estimator, double t)
{
    if (estimator.solver_flag != Estimator::SolverFlag::NON_LINEAR)
        return;
    //printf("position: %f, %f, %f\r", estimator.Ps[WINDOW_SIZE].x(), estimator.Ps[WINDOW_SIZE].y(), estimator.Ps[WINDOW_SIZE].z());
    ROS_DEBUG_STREAM("position: " << estimator.Ps[WINDOW_SIZE].transpose());
    ROS_DEBUG_STREAM("orientation: " << estimator.Vs[WINDOW_SIZE].transpose());
    if (ESTIMATE_EXTRINSIC)
    {
        cv::FileStorage fs(EX_CALIB_RESULT_PATH, cv::FileStorage::WRITE);
        for (int i = 0; i < NUM_OF_CAM; i++)
        {
            //ROS_DEBUG("calibration result for camera %d", i);
            ROS_DEBUG_STREAM("extirnsic tic: " << estimator.tic[i].transpose());
            ROS_DEBUG_STREAM("extrinsic ric: " << Utility::R2ypr(estimator.ric[i]).transpose());

            Eigen::Matrix4d eigen_T = Eigen::Matrix4d::Identity();
            eigen_T.block<3, 3>(0, 0) = estimator.ric[i];
            eigen_T.block<3, 1>(0, 3) = estimator.tic[i];
            cv::Mat cv_T;
            cv::eigen2cv(eigen_T, cv_T);
            if(i == 0)
                fs << "body_T_cam0" << cv_T ;
            else
                fs << "body_T_cam1" << cv_T ;
        }
        fs.release();
    }

    static double sum_of_time = 0;
    static int sum_of_calculation = 0;
    sum_of_time += t;
    sum_of_calculation++;
    ROS_DEBUG("vo solver costs: %f ms", t);
    ROS_DEBUG("average of time %f ms", sum_of_time / sum_of_calculation);

    sum_of_path += (estimator.Ps[WINDOW_SIZE] - last_path).norm();
    last_path = estimator.Ps[WINDOW_SIZE];
    ROS_DEBUG("sum of path %f", sum_of_path);
    // if (ESTIMATE_TD)
        // ROS_INFO("td %f", estimator.td);
}

void pubOdometry(const Estimator &estimator, const std_msgs::Header &header)
{
    if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR)
    {

        nav_msgs::Odometry odometry;
        odometry.header = header;
        odometry.header.frame_id = "world";
        odometry.child_frame_id = "world";
        Quaterniond tmp_Q;
        tmp_Q = Quaterniond(estimator.Rs[WINDOW_SIZE]); 
//以下为自己修改的
        /*
int g=0;
  tf::TransformListener listener(ros::Duration(10));//建立一个缓冲10s钟的监听器
  geometry_msgs::PoseStamped pose;
  geometry_msgs::PoseStamped pose_point;
  pose.header.frame_id = "world";//设置其处于base_laser坐标系下，坐标变换时必须有的操作
  //pose.header.stamp = header.stamp;//设置时间戳
  pose.pose.position.x = estimator.Ps[WINDOW_SIZE].x();
  pose.pose.position.y = estimator.Ps[WINDOW_SIZE].y();
  pose.pose.position.z = estimator.Ps[WINDOW_SIZE].z();
  pose.pose.orientation.w=tmp_Q.w();
  pose.pose.orientation.x=tmp_Q.x();
  pose.pose.orientation.y=tmp_Q.y();
  pose.pose.orientation.z=tmp_Q.z();
while(1)
{
  try{ 
      ros::Time now = ros::Time(0);
    listener.waitForTransform("trans", "world",
                              now, ros::Duration(1.0));
    pose.header.stamp = now;
    listener.transformPose("trans", pose, pose_point);
    g++;
  //  ROS_INFO("world: (%.2f, %.2f. %.2f %.2f %.2f %.2f %.2f) -----> trans: (%.2f, %.2f, %.2f,%.2f,%.2f,%.2f,%.2f) at time %.2f",
   //     pose.pose.position.x, pose.pose.position.y, pose.pose.position.z,pose.pose.orientation.w,pose.pose.orientation.x,pose.pose.orientation.y,pose.pose.orientation.z,
   //     pose_point.pose.position.x, pose_point.pose.position.y, pose_point.pose.position.z,pose_point.pose.orientation.w,pose_point.pose.orientation.x,pose_point.pose.orientation.y,pose_point.pose.orientation.z,pose_point.header.stamp.toSec());//将其发布出来
   if(g==1)
   break;
  }
  catch(tf::TransformException& ex){
      ROS_ERROR("Received an exception trying to transform a point from \"world\" to \"trans\": %s", ex.what());
  }
}
*/
//以上为自己修改的
        odometry.pose.pose.position.x = estimator.Ps[WINDOW_SIZE].x();
        odometry.pose.pose.position.y = estimator.Ps[WINDOW_SIZE].y();
        odometry.pose.pose.position.z = estimator.Ps[WINDOW_SIZE].z();
        odometry.pose.pose.orientation.x = tmp_Q.x();
        odometry.pose.pose.orientation.y = tmp_Q.y();
        odometry.pose.pose.orientation.z = tmp_Q.z();
        odometry.pose.pose.orientation.w = tmp_Q.w();//原本的
//以下为自己修改    
/*    
        odometry.pose.pose.position.x = pose_point.pose.position.x;
        odometry.pose.pose.position.y = pose_point.pose.position.y;
        odometry.pose.pose.position.z = pose_point.pose.position.z;
        odometry.pose.pose.orientation.x = pose_point.pose.orientation.x;
        odometry.pose.pose.orientation.y = pose_point.pose.orientation.y;
        odometry.pose.pose.orientation.z = pose_point.pose.orientation.z;
        odometry.pose.pose.orientation.w = pose_point.pose.orientation.w;
        */
 //以上为自己修改的姿态信息      
        odometry.twist.twist.linear.x = estimator.Vs[WINDOW_SIZE].x();
        odometry.twist.twist.linear.y = estimator.Vs[WINDOW_SIZE].y();
        odometry.twist.twist.linear.z = estimator.Vs[WINDOW_SIZE].z();//原来的
//以下为自己修改
        /*
        odometry.twist.twist.linear.x = estimator.Vs[WINDOW_SIZE].y();
        odometry.twist.twist.linear.y = -estimator.Vs[WINDOW_SIZE].x();
        odometry.twist.twist.linear.z = estimator.Vs[WINDOW_SIZE].z();
        */
//以上为自己修改的
        pub_odometry.publish(odometry);

        

        geometry_msgs::PoseStamped pose_stamped;
        pose_stamped.header = header;
        pose_stamped.header.frame_id = "world";      
        pose_stamped.pose = odometry.pose.pose;

        path.header = header;
        path.header.frame_id = "world";
        path.poses.push_back(pose_stamped);
        pub_path.publish(path);

        // write result to file
        /*ofstream foutC(VINS_RESULT_PATH, ios::app);
        foutC.setf(ios::fixed, ios::floatfield);
        foutC.precision(0);
        foutC << header.stamp.toSec() * 1e9 << ",";
        foutC.precision(5);
        foutC << estimator.Ps[WINDOW_SIZE].x() << ","
              << estimator.Ps[WINDOW_SIZE].y() << ","
              << estimator.Ps[WINDOW_SIZE].z() << ","
              << tmp_Q.w() << ","
              << tmp_Q.x() << ","
              << tmp_Q.y() << ","
              << tmp_Q.z() << ","
              << estimator.Vs[WINDOW_SIZE].x() << ","
              << estimator.Vs[WINDOW_SIZE].y() << ","
              << estimator.Vs[WINDOW_SIZE].z() << "," << endl;*/
	ofstream foutC(VINS_RESULT_PATH, ios::app);
        foutC.setf(ios::fixed, ios::floatfield);
        foutC.precision(0);
        foutC << header.stamp.toSec() << " ";
        foutC.precision(5);
        foutC << estimator.Ps[WINDOW_SIZE].x() << " "
              << estimator.Ps[WINDOW_SIZE].y() << " "
              << estimator.Ps[WINDOW_SIZE].z() << " "
              << tmp_Q.x() << " "
              << tmp_Q.y() << " "
              << tmp_Q.z() << " "
              << tmp_Q.w() << endl;

        foutC.close();
        Eigen::Vector3d tmp_T = estimator.Ps[WINDOW_SIZE];
        // printf("time: %f, t: %f %f %f q: %f %f %f %f \n", header.stamp.toSec(), tmp_T.x(), tmp_T.y(), tmp_T.z(),
        //                                                   tmp_Q.w(), tmp_Q.x(), tmp_Q.y(), tmp_Q.z());
    }
}

void global_pubOdometry(const Estimator &estimator, const std_msgs::Header &header)
{
    if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR)
    {

        nav_msgs::Odometry odometry_global;
        odometry_global.header = header;
        odometry_global.header.frame_id = "world";
        odometry_global.child_frame_id = "world";
        Quaterniond tmp_Q;
        tmp_Q = estimator.final_Q; 

        odometry_global.pose.pose.position.x = estimator.final_P.x();
        odometry_global.pose.pose.position.y = estimator.final_P.y();
        odometry_global.pose.pose.position.z = estimator.final_P.z();
        odometry_global.pose.pose.orientation.x = tmp_Q.x();
        odometry_global.pose.pose.orientation.y = tmp_Q.y();
        odometry_global.pose.pose.orientation.z = tmp_Q.z();
        odometry_global.pose.pose.orientation.w = tmp_Q.w();//原本的
   
        odometry_global.twist.twist.linear.x = estimator.final_V.x();
        odometry_global.twist.twist.linear.y = estimator.final_V.y();
        odometry_global.twist.twist.linear.z = estimator.final_V.z();//原来的

        // std::ofstream foutC("/home/zkz/小论文数据支持/xiyuantiyuchang/final.csv", ios::app);
        // std::ofstream foutC("/home/zkz/小论文数据支持/xiyuantiyuchang/x_final.csv", ios::app);
        // std::ofstream foutC("/home/zk/小论文数据支持/xiyuantiyuchang/ransac_final_RT0.5.csv", ios::app);
        // std::ofstream foutC("/home/zk/小论文数据支持/xiyuanlouti1/withconv_withoutransac.csv", ios::app);

        // foutC.setf(ios::fixed, ios::floatfield);
        // foutC.precision(7);
        // foutC << estimator.t_we_need << " ";
        // foutC.precision(7);
        // foutC << estimator.final_P.x() << " "
        //         << estimator.final_P.y() << " "
        //         << estimator.final_P.z() << " "
        //         << estimator.cov << " "
        //         << 0 << " "
        //         << 0 << " "
        //         << 0 << endl;
        // foutC.close();

        // std::ofstream outfile("/home/zkz/小论文数据支持/xiyuantiyuchang/x_final.txt", std::ios::app);
        // // std::ofstream outfile("/home/zkz/小论文数据支持/xiyuantiyuchang/final.txt", std::ios::app);
        // if (!outfile.is_open()) {
        //     std::cerr << "无法打开文件！" << std::endl;
        //     return; // 返回错误代码
        // }

        // outfile <<estimator.t_we_need<<" "<<estimator.final_P.x()<<" "<<estimator.final_P.y()<<" "<<estimator.final_P.z()<<" "<<estimator.predice_or_not<< std::endl;

        // outfile.close();outfile.close();

        global_pub_odometry.publish(odometry_global);

        geometry_msgs::PoseStamped pose_stamped_global;
        pose_stamped_global.header = header;
        pose_stamped_global.header.frame_id = "world";      
        pose_stamped_global.pose = odometry_global.pose.pose;

        global_path.header = header;
        global_path.header.frame_id = "world";
        global_path.poses.push_back(pose_stamped_global);
        global_pub_path.publish(global_path);

        

	// ofstream foutC(VINS_RESULT_PATH, ios::app);
    //     foutC.setf(ios::fixed, ios::floatfield);
    //     foutC.precision(0);
    //     foutC << header.stamp.toSec() << " ";
    //     foutC.precision(5);
    //     foutC << estimator.Ps[WINDOW_SIZE].x() << " "
    //           << estimator.Ps[WINDOW_SIZE].y() << " "
    //           << estimator.Ps[WINDOW_SIZE].z() << " "
    //           << tmp_Q.x() << " "
    //           << tmp_Q.y() << " "
    //           << tmp_Q.z() << " "
    //           << tmp_Q.w() << endl;

    //     foutC.close();
    //     Eigen::Vector3d tmp_T = estimator.Ps[WINDOW_SIZE];
    //     printf("time: %f, t: %f %f %f q: %f %f %f %f \n", header.stamp.toSec(), tmp_T.x(), tmp_T.y(), tmp_T.z(),
    //                                                       tmp_Q.w(), tmp_Q.x(), tmp_Q.y(), tmp_Q.z());
    }
}

void gps_pubOdometry(const Estimator &estimator, const std_msgs::Header &header)
{
    if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR)
    {

        nav_msgs::Odometry odometry_global;
        odometry_global.header = header;
        odometry_global.header.frame_id = "world";
        odometry_global.child_frame_id = "world";
        Quaterniond tmp_Q;
        tmp_Q = estimator.final_Q; 

        odometry_global.pose.pose.position.x = estimator.gp.x();
        odometry_global.pose.pose.position.y = estimator.gp.y();
        odometry_global.pose.pose.position.z = estimator.gp.z();
        // odometry_global.pose.pose.orientation.x = tmp_Q.x();
        // odometry_global.pose.pose.orientation.y = tmp_Q.y();
        // odometry_global.pose.pose.orientation.z = tmp_Q.z();
        // odometry_global.pose.pose.orientation.w = tmp_Q.w();//原本的
   
        // odometry_global.twist.twist.linear.x = estimator.final_V.x();
        // odometry_global.twist.twist.linear.y = estimator.final_V.y();
        // odometry_global.twist.twist.linear.z = estimator.final_V.z();//原来的

        // global_pub_odometry.publish(odometry_global);

        // std::ofstream outfile("/home/zkz/小论文数据支持/xiyuantiyuchang/gt.txt", std::ios::app);
        // if (!outfile.is_open()) {
        //     std::cerr << "无法打开文件！" << std::endl;
        //     return ; // 返回错误代码
        // }

        // outfile <<estimator.t_we_need<<" "<<estimator.gp.x()<<" "<<estimator.gp.y()<<" "<<estimator.gp.z()<< std::endl;

        // outfile.close();

        std::ofstream foutC("/home/zk/小论文数据支持/xiyuantiyuchang/gt_higher.csv", ios::app);
        foutC.setf(ios::fixed, ios::floatfield);
        foutC.precision(5);
        foutC << estimator.t_we_need << " ";
        foutC.precision(5);
        foutC << estimator.gp.x() << " "
                << estimator.gp.y() << " "
                << estimator.gp.z() << " "
                << 0.5 << " "
                << 0 << " "
                << 0 << " "
                << 0 << endl;
        foutC.close();


        geometry_msgs::PoseStamped pose_stamped_gps;
        pose_stamped_gps.header = header;
        pose_stamped_gps.header.frame_id = "world";      
        pose_stamped_gps.pose = odometry_global.pose.pose;

        gps_path.header = header;
        gps_path.header.frame_id = "world";
        gps_path.poses.push_back(pose_stamped_gps);
        gps_pub_path.publish(gps_path);
    }
}

void pubKeyPoses(const Estimator &estimator, const std_msgs::Header &header)
{
    if (estimator.key_poses.size() == 0)
        return;
    visualization_msgs::Marker key_poses;
    key_poses.header = header;
    key_poses.header.frame_id = "world";
    key_poses.ns = "key_poses";
    key_poses.type = visualization_msgs::Marker::SPHERE_LIST;
    key_poses.action = visualization_msgs::Marker::ADD;
    key_poses.pose.orientation.w = 1.0;
    key_poses.lifetime = ros::Duration();

    //static int key_poses_id = 0;
    key_poses.id = 0; //key_poses_id++;
    key_poses.scale.x = 0.05;
    key_poses.scale.y = 0.05;
    key_poses.scale.z = 0.05;
    key_poses.color.r = 1.0;
    key_poses.color.a = 1.0;

    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        geometry_msgs::Point pose_marker;
        Vector3d correct_pose;
        correct_pose = estimator.key_poses[i];
        pose_marker.x = correct_pose.x();
        pose_marker.y = correct_pose.y();
        pose_marker.z = correct_pose.z();
        key_poses.points.push_back(pose_marker);
    }
    pub_key_poses.publish(key_poses);
}

void pubCameraPose(const Estimator &estimator, const std_msgs::Header &header)
{
    int idx2 = WINDOW_SIZE - 1;

    if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR)
    {
        int i = idx2;
        Vector3d P = estimator.Ps[i] + estimator.Rs[i] * estimator.tic[0];
        Quaterniond R = Quaterniond(estimator.Rs[i] * estimator.ric[0]);

        nav_msgs::Odometry odometry;
        odometry.header = header;
        odometry.header.frame_id = "world";
        odometry.pose.pose.position.x = P.x();
        odometry.pose.pose.position.y = P.y();
        odometry.pose.pose.position.z = P.z();
        odometry.pose.pose.orientation.x = R.x();
        odometry.pose.pose.orientation.y = R.y();
        odometry.pose.pose.orientation.z = R.z();
        odometry.pose.pose.orientation.w = R.w();


        pub_camera_pose.publish(odometry);

        cameraposevisual.reset();
        cameraposevisual.add_pose(P, R);
        if(STEREO)
        {
            Vector3d P = estimator.Ps[i] + estimator.Rs[i] * estimator.tic[1];
            Quaterniond R = Quaterniond(estimator.Rs[i] * estimator.ric[1]);
            cameraposevisual.add_pose(P, R);
        }
        cameraposevisual.publish_by(pub_camera_pose_visual, odometry.header);
    }
}


void pubPointCloud(const Estimator &estimator, const std_msgs::Header &header)
{
    sensor_msgs::PointCloud point_cloud, loop_point_cloud;
    point_cloud.header = header;
    loop_point_cloud.header = header;


    for (auto &it_per_id : estimator.f_manager.feature)
    {
        int used_num;
        used_num = it_per_id.feature_per_frame.size();
        if (!(used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
        if (it_per_id.start_frame > WINDOW_SIZE * 3.0 / 4.0 || it_per_id.solve_flag != 1)
            continue;
        int imu_i = it_per_id.start_frame;
        Vector3d pts_i = it_per_id.feature_per_frame[0].point * it_per_id.estimated_depth;
        Vector3d w_pts_i = estimator.Rs[imu_i] * (estimator.ric[0] * pts_i + estimator.tic[0]) + estimator.Ps[imu_i];

        geometry_msgs::Point32 p;
        p.x = w_pts_i(0);
        p.y = w_pts_i(1);
        p.z = w_pts_i(2);
        point_cloud.points.push_back(p);
    }
    pub_point_cloud.publish(point_cloud);


    // pub margined potin
    sensor_msgs::PointCloud margin_cloud;
    margin_cloud.header = header;

    for (auto &it_per_id : estimator.f_manager.feature)
    { 
        int used_num;
        used_num = it_per_id.feature_per_frame.size();
        if (!(used_num >= 2 && it_per_id.start_frame < WINDOW_SIZE - 2))
            continue;
        //if (it_per_id->start_frame > WINDOW_SIZE * 3.0 / 4.0 || it_per_id->solve_flag != 1)
        //        continue;

        if (it_per_id.start_frame == 0 && it_per_id.feature_per_frame.size() <= 2 
            && it_per_id.solve_flag == 1 )
        {
            int imu_i = it_per_id.start_frame;
            Vector3d pts_i = it_per_id.feature_per_frame[0].point * it_per_id.estimated_depth;
            Vector3d w_pts_i = estimator.Rs[imu_i] * (estimator.ric[0] * pts_i + estimator.tic[0]) + estimator.Ps[imu_i];

            geometry_msgs::Point32 p;
            p.x = w_pts_i(0);
            p.y = w_pts_i(1);
            p.z = w_pts_i(2);
            margin_cloud.points.push_back(p);
        }
    }
    pub_margin_cloud.publish(margin_cloud);
}


void pubTF(const Estimator &estimator, const std_msgs::Header &header)
{
    if( estimator.solver_flag != Estimator::SolverFlag::NON_LINEAR)
        return;
    static tf::TransformBroadcaster br;
    tf::Transform transform;
    tf::Quaternion q;
    // body frame
    Vector3d correct_t;
    Quaterniond correct_q;
    correct_t = estimator.Ps[WINDOW_SIZE];
    correct_q = estimator.Rs[WINDOW_SIZE];

    transform.setOrigin(tf::Vector3(correct_t(0),
                                    correct_t(1),
                                    correct_t(2)));
    q.setW(correct_q.w());
    q.setX(correct_q.x());
    q.setY(correct_q.y());
    q.setZ(correct_q.z());
    transform.setRotation(q);
    br.sendTransform(tf::StampedTransform(transform, header.stamp, "world", "body"));

    // camera frame
    transform.setOrigin(tf::Vector3(estimator.tic[0].x(),
                                    estimator.tic[0].y(),
                                    estimator.tic[0].z()));
    q.setW(Quaterniond(estimator.ric[0]).w());
    q.setX(Quaterniond(estimator.ric[0]).x());
    q.setY(Quaterniond(estimator.ric[0]).y());
    q.setZ(Quaterniond(estimator.ric[0]).z());
    transform.setRotation(q);
    br.sendTransform(tf::StampedTransform(transform, header.stamp, "body", "camera"));
//自己加的tf变换
    transform.setOrigin(tf::Vector3(0,0,0));
    q.setRPY(0, 0,M_PI/2);//以roll,pitch,yaw三个角度给四元数赋值
    transform.setRotation(q);
    br.sendTransform(tf::StampedTransform(transform, ros::Time(0), "world", "trans"));
/************************************************************************************************/

    nav_msgs::Odometry odometry;
    odometry.header = header;
    odometry.header.frame_id = "world";
    odometry.pose.pose.position.x = estimator.tic[0].x();
    odometry.pose.pose.position.y = estimator.tic[0].y();
    odometry.pose.pose.position.z = estimator.tic[0].z();
    Quaterniond tmp_q{estimator.ric[0]};
    odometry.pose.pose.orientation.x = tmp_q.x();
    odometry.pose.pose.orientation.y = tmp_q.y();
    odometry.pose.pose.orientation.z = tmp_q.z();
    odometry.pose.pose.orientation.w = tmp_q.w();
    pub_extrinsic.publish(odometry);

}

void pubKeyframe(const Estimator &estimator)
{
    // pub camera pose, 2D-3D points of keyframe
    if (estimator.solver_flag == Estimator::SolverFlag::NON_LINEAR && estimator.marginalization_flag == 0)
    {
        int i = WINDOW_SIZE - 2;
        //Vector3d P = estimator.Ps[i] + estimator.Rs[i] * estimator.tic[0];
        Vector3d P = estimator.Ps[i];
        Quaterniond R = Quaterniond(estimator.Rs[i]);

        nav_msgs::Odometry odometry;
        odometry.header.stamp = ros::Time(estimator.Headers[WINDOW_SIZE - 2]);
        odometry.header.frame_id = "world";
        odometry.pose.pose.position.x = P.x();
        odometry.pose.pose.position.y = P.y();
        odometry.pose.pose.position.z = P.z();
        odometry.pose.pose.orientation.x = R.x();
        odometry.pose.pose.orientation.y = R.y();
        odometry.pose.pose.orientation.z = R.z();
        odometry.pose.pose.orientation.w = R.w();
        //printf("time: %f t: %f %f %f r: %f %f %f %f\n", odometry.header.stamp.toSec(), P.x(), P.y(), P.z(), R.w(), R.x(), R.y(), R.z());

        pub_keyframe_pose.publish(odometry);


        sensor_msgs::PointCloud point_cloud;
        point_cloud.header.stamp = ros::Time(estimator.Headers[WINDOW_SIZE - 2]);
        point_cloud.header.frame_id = "world";
        for (auto &it_per_id : estimator.f_manager.feature)
        {
            int frame_size = it_per_id.feature_per_frame.size();
            if(it_per_id.start_frame < WINDOW_SIZE - 2 && it_per_id.start_frame + frame_size - 1 >= WINDOW_SIZE - 2 && it_per_id.solve_flag == 1)
            {

                int imu_i = it_per_id.start_frame;
                Vector3d pts_i = it_per_id.feature_per_frame[0].point * it_per_id.estimated_depth;
                Vector3d w_pts_i = estimator.Rs[imu_i] * (estimator.ric[0] * pts_i + estimator.tic[0])
                                      + estimator.Ps[imu_i];
                geometry_msgs::Point32 p;
                p.x = w_pts_i(0);
                p.y = w_pts_i(1);
                p.z = w_pts_i(2);
                point_cloud.points.push_back(p);

                int imu_j = WINDOW_SIZE - 2 - it_per_id.start_frame;
                sensor_msgs::ChannelFloat32 p_2d;
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].point.x());
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].point.y());
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].uv.x());
                p_2d.values.push_back(it_per_id.feature_per_frame[imu_j].uv.y());
                p_2d.values.push_back(it_per_id.feature_id);
                point_cloud.channels.push_back(p_2d);
            }

        }
        pub_keyframe_point.publish(point_cloud);
    }
}
