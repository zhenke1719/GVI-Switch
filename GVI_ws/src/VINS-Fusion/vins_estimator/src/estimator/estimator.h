/*******************************************************
 * Copyright (C) 2025, IUSL
 * Author: zhenke(zhenke@mail.nwpu.edu.cn)
 *  This file is based on VINS.
*******************************************************/

#pragma once
 
#include <thread>
#include <mutex>
#include <std_msgs/Header.h>
#include <std_msgs/Float32.h>
#include <ceres/ceres.h>
#include <unordered_map>
#include <queue>
#include <opencv2/core/eigen.hpp>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Geometry>
#include <unsupported/Eigen/MatrixFunctions>
#include "LocalCartesian.hpp"
#include "Geocentric.hpp"
#include "parameters.h"
#include "feature_manager.h"
#include "../utility/utility.h"
#include "../utility/tic_toc.h"
#include "../initial/solve_5pts.h"
#include "../initial/initial_sfm.h"
#include "../initial/initial_alignment.h"
#include "../initial/initial_ex_rotation.h"
#include "../factor/imu_factor.h"
#include "../factor/pose_local_parameterization.h"
#include "../factor/marginalization_factor.h"
#include "../factor/projectionTwoFrameOneCamFactor.h"
#include "../factor/projectionTwoFrameTwoCamFactor.h"
#include "../factor/projectionOneFrameTwoCamFactor.h"
#include "../featureTracker/feature_tracker.h"
#include <geometry_msgs/PoseStamped.h>
#include <ros/assert.h>
#include <random> 
#include <cmath>


class Estimator;
extern Estimator estimator;
class KalmanFilter
{
public:
KalmanFilter() = default;
GeographicLib::LocalCartesian geoConverter;

Eigen::Matrix3d error_p;
int dim_x;
int dim_z;
Eigen::MatrixXd F;
Eigen::MatrixXd Fk;
Eigen::MatrixXd Phikk_1;
Eigen::MatrixXd H;

Eigen::MatrixXd Q;
Eigen::MatrixXd P;
Eigen::MatrixXd kk;
Eigen::VectorXd x;

std::map<double, Eigen::VectorXd> history_state;
int time_num;

double update_time = -1;
double m_update_time = -1;
Eigen::Vector3d a_m, g_m;
void initState(Eigen::Vector3d yk, Eigen::Matrix3d const_est_R, double time);
void update(const Eigen::Vector3d xyz,double time,double cov,Eigen::MatrixXd covariance);
void updateESIKF(const Eigen::Vector3d yk, double time, double cov, Eigen::MatrixXd covariance);


void predict(double time);

//输出为状态一步预测结果，一步预测之后的xyz
Eigen::VectorXd getPrediction();
void keep(double time);
Eigen::Matrix3d get_hat(Vector3d v);
};

class Estimator
{
  public:
   
    Estimator();//
    ~Estimator();
    void setParameter();
    
    

    // interface
    bool extractData(double local_time, std::map<double, pair<double,Eigen::Vector3d >> &datas, std::pair<double , Eigen::Vector3d> &output, double time_threshold);//
    void initFirstPose(Eigen::Vector3d p, Eigen::Matrix3d r);
    void inputIMU(double t, const Vector3d &linearAcceleration, const Vector3d &angularVelocity);
    void GPS2XYZ(double longitude, double latitude, double altitude, Eigen::Vector3d &xyz);
    Eigen::Vector3d average(std::queue<pair<double,Eigen::Vector3d>> &gps_pre);
    double variance(std::queue<pair<double,Eigen::Vector3d>> &gps_pre);
    // double get_stimulate(std::queue<pair<double,Eigen::Vector3d>> &acc);
    double get_stimulate(std::queue<pair<double,Eigen::Vector3d>> &acc, std::queue<pair<double,Eigen::Vector3d>> &gyr);
    void inputGPS(double t, const Eigen::Vector3d &lon_lat_alt_temp, const Eigen::MatrixXd &conv);
    void inputFeature(double t, const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> &featureFrame);
    void inputImage(double t, const cv::Mat &_img, const cv::Mat &_img1 = cv::Mat());
    void processIMU(double t, double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity);
    void processImage(const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> &image, const double header);
    void processMeasurements();
    void changeSensorType(int use_imu, int use_stereo);

    // internal
    void clearState();
    bool initialStructure();
    bool visualInitialAlign();
    bool relativePose(Matrix3d &relative_R, Vector3d &relative_T, int &l);
    void slideWindow();
    void slideWindowNew();
    void slideWindowOld();
    void optimization();
    void vector2double();
    void double2vector();
    bool failureDetection();
    bool getIMUInterval(double t0, double t1, vector<pair<double, Eigen::Vector3d>> &accVector, 
                                              vector<pair<double, Eigen::Vector3d>> &gyrVector);
    void getPoseInWorldFrame(Eigen::Matrix4d &T);
    void getPoseInWorldFrame(int index, Eigen::Matrix4d &T);
    void predictPtsInNextFrame();
    void outliersRejection(set<int> &removeIndex);
    double reprojectionError(Matrix3d &Ri, Vector3d &Pi, Matrix3d &rici, Vector3d &tici,
                                     Matrix3d &Rj, Vector3d &Pj, Matrix3d &ricj, Vector3d &ticj, 
                                     double depth, Vector3d &uvi, Vector3d &uvj);
    void updateLatestStates();
    void fastPredictIMU(double t, Eigen::Vector3d linear_acceleration, Eigen::Vector3d angular_velocity);
    bool IMUAvailable(double t);
    void initFirstIMUPose(vector<pair<double, Eigen::Vector3d>> &accVector);

    void make_pair_to_get_R();
    bool optimize_R(std::vector<std::pair<Eigen::Vector3d, std::pair<double , Eigen::Vector3d>>> &point_vector,
                                    Eigen::Matrix3d &R,
                                    Eigen::Vector3d &t);
    //2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
    // Eigen::Vector3d processNewData(std::pair<double, Eigen::Vector3d>& new_data, std::vector<std::pair<double, Eigen::Vector3d>>& buffer);
    // void removeOutdatedData(double current_time, std::vector<std::pair<double, Eigen::Vector3d>>& buffer);
    // void validateModel(const std::vector<std::pair<double, Eigen::Vector3d>>& buffer);
    // void maintainBufferSize(std::vector<std::pair<double, Eigen::Vector3d>>& buffer, const std::pair<double, Eigen::Vector3d>& new_data);
    // void updateInliers(std::vector<std::pair<double, Eigen::Vector3d>>& buffer);
    // void initializeInliers(const std::vector<std::pair<double, Eigen::Vector3d>>& buffer);
    // void updateInliersWithNewData(const std::vector<std::pair<double, Eigen::Vector3d>>& buffer);
    // Eigen::Vector3d computeModel();
    // double distanceToModel(const Eigen::Vector3d& point, const Eigen::Vector3d& model);
    // Eigen::Vector3d getLatestValidPosition(const std::vector<std::pair<double, Eigen::Vector3d>>& buffer);
    // std::vector<std::pair<double, Eigen::Vector3d>> inliers_;
    //2222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
    double k_estimate, b_estimate, best_parameters_k, best_parameters_b;
    std::vector<std::pair<double, Eigen::Vector3d>> global_odom_before_buf, inliners_points;
    int all_in_ransac, inline_ransac;
    double inline_rate;
    // std::vector<std::pair<double, Eigen::Vector3d>> inliner_points;
    Eigen::Quaterniond euler2quaternion(Vector3d euler);

    enum SolverFlag
    {
        INITIAL,
        NON_LINEAR
    };

    enum MarginalizationFlag
    {
        MARGIN_OLD = 0,
        MARGIN_SECOND_NEW = 1
    };
    KalmanFilter kf;
    GeographicLib::LocalCartesian geoConverter;
    std::mutex mProcess;
    std::mutex mBuf;
    std::mutex mPropagate;
    queue<pair<double, Eigen::Vector3d>> accBuf;
    queue<pair<double, Eigen::Vector3d>> gyrBuf;
//------------------------------------------------------------------------
    std::map<double, Eigen::Vector3d> local_odom_buf;
    // std::map<double, Eigen::Vector3d> global_odom_buf;
    std::map<double, pair<double,Eigen::Vector3d >> global_odom_buf;
//------------------------------------------------------------------------
    queue<pair<double, map<int, vector<pair<int, Eigen::Matrix<double, 7, 1> > > > > > featureBuf;
    double prevTime, curTime;
    bool openExEstimation;

    std::thread trackThread;
    std::thread processThread;

    FeatureTracker featureTracker;

    SolverFlag solver_flag;
    MarginalizationFlag  marginalization_flag;
    Vector3d g;

    Matrix3d ric[2];
    Vector3d tic[2];

    Vector3d        Ps[(WINDOW_SIZE + 1)];
    Vector3d        Vs[(WINDOW_SIZE + 1)];
    Matrix3d        Rs[(WINDOW_SIZE + 1)];
    Vector3d        Bas[(WINDOW_SIZE + 1)];
    Vector3d        Bgs[(WINDOW_SIZE + 1)];
    double td;

    Matrix3d back_R0, last_R, last_R0;
    Vector3d back_P0, last_P, last_P0;
    double Headers[(WINDOW_SIZE + 1)];

    IntegrationBase *pre_integrations[(WINDOW_SIZE + 1)];
    Vector3d acc_0, gyr_0;

    vector<double> dt_buf[(WINDOW_SIZE + 1)];
    vector<Vector3d> linear_acceleration_buf[(WINDOW_SIZE + 1)];
    vector<Vector3d> angular_velocity_buf[(WINDOW_SIZE + 1)];

    int frame_count;
    int sum_of_outlier, sum_of_back, sum_of_front, sum_of_invalid;
    int inputImageCnt;

    FeatureManager f_manager;
    MotionEstimator m_estimator;
    InitialEXRotation initial_ex_rotation;

    bool first_imu;
    bool is_valid, is_key;
    bool failure_occur;
    bool we_loss_it;

    vector<Vector3d> point_cloud;
    vector<Vector3d> margin_cloud;
    vector<Vector3d> key_poses;
    vector<Vector3d> gp_pass;
    double initial_timestamp;
    double norm;
    bool get_ransac;


    double para_Pose[WINDOW_SIZE + 1][SIZE_POSE];
    double para_SpeedBias[WINDOW_SIZE + 1][SIZE_SPEEDBIAS];
    double para_Feature[NUM_OF_F][SIZE_FEATURE];
    double para_Ex_Pose[2][SIZE_POSE];
    double para_Retrive_Pose[SIZE_POSE];
    double para_Td[1][1];
    double para_Tr[1][1];
    double latest_time , stimulate ,predice_or_not;

    int loop_window_index;

    MarginalizationInfo *last_marginalization_info;
    vector<double *> last_marginalization_parameter_blocks;

    map<double, ImageFrame> all_image_frame;
    IntegrationBase *tmp_pre_integration;

    Eigen::Vector3d initP , V_n;
    Eigen::Matrix3d initR;
    Eigen::Vector3d est_t , const_est_t;
    Eigen::Matrix3d est_R , const_est_R;

    Eigen::VectorXd state;
    Eigen::Vector3d delta_P , delta_P_divice;
    
    Eigen::Vector3d gp;
    Eigen::Vector3d latest_P, latest_V, latest_Ba, latest_Bg, latest_acc_0, latest_gyr_0;
    Eigen::Quaterniond latest_Q , Q_n_kf , delta_q_to_fix , final_Q;
    Eigen::Vector3d global_odom_init , final_P , final_V;
    Eigen::Matrix3d R_s_c_enu , latest_R;
    Eigen::Vector3d l_to_g_p;
    Eigen::Vector3d yk_still;
    std::map<double,Eigen::Vector3d> pubodom_raw;
    std::pair<double,Eigen::Vector3d> temp_t_l_to_g_p;
    bool initFirstPoseFlag;
    bool initThreadFlag;

    double t_we_need , t_init_tag;
    double cov;
    // bool initGPS;
    private:
    struct PlaneModel {
        double a, b, c, d;          // 
        double inlier_thresh;       // 
        int min_inliers;            // 
        std::mt19937 rng;           // 
        PlaneModel(double thresh = 1.0, int min_in = 5) 
            : inlier_thresh(thresh), min_inliers(min_in) {
            std::random_device rd;
            rng.seed(rd());
        }
    };

    PlaneModel plane_model_;                // 
    std::vector<std::pair<double, Eigen::Vector3d>> pos_buffer_;  // 
    std::vector<std::pair<double, Eigen::Vector3d>> inliers_;       // 
    double max_time_diff_;                // 
    size_t buf_capacity_; 

    void initializePlaneModel(const std::vector<std::pair<double, Eigen::Vector3d>>& buffer);

    void updateInliersWithCurrentModel(const std::vector<std::pair<double, Eigen::Vector3d>>& buffer);

    void reestimatePlaneModelIfNeeded(const std::vector<std::pair<double, Eigen::Vector3d>>& buffer);

    double distanceToPlane(const Eigen::Vector3d& point);

    Eigen::Vector3d computeCentroid(const std::vector<std::pair<double, Eigen::Vector3d>>& points);

    std::vector<std::pair<double, Eigen::Vector3d>> randomSample(
        const std::vector<std::pair<double, Eigen::Vector3d>>& buffer, int n);

    void removeOutdatedData(double current_time);

    void maintainBufferSize(std::vector<std::pair<double, Eigen::Vector3d>>& buffer, 
                           const std::pair<double, Eigen::Vector3d>& new_data);

    Eigen::Vector3d processGPS(const double& timestamp, const Eigen::Vector3d& gps_pos);
    void updatePlaneModel(const std::vector<std::pair<double, Eigen::Vector3d>>& buffer);
    Eigen::Vector3d getLatestInlierPosition();

    
};


