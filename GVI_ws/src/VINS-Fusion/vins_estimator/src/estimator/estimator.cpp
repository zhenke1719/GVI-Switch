/*******************************************************
 * Copyright (C) 2025, IUSL
 * Author: zhenke(zhenke@mail.nwpu.edu.cn)
 *  This file is based on VINS.
*******************************************************/

#include "estimator.h"
#include "../utility/visualization.h"
#include <cmath>
bool initGPS = false;
bool optimized = false;

bool get_Q = false;
bool get_P = false;
bool has_been_indoors = false;
bool has_been_in = false;
bool optimize_twice = false;
int pat = 1 , times_in_threshold = 0;
double threshold = 0.12;

bool stimulating = false;

static double last_valid_gnss_time = -1.0;
// double  time_of_vio_used, time_of_gps_used;
// Estimator estimator1;


bool comparePairs(const std::pair<Eigen::Vector3d, std::pair<double, Eigen::Vector3d>>& a,
                  const std::pair<Eigen::Vector3d, std::pair<double, Eigen::Vector3d>>& b) 
{
    return a.second.first < b.second.first; // 
}

void sortByDouble(std::vector<std::pair<Eigen::Vector3d, std::pair<double, Eigen::Vector3d>>>& point_vector) 
{
    std::sort(point_vector.begin(), point_vector.end(), comparePairs);
}

void KalmanFilter::initState(Eigen::Vector3d yk, Eigen::Matrix3d const_est_R, double time)
{
    // cout<<"start_to_init111111111111111111111111111111111111111111111111111111111111111111"<<endl;
    double dt = 0.03;
    x = Eigen::VectorXd(15);
    x << yk, Eigen::VectorXd::Zero(12,1);
    F = Eigen::MatrixXd::Zero(15, 15);
    // Eigen::Matrix3d a_m_hat , g_m_hat;
    // a_m_hat << 0,-a_m(2),a_m(1).
    //            a_m(2),0,-a_m(0),
    //            -a_m(1),a_m(0),0;
    // g_m_hat << 0,-g_m(2),g_m(1).
    //            g_m(2),0,-g_m(0),
    //            -g_m(1),g_m(0),0;
    F << Eigen::MatrixXd::Zero(3,3), Eigen::MatrixXd::Identity(3,3), Eigen::MatrixXd::Zero(3,9),
        Eigen::MatrixXd::Zero(3,6), -const_est_R*get_hat(a_m), -const_est_R, Eigen::MatrixXd::Zero(3,3),
        Eigen::MatrixXd::Zero(3,6), -get_hat(g_m), Eigen::MatrixXd::Zero(3,3), -Eigen::MatrixXd::Identity(3,3),
        Eigen::MatrixXd::Zero(6,15);

    Fk = F*dt;
    Phikk_1 = Eigen::MatrixXd::Identity(15, 15) + Fk;
    Phikk_1.block(6,6,3,3) = -get_hat(g_m)*dt;
    H = Eigen::MatrixXd::Zero(3, 15);
    H << 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
        
    error_p = 0*Eigen::Matrix3d::Identity();
    Q = 1e-2*Eigen::MatrixXd::Identity(15, 15);
    P = Eigen::MatrixXd::Identity(15, 15);
    // kk.setZero();
    update_time = time;  //
    // history_state[time] = x;  //

    // time_num = history_state.size();
    // cout<<"init_successfully22222222222222222222222222222222222222222222222222222222222222222"<<endl;
}

void KalmanFilter::update(const Eigen::Vector3d xyz,double time,double cov,Eigen::MatrixXd covariance)
{
    // cout<<"start_to_update5555555555555555555555555555555555555555555555555555555555555555"<<endl;
    // TicToc timer2;
    // timer2.tic();
    //-----------------------------------------------------------------------------------------------
    // double dt = time - update_time;
    // cout<<"delta_t ======update======= "<<dt<<endl;
    // F.block(3,6,3,3) = -estimator1.const_est_R*get_hat(a_m);
    // F.block(6,6,3,3) = -get_hat(g_m);
    // Fk = F*dt;
    // Phikk_1 = Eigen::MatrixXd::Identity(15, 15) + Fk;
    // Phikk_1.block(6,6,3,3) = -get_hat(g_m)*dt;
    //------------------------------------------------------------------------------------------------
    // x = Phikk_1 * x;    //状态一步预测    X_(k/k-1) = F_(k/k-1)*X_(k-1)   第一个方程
    // P = F * P * F.transpose() + Q;  //状态一步预测均方误差   P方差阵递推   P_(k/k-1) = F_(k/k-1) * P_(k-1) * F_(k/k-1)^(T) + Q_(k-1)   第二个方程
    Eigen::Vector3d y = xyz - H * x;   // y = Z_k - H_k * X_(k/k-1)
    // Eigen::MatrixXd S = H * P * H.transpose() + cov*Eigen::MatrixXd::Identity(3, 3);
    // Eigen::MatrixXd S = H * P * H.transpose() + cov*covariance;
    Eigen::MatrixXd S = H * P * H.transpose() + cov*covariance + error_p;
    // cout<<"start_to_update6666666666666666666666666666666666666666666666666666666666666666"<<endl;

       //S = H_k * P_(k/k-1) * H_k^(T) + R_k
    // cout<<"kk kk00000000000000000000000000000000000000000000000000000= "<<kk.block(3,0,3,3)<<endl; 
    Eigen::MatrixXd K = P * H.transpose() * S.inverse();  //滤波增益   K增益矩阵   K_k = P_(k-1) * H_k^(T) * S^(-1)   第三个方程
    kk =K;
    // cout<<"增增又益益k = "<<kk<<endl; 
    // cout<<"V_x bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb= "<<x(3)<<endl; 
    // cout<<"kk kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk= "<<kk.block(3,0,3,3)<<endl; 
    x = x + K * y;  //状态估计   X_K = X_(k/k-1) + K_k * y     第四个方程

    // cout<<"zzz_y 000000000000000000000000000000000000000000000000000000= "<<y<<endl; 
    // cout<<"zzz_x aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa= "<<(K * y)(3)<<endl; //这个值4000+不对
    // cout<<"V_x ccccccccccccccccccccccccccccccccccccccccccccccccccccccc= "<<x(3)<<endl; 
    P = P - K * H * P;  //状态估计均方误差  P_k = P_(k/k-1) - K_k * H_k * P_(k/k-1)    第五个方程
    // cout<<"估计误差p"<<H*P*H.transpose()<<endl;
    // cout<<"测量误差R"<<covariance<<endl;
    update_time = time; //把当前time更新为update_time
    // sleep(100);
    // double temp_time = timer2.toc();
    // std::ofstream foutC("/home/zk/小论文数据支持/cheku/all_upd_tlist.csv", ios::app);
    // foutC.setf(ios::fixed, ios::floatfield);
    // foutC.precision(6);
    // foutC << temp_time << endl;
    // foutC.close();
    // history_state[time] = x;  //把当前xyz更新为history_state

    // time_num = history_state.size();  //time_num更新次数
    // cout<<"update_successfully66666666666666666666666666666666666666666666666666666666666666666666"<<endl;
}

void KalmanFilter::updateESIKF(const Eigen::Vector3d yk, double time, double cov, Eigen::MatrixXd covariance)
{
    // ESIKF 迭代更新
    // yk: 观测残差 (VIO位置 - GPS位置)
    // x: 误差状态向量 (Error State)
    
    int max_iters = 5; // 最大迭代次数
    Eigen::VectorXd x_prior = x; // 保存先验状态 x_{k|k-1}
    Eigen::VectorXd x_curr = x;  // 当前迭代的估计状态
    Eigen::MatrixXd K_curr;      // 当前卡尔曼增益
    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(15, 15);

    for(int i = 0; i < max_iters; i++)
    {
        // 1. 计算当前估计下的残差
        // 观测方程: y = H * x + n
        // residual = measurement - H * x_curr
        Eigen::Vector3d residual = yk - H * x_curr;

        // 2. 计算 S 矩阵 (H 在此处是常数矩阵，若为非线性需在此处重新线性化)
        // R = cov * covariance + error_p
        Eigen::MatrixXd R = cov * covariance;
        // Eigen::MatrixXd R = cov * covariance + error_p;
        Eigen::MatrixXd S = H * P * H.transpose() + R + error_p;
        
        // 3. 计算卡尔曼增益 K
        K_curr = P * H.transpose() * S.inverse();

        // 4. 更新误差状态 (IEKF 公式)
        // x_{i+1} = x_prior + K * (residual + H * (x_curr - x_prior))
        // 对于线性系统 H，这项其实简化为 x_prior + K * (yk - H * x_prior)
        // 但保留完整形式以符合 ESIKF 结构
        
        Eigen::VectorXd innovation = residual + H * (x_curr - x_prior); 
        // 注意：标准 IEKF 推导中这里通常处理为 effective residual
        // 这里使用简化形式：x_new = x_prior + K * (y_meas - H * x_prior) 在线性下成立
        // 为了体现迭代过程，我们计算步长：
        
        Eigen::VectorXd step = residual - H * (x_prior - x_curr); // 这一项在线性下等于 yk - H * x_prior
        Eigen::VectorXd x_next = x_prior + K_curr * step;

        // 5. 收敛性检查
        if((x_next - x_curr).norm() < 1e-2) 
        {
            x_curr = x_next;
            cout<<"shoulian "<<endl;
            break;
        }
        x_curr = x_next;
    }

    // 6. 更新最终状态和协方差
    kk = K_curr;
    x = x_curr; 
    // 更新后验协方差 P = (I - KH)P
    P = (I - kk * H) * P; 
    update_time = time;
}

void KalmanFilter::predict(double time)
{
    // TicToc t_predict;
    // t_predict.tic();
    // cout<<"delta_time"<<time - ros::Time::now().toSec()<<endl;
    double dt = time - update_time;
    // cout<<"delta_t ======pre======= "<<dt<<endl;
    F.block(3,6,3,3) = -estimator.const_est_R*estimator.latest_R*get_hat(a_m);
    F.block(3,9,3,3) = -estimator.const_est_R*estimator.latest_R;
    F.block(6,6,3,3) = -get_hat(g_m);
    Fk = F * dt;
    // cout<<"Fk ======pre======= "<<F.block(3,6,3,3)<<endl;
    Phikk_1 = Eigen::MatrixXd::Identity(15, 15) + Fk;
    Phikk_1.block(6,6,3,3) = -get_hat(g_m)*dt;

    x = Phikk_1 * x;  //状态一步预测    X_(k/k-1) = F_(k/k-1)*X_(k-1)   第一个方程
    // cout<<"delta_vx ======pre======= "<<x(3)<<endl;
    P = Phikk_1 * P * Phikk_1.transpose() + Q; 

}

Eigen::VectorXd KalmanFilter::getPrediction()
{
    return x;
}

void KalmanFilter::keep(double time)
{
    double dt = time - update_time;
    if(dt<0.005)
    {
        dt = 0.005;
    }
    // cout<<"delta_t ======keep======= "<<dt<<endl;
    update_time = time;
    // F.block(3,6,3,3) = -estimator.const_est_R*get_hat(a_m);//
    F.block(3,6,3,3) = -estimator.const_est_R*estimator.latest_R*get_hat(a_m);//
    F.block(3,9,3,3) = -estimator.const_est_R*estimator.latest_R;//
    F.block(6,6,3,3) = -get_hat(g_m);
    Fk = F * dt;
    // cout<<"delta_vx ======keeep======= "<<x(3)<<endl;
    // cout<<"keep转移矩阵Fk"<<Fk<<endl;
    Phikk_1 = Eigen::MatrixXd::Identity(15, 15) + Fk;
    Phikk_1.block(6,6,3,3) = -get_hat(g_m)*dt;
    // P = Phikk_1 * P * Phikk_1.transpose() + Q;//这里出现问题？？？？？？？？？？S 
    // P = P - kk * H * P; //保持以下
    // cout<<"keep太转移矩阵Fk"<<Fk<<endl; 
}

Estimator::Estimator(): f_manager{Rs}
{
    ROS_INFO("init begins");
    initThreadFlag = false;
    plane_model_ = PlaneModel(1.0, 5);  // 内点阈值0.5米，最小内点数10
    max_time_diff_ = 3.0;                // 最大时间差1秒
    buf_capacity_ = 10; 
    clearState();
}

Estimator::~Estimator()
{
    if (MULTIPLE_THREAD)
    {
        processThread.join();
        printf("join thread \n");
    }
}

void Estimator::clearState()
{
    mProcess.lock();
    while(!accBuf.empty())
        accBuf.pop();
    while(!gyrBuf.empty())
        gyrBuf.pop();
    while(!featureBuf.empty())
        featureBuf.pop();

    prevTime = -1;
    curTime = 0;
    openExEstimation = 0;
    initP = Eigen::Vector3d(0, 0, 0);
    initR = Eigen::Matrix3d::Identity();
    inputImageCnt = 0;
    initFirstPoseFlag = false;

    get_ransac = false;

    for (int i = 0; i < WINDOW_SIZE + 1; i++)
    {
        Rs[i].setIdentity();
        Ps[i].setZero();
        Vs[i].setZero();
        Bas[i].setZero();
        Bgs[i].setZero();
        dt_buf[i].clear();
        linear_acceleration_buf[i].clear();
        angular_velocity_buf[i].clear();

        if (pre_integrations[i] != nullptr)
        {
            delete pre_integrations[i];
        }
        pre_integrations[i] = nullptr;
    }

    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        tic[i] = Vector3d::Zero();
        ric[i] = Matrix3d::Identity();
    }

    first_imu = false,
    sum_of_back = 0;
    sum_of_front = 0;
    frame_count = 0;
    solver_flag = INITIAL;
    initial_timestamp = 0;
    all_image_frame.clear();

    if (tmp_pre_integration != nullptr)
        delete tmp_pre_integration;
    if (last_marginalization_info != nullptr)
        delete last_marginalization_info;

    tmp_pre_integration = nullptr;
    last_marginalization_info = nullptr;
    last_marginalization_parameter_blocks.clear();

    f_manager.clearState();

    failure_occur = 0;

    mProcess.unlock();
}

void Estimator::setParameter()
{
    mProcess.lock();
    for (int i = 0; i < NUM_OF_CAM; i++)//NUM_OF_CAM==2
    {
        tic[i] = TIC[i];
        ric[i] = RIC[i];
        cout << " exitrinsic cam " << i << endl  << ric[i] << endl << tic[i].transpose() << endl;
    }
    f_manager.setRic(ric);
    ProjectionTwoFrameOneCamFactor::sqrt_info = FOCAL_LENGTH / 1.5 * Matrix2d::Identity();
    ProjectionTwoFrameTwoCamFactor::sqrt_info = FOCAL_LENGTH / 1.5 * Matrix2d::Identity();
    ProjectionOneFrameTwoCamFactor::sqrt_info = FOCAL_LENGTH / 1.5 * Matrix2d::Identity();
    td = TD;
    g = G;
    cout << "set g " << g.transpose() << endl;
    featureTracker.readIntrinsicParameter(CAM_NAMES);
    cov = -1.0;

    std::cout << "MULTIPLE_THREAD is " << MULTIPLE_THREAD << '\n';//值为1
    if (MULTIPLE_THREAD && !initThreadFlag)
    {
        initThreadFlag = true;
        processThread = std::thread(&Estimator::processMeasurements, this);
    }
    mProcess.unlock();
}

void Estimator::changeSensorType(int use_imu, int use_stereo)
{
    bool restart = false;
    mProcess.lock();
    if(!use_imu && !use_stereo)
        printf("at least use two sensors! \n");
    else
    {
        if(USE_IMU != use_imu)
        {
            USE_IMU = use_imu;
            if(USE_IMU)
            {
                // reuse imu; restart system
                restart = true;
            }
            else
            {
                if (last_marginalization_info != nullptr)
                    delete last_marginalization_info;

                tmp_pre_integration = nullptr;
                last_marginalization_info = nullptr;
                last_marginalization_parameter_blocks.clear();
            }
        }
        
        STEREO = use_stereo;
        printf("use imu %d use stereo %d\n", USE_IMU, STEREO);
    }
    mProcess.unlock();
    if(restart)
    {
        clearState();
        setParameter();
    }
}

void Estimator::inputImage(double t, const cv::Mat &_img, const cv::Mat &_img1)
{
    inputImageCnt++;
    map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> featureFrame;
    TicToc featureTrackerTime;

    if(_img1.empty())
        featureFrame = featureTracker.trackImage(t, _img);
    else
        featureFrame = featureTracker.trackImage(t, _img, _img1);//

    if (SHOW_TRACK)//SHOW_TRACK==1
    {
        cv::Mat imgTrack = featureTracker.getTrackImage();
        pubTrackImage(imgTrack, t);//
    }
    
    if(MULTIPLE_THREAD)  //MULTIPLE_THREAD==1
    {     
        if(inputImageCnt % 2 == 0)
        {
            mBuf.lock();
            featureBuf.push(make_pair(t, featureFrame));
            mBuf.unlock();
        }
    }
    else
    {
        mBuf.lock();
        featureBuf.push(make_pair(t, featureFrame));
        mBuf.unlock();
        TicToc processTime;
        processMeasurements();
        printf("process time: %f\n", processTime.toc());
    }
    
}

void Estimator::inputIMU(double t, const Vector3d &linearAcceleration, const Vector3d &angularVelocity)
{
    kf.a_m = linearAcceleration;
    kf.g_m = angularVelocity;
    mBuf.lock();
    accBuf.push(make_pair(t, linearAcceleration));
    gyrBuf.push(make_pair(t, angularVelocity));
    //printf("input imu with time %f \n", t);
    mBuf.unlock();

    if (solver_flag == NON_LINEAR)//
    {
        mPropagate.lock();//
        fastPredictIMU(t, linearAcceleration, angularVelocity);//
        pubLatestOdometry(latest_P, latest_Q, latest_V, t);//
        global_pubLatestOdometry(final_P, final_Q, final_V, t);//
        mPropagate.unlock();
    }
}

bool Estimator::optimize_R(std::vector<std::pair<Eigen::Vector3d, std::pair<double , Eigen::Vector3d>>> &point_vector,
                                    Eigen::Matrix3d &R,
                                    Eigen::Vector3d &t)
    {
        //---------------------------------------------------------------------------------------------------------
        if(!optimized)
        {
            if (point_vector.size() < 100)
            // if (point_vector.size() < 20)
            // if (point_vector.size() < 3)
                return false;
            R = Eigen::MatrixXd::Identity(3,3); 
            Eigen::Vector3d p_hat, q_hat;
            p_hat << 0, 0, 0;
            q_hat << 0, 0, 0;
            // SUM( R * VIO + t - ENU)
            for (auto &pair : point_vector)
            {
                // pair.second.z() = pair.first.z();
                p_hat += pair.first;//vio_t
                q_hat += pair.second.second;//gps_t
            }
            p_hat /= point_vector.size();
            q_hat /= point_vector.size();

            Eigen::MatrixXd X(3, point_vector.size());
            Eigen::MatrixXd Y(3, point_vector.size());
            Eigen::MatrixXd X_DOUBLE(2, point_vector.size());
            Eigen::MatrixXd Y_DOUBLE(2, point_vector.size());
            int i = 0;
            for (auto &pair : point_vector)//
            {
                // pair.second.z() = pair.first.z();
                X.col(i) = pair.first - p_hat;
                Y.col(i) = pair.second.second - q_hat;//
                i++;
            }
            X_DOUBLE = X.block(0,0,2,point_vector.size());
            Y_DOUBLE = Y.block(0,0,2,point_vector.size());
            Eigen::MatrixXd S = X_DOUBLE * Y_DOUBLE.transpose();//
            JacobiSVD<MatrixXd> svd(S, ComputeFullU | ComputeFullV);
            Eigen::MatrixXd U = svd.matrixU();
            Eigen::MatrixXd V = svd.matrixV();
            double det = (V * U.transpose()).determinant();
            int rows = V.cols();
            Eigen::MatrixXd Lambda = Eigen::MatrixXd::Identity(rows, rows);
            Lambda(rows - 1, rows - 1) = det;
            // Eigen::Matrix3d RR =  Eigen::Matrix3d::Identity();


            R.block<2,2>(0,0) = V * Lambda * U.transpose();
            R.bottomRightCorner<1,1>() = Eigen::Matrix<double, 1, 1>::Constant(1.0);

            t = q_hat - R * p_hat;//
            
            cout << "R" << endl;
            cout << R << endl;
            cout << "t" << endl;
            cout << t << endl;
            
            Eigen::Vector3d error ;
            for (auto &pair : point_vector)
            {
                pair.second.second.z() = pair.first.z();
                error += (R * pair.first + t - pair.second.second);
                // cout<<"error_th"<<error<<endl;
            }
                // error += (R * pair.first + t - pair.second).norm();
            error /= point_vector.size();
            cout << "error:"  <<  error << endl;//


            
            if (abs(R.determinant() - 1) < 1e-3)
                return true;//
            else
            {
                cout<<"tooooooooo_huge_error"<<endl;
                return false;
            }
        }
        else if(0)
        // else if(optimized)
        {
            if (point_vector.size() < 48)
            // if (point_vector.size() < 3)
                return false;
            R = Eigen::MatrixXd::Identity(3,3); 
            Eigen::Vector3d p_hat, q_hat;
            p_hat << 0, 0, 0;
            q_hat << 0, 0, 0;
            // SUM( R * VIO + t - ENU)
            for (auto &pair : point_vector)
            {
                // pair.second.z() = pair.first.z();
                p_hat += pair.first;//vio_t
                q_hat += pair.second.second;//gps_t
            }
            p_hat /= point_vector.size();
            q_hat /= point_vector.size();

            Eigen::MatrixXd X(3, point_vector.size());
            Eigen::MatrixXd Y(3, point_vector.size());
            Eigen::MatrixXd X_DOUBLE(2, point_vector.size());
            Eigen::MatrixXd Y_DOUBLE(2, point_vector.size());
            int i = 0;
            for (auto &pair : point_vector)//
            {
                // pair.second.z() = pair.first.z();
                X.col(i) = pair.first - p_hat;
                Y.col(i) = pair.second.second - q_hat;//
                i++;
            }
            X_DOUBLE = X.block(0,0,2,point_vector.size());
            Y_DOUBLE = Y.block(0,0,2,point_vector.size());
            Eigen::MatrixXd S = X_DOUBLE * Y_DOUBLE.transpose();//
            JacobiSVD<MatrixXd> svd(S, ComputeFullU | ComputeFullV);
            Eigen::MatrixXd U = svd.matrixU();
            Eigen::MatrixXd V = svd.matrixV();
            double det = (V * U.transpose()).determinant();
            int rows = V.cols();
            Eigen::MatrixXd Lambda = Eigen::MatrixXd::Identity(rows, rows);
            Lambda(rows - 1, rows - 1) = det;
            // Eigen::Matrix3d RR =  Eigen::Matrix3d::Identity();


            R.block<2,2>(0,0) = V * Lambda * U.transpose();
            R.bottomRightCorner<1,1>() = Eigen::Matrix<double, 1, 1>::Constant(1.0);

            t = q_hat - R * p_hat;//
            
            cout << "R" << endl;
            cout << R << endl;
            cout << "t" << endl;
            cout << t << endl;
            
            Eigen::Vector3d error ;
            for (auto &pair : point_vector)
            {
                pair.second.second.z() = pair.first.z();
                error += (R * pair.first + t - pair.second.second);
                // cout<<"error_th"<<error<<endl;
            }
                // error += (R * pair.first + t - pair.second).norm();
            error /= point_vector.size();
            cout << "error:"  <<  error << endl;//

            if (abs(R.determinant() - 1) < 1e-3)
                return true;
            else
            {
                cout<<"tooooooooo_huge_error"<<endl;
                return false;
            }
        }
    }

  bool Estimator::extractData(double local_time, std::map<double, pair<double,Eigen::Vector3d >> &datas, std::pair<double , Eigen::Vector3d> &output, double time_threshold)
  {
    try
    {
      if (datas.size() == 0)
        return false;
      std::map<double, pair<double,Eigen::Vector3d >>::iterator iter, data_iter;
      if (datas.find(local_time) != datas.end())//
      {
        output = datas.find(local_time)->second;
        return true;
      }
      else//到最后了
      {
        datas[local_time] = std::make_pair(0, Eigen::Vector3d(0, 0, 0));
        iter = datas.find(local_time);//
        assert(iter != datas.end());//
        if (iter == datas.begin())
          // ?because if datas.size() == 0 will continue. Thus it's allowed to add 1 to iter;
          data_iter = std::next(iter, 1);
        else if (iter == std::prev(datas.end(), 1))
          // ?because if datas.size() == 0 will continue. Thus it's allowed to add -1 to iter;
          data_iter = std::prev(iter, 1);
        else
          data_iter = (std::next(iter, 1)->first - local_time) < (local_time - std::prev(iter, 1)->first) ? std::next(iter, 1) : std::prev(iter, 1);//?

        if (std::abs(data_iter->first - local_time) < time_threshold)
        {
          datas.erase(datas.find(local_time));//
          output = data_iter->second;
          return true;
        }
        else
        {
          datas.erase(datas.find(local_time));
          printf("begin:%lf, end:%lf, local:%lf\n", datas.begin()->first, std::prev(datas.end(), 1)->first, local_time);
          printf("too big time_threshoold\n");
          return false;
        }
      }
    }
    catch (const std::exception &e)
    {
      std::cerr << e.what() << '\n';
    }//
    return false;
  }

Eigen::Matrix3d KalmanFilter::get_hat(Vector3d v)
  {
    Eigen::Matrix3d skew;
        skew << 0, -v(2), v(1),
            v(2), 0, -v(0),
            -v(1), v(0), 0;
    return skew;
  }

void Estimator::make_pair_to_get_R()
    {
        // cout<<"length_of_local = "<<local_odom_buf.size()<<"length_of_global"<<global_odom_buf.size()<<endl;
        if (!optimized)
        {
            // cout<<"start_to_optimize!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1"<<endl;
            // std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> point_vector;
            std::vector<std::pair<Eigen::Vector3d, std::pair<double , Eigen::Vector3d>>> point_vector;
            auto iter = local_odom_buf.begin();
            auto iterg = global_odom_buf.begin();
            while (iter != local_odom_buf.end())
            {
                double t = iter->first;
                Eigen::Vector3d lp = iter->second;
                // Eigen::Vector3d gp;
                std::pair<double , Eigen::Vector3d> gp;
                if (extractData(t, global_odom_buf, gp, 1.2))//
                {
                    point_vector.push_back(std::make_pair(lp, gp));
                }     
                iter++;
            }
            optimized = optimize_R(point_vector, est_R, est_t);
            const_est_R = est_R;
            const_est_t = est_t;
        }
        else
        {
            // if(optimized)//
            if(0)
            {
                // std::vector<std::pair<Eigen::Vector3d, Eigen::Vector3d>> point_vector;
                std::vector<std::pair<Eigen::Vector3d, std::pair<double , Eigen::Vector3d>>> point_vector;
                auto iter = local_odom_buf.begin();
                auto iterg = global_odom_buf.begin();
                while (iter != local_odom_buf.end())
                {
                    double t = iter->first;
                    Eigen::Vector3d lp = iter->second;
                    std::pair<double , Eigen::Vector3d> gp;
                    if (extractData(t, global_odom_buf, gp, 0.1))//
                        {
                            if(point_vector.size()<=60)
                            point_vector.push_back(std::make_pair(lp, gp));
                            // else if(point_vector.size()>200)
                        }     
                    iter++;
            }
            sortByDouble(point_vector);
            while(point_vector.size()>60)
            {

                point_vector.erase(--point_vector.end());
            }
            bool optimized_to_fix = optimize_R(point_vector, est_R, est_t);
            const_est_R = est_R;
            const_est_t = est_t;
            }
        }
    }

void Estimator::GPS2XYZ(double longitude, double latitude, double altitude, Eigen::Vector3d &xyz)
    {
        if(!initGPS)
        {
            geoConverter.Reset(latitude, longitude, altitude);//
            initGPS = true;
            // origin_height = altitude;
        }
        geoConverter.Forward(latitude, longitude, altitude, xyz(0), xyz(1), xyz(2));//
    }

Eigen::Vector3d Estimator::average(std::queue<pair<double,Eigen::Vector3d>> &gps_pre)
    {
        // cout<<"start_to_get_ave"<<endl;
        std::queue<pair<double,Eigen::Vector3d>> temp_gps_pre = gps_pre;

        Eigen::Vector3d sum(0,0,0);
        Eigen::Vector3d average_final;
        int len = gps_pre.size();
        // sum<<(0,0,0);
        //  cout<<sum<<endl;
        while (!temp_gps_pre.empty())
        {
        sum += temp_gps_pre.front().second;
        temp_gps_pre.pop();
        }
        average_final = sum/len;
        // cout<<"ave ="<<average_final<<endl;
        return average_final;
    }
double Estimator::variance(std::queue<pair<double,Eigen::Vector3d>> &gps_pre)
{
    std::queue<pair<double,Eigen::Vector3d>> temp_gps_pre = gps_pre;
    Eigen::Vector3d sum_variance(0,0,0);
    Eigen::Vector3d average_temp(0,0,0);
    Eigen::Vector3d seigama_temp(0,0,0);
    int len = gps_pre.size();
    // cout<<"length_variance"<<len<<endl;
    average_temp = average(gps_pre);
    while (!temp_gps_pre.empty())
    {
        sum_variance.x() += pow(temp_gps_pre.front().second.x() - average_temp.x() , 2);
        sum_variance.y() += pow(temp_gps_pre.front().second.y() - average_temp.y() , 2);
        sum_variance.z() += pow(temp_gps_pre.front().second.z() - average_temp.z() , 2);
        temp_gps_pre.pop();
    }
    seigama_temp.x() = sum_variance.x()/len;
    seigama_temp.y() = sum_variance.y()/len;
    seigama_temp.z() = sum_variance.z()/len;
    // cout<<"variance"<<seigama_temp<<endl;
    return seigama_temp.norm();
}

Eigen::Quaterniond Estimator::euler2quaternion(Vector3d euler)//
  {
    double cr = cos(euler(0)/2);//x
    double sr = sin(euler(0)/2);
    double cp = cos(euler(1)/2);//y
    double sp = sin(euler(1)/2);
    double cy = cos(euler(2)/2);//z
    double sy = sin(euler(2)/2);
    // double cr = cos(euler(2)/2);
    // double sr = sin(euler(2)/2);
    // double cp = cos(euler(1)/2);
    // double sp = sin(euler(1)/2);
    // double cy = cos(euler(0)/2);
    // double sy = sin(euler(0)/2);
    Quaterniond q;
    q.w() = cr*cp*cy + sr*sp*sy;
    q.x() = sr*cp*cy - cr*sp*sy;
    q.y() = cr*sp*cy + sr*cp*sy;
    q.z() = cr*cp*sy - sr*sp*cy;
    return q; 
  }

double Estimator::get_stimulate(std::queue<pair<double,Eigen::Vector3d>> &a, std::queue<pair<double,Eigen::Vector3d>> &g)
// double Estimator::get_stimulate(std::queue<pair<double,Eigen::Vector3d>> &a)
    {
        // cout<<"length_of_acc&gyr+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++s"<<a.size()<<endl;
        Eigen::Vector3d a_mean = average(a);
        // Eigen::Vector3d g_mean = average(g);
        // cout<<"a_mean"<<a_mean<<endl;
        double a_var = variance(a);
        // cout<<"a_var"<<a_var<<endl;
        double g_var = variance(g);
        double length_of_stimulate = a.size();
        // cout<<"length_of_stimulate"<<length_of_stimulate<<endl;
        std::queue<pair<double,Eigen::Vector3d>> a_temp = a;
        std::queue<pair<double,Eigen::Vector3d>> g_temp = g;
        double sum = 0;
        std::vector<Eigen::Vector3d>::iterator itor_a , itor_g;
        // int i = 0;
        while (!a_temp.empty())
        {
            while(!g_temp.empty())
            {
                sum+=pow(g_temp.front().second.norm(),2)/g_var;
                g_temp.pop();
            }
            sum+=pow((a_temp.front().second - 9.805*a_mean/a_mean.norm()).norm(),2)/a_var;
            a_temp.pop();
            // i++;
            // cout<<"num"<<i<<endl;
        }
        // cout<<"stimulate = "<<sum<<endl;
        sum = sum/length_of_stimulate;
        // cout<<"stimulate = "<<sum<<endl;
        return sum;
    }
//---------------------------------------------------------------------------------------------------------------------

void Estimator::initializePlaneModel(const std::vector<std::pair<double, Eigen::Vector3d>>& buffer) {

    int num_samples = std::min(plane_model_.min_inliers, (int)buffer.size());
    auto samples = randomSample(buffer, num_samples);

    for (int i = 0; i < samples.size() - 2; ++i) {
        Eigen::Vector3d p1 = samples[i].second;
        Eigen::Vector3d p2 = samples[i+1].second;
        Eigen::Vector3d p3 = samples[i+2].second;


        Eigen::Vector3d v1 = p2 - p1;
        Eigen::Vector3d v2 = p3 - p1;
        if (v1.cross(v2).norm() < 1e-6) continue; 


        Eigen::Vector3d normal = v1.cross(v2).normalized();
        plane_model_.a = normal.x();
        plane_model_.b = normal.y();
        plane_model_.c = normal.z();
        plane_model_.d = -(normal.dot(p1));  


        std::vector<std::pair<double, Eigen::Vector3d>> temp_inliers;
        for (const auto& data : buffer) {
            if (distanceToPlane(data.second) < plane_model_.inlier_thresh) {
                temp_inliers.push_back(data);
            }
        }


        if (temp_inliers.size() >= plane_model_.min_inliers) {
            inliers_ = temp_inliers;
            return;
        }
    }


    inliers_.clear();
}

void Estimator::updateInliersWithCurrentModel(const std::vector<std::pair<double, Eigen::Vector3d>>& buffer) {
    std::vector<std::pair<double, Eigen::Vector3d>> new_inliers;
    for (const auto& data : buffer) {
        if (distanceToPlane(data.second) < plane_model_.inlier_thresh) {
            cout<<"threshold =========="<<distanceToPlane(data.second)<<endl;
            new_inliers.push_back(data);
        }
    }
    inliers_.swap(new_inliers); 
}

void Estimator::reestimatePlaneModelIfNeeded(const std::vector<std::pair<double, Eigen::Vector3d>>& buffer) {
    if (inliers_.size() >= plane_model_.min_inliers * 2) {  
        Eigen::Vector3d centroid = computeCentroid(inliers_);
        Eigen::Matrix3d A = Eigen::Matrix3d::Zero();
        for (const auto& data : inliers_) {
            Eigen::Vector3d p = data.second - centroid;
            A += p * p.transpose();  // 
        }
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> solver(A);
        Eigen::Vector3d normal = solver.eigenvectors().col(0);  // 
        
        plane_model_.a = normal.x();
        plane_model_.b = normal.y();
        plane_model_.c = normal.z();
        plane_model_.d = -(normal.dot(centroid));  // 

        // 重新筛选内点（基于新模型）
        std::vector<std::pair<double, Eigen::Vector3d>> refined_inliers;
        for (const auto& data : buffer) {
            if (distanceToPlane(data.second) < plane_model_.inlier_thresh) {
                cout<<"threshold =========="<<distanceToPlane(data.second)<<endl;
                refined_inliers.push_back(data);
            }
        }
        inliers_.swap(refined_inliers);  //
    }
}

double Estimator::distanceToPlane(const Eigen::Vector3d& point) {
    return std::abs(plane_model_.a * point.x() + 
                   plane_model_.b * point.y() + 
                   plane_model_.c * point.z() + 
                   plane_model_.d) / 
           std::sqrt(plane_model_.a*plane_model_.a + 
                     plane_model_.b*plane_model_.b + 
                     plane_model_.c*plane_model_.c);
}

Eigen::Vector3d Estimator::computeCentroid(const std::vector<std::pair<double, Eigen::Vector3d>>& points) {
    Eigen::Vector3d centroid = Eigen::Vector3d::Zero();
    for (const auto& p : points) {
        centroid += p.second;
    }
    return centroid / points.size();
}

std::vector<std::pair<double, Eigen::Vector3d>> Estimator::randomSample(
    const std::vector<std::pair<double, Eigen::Vector3d>>& buffer, int n) {
    std::vector<std::pair<double, Eigen::Vector3d>> samples;
    std::uniform_int_distribution<size_t> dist(0, buffer.size() - 1);
    for (int i = 0; i < n; ++i) {
        size_t idx = dist(plane_model_.rng);
        samples.push_back(buffer[idx]);
    }
    return samples;
}

void Estimator::removeOutdatedData(double current_time) {
    while (!pos_buffer_.empty() && 
           (current_time - pos_buffer_.front().first) > max_time_diff_) {
        cout<<"too_huge_time????"<<endl;
        pos_buffer_.erase(pos_buffer_.begin());
    }
}

void Estimator::maintainBufferSize(std::vector<std::pair<double, Eigen::Vector3d>>& buffer, 
                                 const std::pair<double, Eigen::Vector3d>& new_data) {
    buffer.push_back(new_data);
    if (buffer.size() > buf_capacity_) {
        buffer.erase(buffer.begin());  // 
    }
}
Eigen::Vector3d Estimator::processGPS(const double& timestamp, const Eigen::Vector3d& gps_pos) {
    std::pair<double, Eigen::Vector3d> new_data = {timestamp, gps_pos};
    maintainBufferSize(pos_buffer_, new_data);

    removeOutdatedData(timestamp);

    if (pos_buffer_.size() >= plane_model_.min_inliers) {
        updatePlaneModel(pos_buffer_);
    }

    Eigen::Vector3d latest_inlier = getLatestInlierPosition();
    if (!latest_inlier.isZero()) {

        return latest_inlier;
    }
}

Eigen::Vector3d Estimator::getLatestInlierPosition() {
    if (inliers_.empty()) return Eigen::Vector3d::Zero();

    // 查找时间最新的内点
    const auto* latest = &inliers_[0];
    for (const auto& data : inliers_) {
        if (data.first > latest->first) {
            latest = &data;
        }
    }
    return latest->second;
}

void Estimator::updatePlaneModel(const std::vector<std::pair<double, Eigen::Vector3d>>& buffer) {
    if (buffer.size() < plane_model_.min_inliers) {
        inliers_.clear();
        return;
    }

    if (inliers_.empty()) {
        initializePlaneModel(buffer);
        return;
    }

    updateInliersWithCurrentModel(buffer);
    reestimatePlaneModelIfNeeded(buffer);
}
//11-----------------------------------------------------------------------------------------------------------------
void Estimator::inputGPS(double t, const Eigen::Vector3d &lon_lat_alt_temp, const Eigen::MatrixXd &conv)
{
    // time_of_gps_used = ros::Time::now().toSec();
    // cout<<"22222222222222222222222222222222-vio"<<endl;
    // TicToc timer9;
    // timer9.tic();
    int count_nei = 0;
    int total_nei = 0;
    inliners_points.clear();
    bool gps_pretreatment = false;
    Eigen::Vector3d global_odom, global_odom_before;
    
    double t_global = t;
    double t_buf = ros::Time::now().toSec();
    int init_gps_cao = 0;
    predice_or_not = conv(0,0);
    cov = conv(0,0);
    // cout<<conv(0,0)<<endl;
    if(initGPS)
    {
        GPS2XYZ(lon_lat_alt_temp.x(), lon_lat_alt_temp.y(), Ps[WINDOW_SIZE].z(), global_odom_before);
        // GPS2XYZ(lon_lat_alt_temp.x(), lon_lat_alt_temp.y(), Ps[WINDOW_SIZE].z(), global_odom);
        init_gps_cao = 2;
    }
    else if(!initGPS)
    {
        GPS2XYZ(lon_lat_alt_temp.x(), lon_lat_alt_temp.y(), Ps[WINDOW_SIZE].z(), global_odom_before);//
        // GPS2XYZ(lon_lat_alt_temp.x(), lon_lat_alt_temp.y(), Ps[WINDOW_SIZE].z(), global_odom);
        init_gps_cao = 1;//
    }
    // cout << "Elapsed time: " << elapsed_time << " seconds一切正常" << endl;
    // cout<<"init_gps_cao = "<<init_gps_cao<<endl;
    // if(global_odom_buf.size() == 0)
    //222222222222222222222222222222222222222222222222222222222222222222222222222222222222025.6.11
    // if(has_been_indoors)
    // {
    //     pos_buffer_.clear();
    // }
    // global_odom = processGPS(t_buf, global_odom_before);
    

    //22222222222222222222222222222222222222222222222222222222222222222222222222222222222222222222
    //1111111111111111111111111111111111111111111111111111111111111111111111111111111112025.2.26
    global_odom_before_buf.push_back(std::make_pair(t_buf,global_odom_before));
    while(global_odom_before_buf.size()>10)
    {
        global_odom_before_buf.erase(global_odom_before_buf.begin());
        // cout<<"length_of_gps_before"<<gps_before.size()<<endl;//没有进来过？
    }
    // if((global_odom_before_buf.size()==10))
    if((global_odom_before_buf.size()==10))
    {
        if(!get_ransac)
        {
            for(int i=0;i<8;i++)
            {
                int first_shot=rand()%8;//产生0-10的随机数
                int second_shot=rand()%8;//产生0-10的随机数

                std::pair<double, Eigen::Vector3d> sample_first , best_sample_first;
                std::pair<double, Eigen::Vector3d> sample_second , best_sample_second;

                sample_first = global_odom_before_buf[first_shot];
                sample_second = global_odom_before_buf[second_shot];

                if(first_shot == second_shot)
                {
                    continue;
                }

                k_estimate = (sample_second.second.y() - sample_first.second.y())/(sample_second.second.x() - sample_first.second.x());
                b_estimate = sample_first.second.y() - k_estimate*sample_first.second.x();

                for (int i = 0; i < global_odom_before_buf.size(); i++)
                {
                    // delta = k * x  + b - y
                    double y_error = abs(k_estimate * global_odom_before_buf[i].second.x() + b_estimate - global_odom_before_buf[i].second.y()) / sqrt((pow(k_estimate, 2) + 1));
                    // if(y_error<=0.001)
                    // {
                    //     y_error = 0.02;
                    // }
                    // cout << "y_error = " << y_error << endl;
                    if (y_error < threshold)
                    {
                        count_nei++;
                    }
                }
                if (count_nei > total_nei)
                {
                    total_nei = count_nei;
        
                    best_sample_first = sample_first;
                    best_sample_second = sample_second;
        
                    best_parameters_k = k_estimate;
                    best_parameters_b = b_estimate;
                }
                count_nei = 0;
            }
        }
        
        std::vector<std::pair<double, Eigen::Vector3d>> inliner_points;

        if(!get_ransac)
        {
            inline_ransac = 0;
            all_in_ransac = 0;
            for (int i = 0; i < global_odom_before_buf.size(); ++i)
            {
                all_in_ransac++;
                // delta = k * x  + b - y
                double y_error = abs(best_parameters_k * global_odom_before_buf[i].second.x() + best_parameters_b - global_odom_before_buf[i].second.y()) / sqrt((pow(best_parameters_k, 2) + 1));
                cout << "y_error = " << y_error << endl;
                if (y_error < threshold)
                {
                    // cout<<"find_inliner"<<i<<endl;
                    inliners_points.push_back(global_odom_before_buf[i]);
                    inline_ransac++;
                    // cout<<"gps_before"<<i<<" = "<<gps_before[i].second<<endl;//这个值一切正常
                }
            }
            // get_ransac = true;//始终将其置为false，即增量式ransac算法的不进行增量，认为每次判定内点数量均达不到阈值，每次都会重新建立ransac拟合模型，这种情况下，程序正常运行。
        }
        else
        {
            double y_error_end = abs(best_parameters_k * global_odom_before_buf.back().second.x() + best_parameters_b - global_odom_before_buf.back().second.y()) / sqrt((pow(best_parameters_k, 2) + 1));
            double y_error_begin = abs(best_parameters_k * global_odom_before_buf[0].second.x() + best_parameters_b - global_odom_before_buf[0].second.y()) / sqrt((pow(best_parameters_k, 2) + 1));
            cout<<"end = "<<y_error_end<<endl;
            if(y_error_end <=  threshold)
            {
                inliners_points.push_back(global_odom_before_buf.back());
                inline_ransac++;
            }
            else{get_ransac = false;}
            inline_ransac--;
            cout<<"rate = "<<(double)inline_ransac/(double)all_in_ransac<<endl;
            if((double)inline_ransac/(double)all_in_ransac < 0.9)
            {
                get_ransac = false;
            }
        }
        while(inliners_points.size()>10)
        {
            inliners_points.erase(inliners_points.begin());
        }
        int len_of_gps_before = inliners_points.size();
        // cout<<"after_ransac"<<inliners_points.size()<<endl;
        // global_odom = inliners_points[len_of_gps_before-1].second;
        global_odom = inliners_points.back().second;
        gps_pretreatment = true;
    }
    //11111111111111111111111111111111111111111111111111111111111111111111111111111111112025.2.26

    if(init_gps_cao == 1)
    {
        if(conv(0,0)<=CONV_THRESHOLD&&!global_odom.isZero())//0.04,0.6?0.5?
        {
            global_odom_init = global_odom;//
            // gnss_time_offset = (ros::Time::now() - GPS_msg->header.stamp).toSec();
        }
    }
    // cout<<"init_gps_pose"<<global_odom_init<<endl;
    if(conv(0,0)<=CONV_THRESHOLD&&!global_odom.isZero())//0.04-0.6-0.5---------------------------第一个大if----------------------------------------
    // if(conv(0,0)<=CONV_THRESHOLD)
    {
        // cout<<"22222222222222222222222222222222-vio"<<endl;
        gp = global_odom - global_odom_init;
        gp_pass.push_back(gp);
        // cout<<"gp = "<<gp<<endl;
        // if(delta_p.norm()>0.8)
        // {
        //     kf.error_p << delta_p.x(),0,0,
        //                   0,delta_p.y(),0,
        //                   0,0,delta_p.z();
        // }
        // else if(delta_p.norm()<0.8)
        // {
        //     kf.error_p = 0.000001*Eigen::Matrix3d::Identity();
        // }

        if(!optimized)
        {
            // cout<<"333333333333333333333333333333333333-vio"<<endl;
            norm = 0.1;//初始化一个值试试？？？？

            if(accBuf.empty()&&gyrBuf.empty())
            {
                printf("not receive imu\n");
                return;
            }
            stimulate = get_stimulate(accBuf , gyrBuf);
            cout<<"stimulate =================================== "<<stimulate<<endl;
            // stimulate = get_stimulate(accBuf);//std::queue<pair<double, Eigen::Vector3d>>
            if(stimulate>0.6)
            {
                stimulating  = true;
            }
            else stimulating = false;
            if(stimulating)
            {
                global_odom_buf[t_buf].first = conv(0,0);
                global_odom_buf[t_buf].second = gp;//
            }
        }
        if(optimized)
        {
            // cout<<"4444444444444444444444444444444444444-vio"<<endl;
            global_odom_buf[t_buf].first = conv(0,0);
            global_odom_buf[t_buf].second = gp;//
            if (global_odom_buf.size() > 400) 
            {
                // 
                std::map<double, pair<double,Eigen::Vector3d >>::iterator it = global_odom_buf.begin();
                while (global_odom_buf.size() > 400) 
                {
                    global_odom_buf.erase(it); // 
                }
            }
            // cout<<"delta_t_of_gps&zitai"<<t_buf - temp_t_l_to_g_p.first<<endl;
            Eigen::Vector3d delta_p = final_P - gp;
            // cout<<"final_p = "<<final_P<<"gp = "<<gp<<endl;
            // cout<<"33333333333333333333333333333333333333333333333333"<<delta_p.norm()<<endl;
            Eigen::Vector3d yk = l_to_g_p - gp;// l_to_g_p在localcallback里面
            // cout<<"yk1 = "<<yk<<endl;
            yk.z() = 0.0;


            

            if(has_been_indoors)
            {
                if(delta_p.norm()<0.4)
                {
                    times_in_threshold++;
                    if(times_in_threshold>5)//15
                    {
                        kf.error_p = 0*Eigen::Matrix3d::Identity();
                        has_been_indoors = false;
                        has_been_in = false;
                        norm = 0.1;
                        cout<<"allright--------------------------------------------------"<<endl;
                    }
                }
                else if(delta_p.norm()>=0.4)
                {
                    times_in_threshold = 0;

                    double t_lost = 0.0;
                    if(last_valid_gnss_time > 0.01)
                    {
                        t_lost = t_global - last_valid_gnss_time;
                        last_valid_gnss_time = -1;
                    }

                    kf.error_p << abs(delta_p.x())/delta_p.norm(),0,0,
                                0,abs(delta_p.y())/delta_p.norm(),0,
                                0,0,abs(delta_p.z())/delta_p.norm();
                    cout<<"t_lost = "<<t_lost<<endl;
                    kf.error_p = std::log10(t_lost+1.0) * 1.0 * kf.error_p;
                    t_lost*=0.5;

                    // kf.error_p << delta_p.x(),0,0,
                    //             0,delta_p.y(),0,
                    //             0,0,delta_p.z();
                    cout<<"delta_p"<<delta_p<<endl;
                    if(gp_pass.size()>2)
                    {
                        Eigen::Vector3d end = gp_pass.back();
                        Eigen::Vector3d before_end = gp_pass[gp_pass.size()-2];
                        norm = (end - before_end).norm();
                        cout<<"norm_="<<norm<<endl;
                        if(norm > 1.2)//1.0
                        {
                            norm = 0.15;
                            gp_pass.erase(gp_pass.begin(), gp_pass.end() - 2);
                        }
                        else
                        {
                            norm *= 1.4;
                        }
                    }
                    yk = yk_still+norm*((yk-yk_still)/(yk-yk_still).norm());
                    cout<<"yk_when_error = "<<yk<<endl;
                    yk_still = yk;
                    // cout<<"yk2 = "<<yk<<endl;
                }
            }

            
            //用GPS来的时间试试？？？
            if(kf.update_time == -1)
            {
                kf.initState(yk,const_est_R,t_global);
                kf.error_p.setZero();
                
            }
            else
            {
                // cout<<"888888888888888888888888888888888888888888-vio"<<endl;
                Eigen::VectorXd state_in_last_phase;
                Eigen::MatrixXd Fk_in_last_phase;
                Eigen::MatrixXd P_in_last_phase;

                
                
                if(conv(0,0) < CONV_THRESHOLD&&!global_odom.isZero())//0.04,0.6,0.5
                {
                    // cout<<"777777777777777777777777777777777777-vio"<<endl;
                    // cout<<"bool = "<<optimize_twice<<endl;
                    
                    // timer9.tic();
                    state_in_last_phase = kf.getPrediction();
                    // double elapsed_time9 = timer9.toc(); // 
                    // state_in_last_phase = kf.x;
                    // cout<<"state_before_update-------------------------------------------"<<state_in_last_phase<<endl;//
                    Fk_in_last_phase = kf.Fk;
                    P_in_last_phase = kf.P;
                    // cout<<"ykx = "<<yk<<endl;
                    // TicToc timer9;
                    // timer9.tic();
                    
                    // kf.update(yk, t_global, 1,conv);
                    kf.updateESIKF(yk, t_global, 1, conv);
                    // cout<<"kf.x22222222222222222222222222222222222222222222222222222222222222"<<kf.x<<endl;
                    state = kf.getPrediction();
                    // double elapsed_time9 = timer9.toc(); // 
                    // state(2) = global_z.z();
                    Eigen::Vector3d position_after_filter;
                    position_after_filter.x() = l_to_g_p.x()-state(0);
                    position_after_filter.y() = l_to_g_p.y()-state(1);
                    position_after_filter.z() = l_to_g_p.z()-state(2);
                    // cout<<"state良好 = "<<state(3)<<endl;
                    // cout<<"assert_maybe良好 = "<<position_after_filter<<endl;
                    // if(abs(latest_P.z())>0.3)//再考虑用什么标志位？
                    
                    // if(get_Q&&get_P)
                    // {
                    //     // cout<<"boat_new_threshold"<<(latest_P - position_after_filter).norm()<<endl;
                    //     if(has_been_indoors&&((final_P - position_after_filter).norm() > 0.8))//latest_P这里是vio的结果，————检查
                    //     {
                    //         delta_P = final_P - position_after_filter;
                    //         cout<<"error_jump"<<delta_P<<endl;
                    //         optimize_twice =  true;
                    //         times_in_threshold = 0;
                    //         cout<<"wo_kan_ni_you_wen_ticccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc"<<endl;
                    //         cout<<"pat = "<<pat<<endl;
                    //     }
                    //     else if(has_been_indoors&&((final_P - position_after_filter).norm() <= 0.8))
                    //     {
                    //         times_in_threshold ++;
                    //         optimize_twice = false;
                    //         cout<<"better"<<endl;
                    //         // const_est_R = est_R;
                    //     }

                    //     if(times_in_threshold > 3)
                    //     {
                    //         times_in_threshold = 0;
                    //         pat = 1;
                    //         ready_to_predict = true;
                    //         has_been_indoors = false;
                    //         we_loss_it = false;
                    //         kf.m_update_time = -1;
                    //     }
                    // }
                    // if(!optimize_twice)//如果正常
                    // {
                        pubodom_raw[t_global] = position_after_filter;//
                        // cout<<"final_P000000000000000000000000000000000000"<<position_after_filter<<endl;
                        Eigen::Vector3d rotate_after_filter(state(6),state(7),state(8));
                        delta_q_to_fix = euler2quaternion(rotate_after_filter);
                        
                        if(pubodom_raw.size()>6)
                        {
                            double i = -1.0;
                            double total = 0;
                            Eigen::Vector3d sum_temp(0,0,0);
                            std::map<double, Eigen::Vector3d>::iterator iter = pubodom_raw.begin(); 
                            while(iter!=pubodom_raw.end())
                            {
                                // cout<<"final_P22222222222222222222222222222222222="<<iter->second<<endl;
                                // cout<<"final_P33333333333333333333333333333333333="<<i<<endl;
                                // cout<<"final_P44444444444444444444444444444444444="<<sum_temp<<endl;
                                sum_temp.x() += i*iter->second.x();
                                sum_temp.y() += i*iter->second.y();
                                sum_temp.z() += i*iter->second.z();
                                total +=i;
                                // cout<<"final_P55555555555555555555555555555555555="<<total<<endl;
                                i += 1.0;
                                iter++;
                            }
                            // cout<<"total = "<<total<<endl;
                            final_P = sum_temp/total;
                            // cout<<"final_P111111111111111111111111"<<final_P<<endl;
                            final_V.x() = V_n.x() -state(3);
                            final_V.y() = V_n.y() -state(4);
                            final_V.z() = V_n.z() -state(5);
                            final_Q = delta_q_to_fix*Q_n_kf;
                            get_Q = true;
                            get_P = true;
                            // cout<<"P&Q——————————————————————————————GOT"<<endl;
                            
                        }
                        while(pubodom_raw.size()>6)
                        {
                            pubodom_raw.erase(pubodom_raw.begin());
                            // cout<<"length_of_pubodom_raw"<<pubodom_raw.size()<<endl;
                        }
                }
                else 
                {
                    Eigen::Vector3d position_after_filter;
                    position_after_filter.x() = l_to_g_p.x()-state(0);
                    position_after_filter.y() = l_to_g_p.y()-state(1);
                    position_after_filter.z() = l_to_g_p.z()-state(2);
                    // cout << "Elapsed time: " << elapsed_time << " seconds一切正常3" << endl;
                    // cout<<"assert_maybe奇辣 = "<<position_after_filter<<endl;
                    pubodom_raw[t_global] = position_after_filter;
                    Eigen::Vector3d rotate_after_filter(state(6),state(7),state(8));
                    delta_q_to_fix = euler2quaternion(rotate_after_filter);
                    if(pubodom_raw.size()>6)//6
                        {
                            double i = -1.0;
                            double total = 0;
                            Eigen::Vector3d sum_temp(0,0,0);
                            std::map<double, Eigen::Vector3d>::iterator iter = pubodom_raw.begin(); 
                            while(iter!=pubodom_raw.end())
                            {
                                // cout<<"final_P22222222222222222222222222222222222="<<iter->second<<endl;
                                // cout<<"final_P33333333333333333333333333333333333="<<i<<endl;
                                // cout<<"final_P44444444444444444444444444444444444="<<sum_temp<<endl;
                                sum_temp.x() += i*iter->second.x();
                                sum_temp.y() += i*iter->second.y();
                                sum_temp.z() += i*iter->second.z();
                                total +=i;
                                // cout<<"final_P55555555555555555555555555555555555="<<total<<endl;
                                i += 1.0;
                                iter++;
                            }
                            // cout<<"total = "<<total<<endl;
                            final_P = sum_temp/total;
                            final_V.x() = V_n.x() -state(3);
                            final_V.y() = V_n.y() -state(4);
                            final_V.z() = V_n.z() -state(5);
                            final_Q = delta_q_to_fix*Q_n_kf;
                            get_Q = true;
                            get_P = true;
                            // double elapsed_time9 = timer9.toc(); // 
                            // cout << "Elapsed time: " << elapsed_time9 << " seconds一切正常1" << endl;
                        }
                        while(pubodom_raw.size()>6)
                        {
                            pubodom_raw.erase(pubodom_raw.begin());
                            // cout<<"length_of_pubodom_raw"<<pubodom_raw.size()<<endl;
                        }
                } 
            }
        }
    }
    else if(conv(0,0)>CONV_THRESHOLD&&!global_odom.isZero())//0.04 0.5 0.6
    {

        

        pos_buffer_.clear();
        
        if(optimized)
        {
            if(!has_been_indoors)
            {
                last_valid_gnss_time = t_global; // 初始化时间
            }
            has_been_indoors = true;
            if(!has_been_in)
            {
                yk_still = l_to_g_p-gp;
                has_been_in = true;
            }
            
            Eigen::Vector3d position_after_filter;
            position_after_filter.x() = l_to_g_p.x()-state(0);
            position_after_filter.y() = l_to_g_p.y()-state(1);
            position_after_filter.z() = l_to_g_p.z()-state(2);
            // cout<<"state吉拉2 = "<<state(3)<<endl;
            // cout<<"assert_maybe奇辣2 = "<<position_after_filter<<endl;
            // cout << "Elapsed time: " << elapsed_time << " seconds一切正常" << endl;
            pubodom_raw[t_global] = position_after_filter;
            Eigen::Vector3d rotate_after_filter(state(6),state(7),state(8));
            delta_q_to_fix = euler2quaternion(rotate_after_filter);

                // if(conv(0,0)<1.5)
                // {
                //     cout<<"conv(0,0)<1.5"<<endl;
                //     if(kf.m_update_time==-1)
                //     {
                //         kf.mem_init_state(position_after_filter , t_global);
                //         cout<<"ready_to_loss"<<endl;
                //     }
                //     else
                //     {
                //         Eigen::Vector3d g_xyz = global_odom - global_odom_init;
                //         kf.mem_update(g_xyz,t_global,1,conv);
                //         if(we_loss_it)
                //         {
                //             Eigen::Vector3d temp_g_p;
                //             temp_g_p.x() = kf.m_x(0);
                //             temp_g_p.y() = kf.m_x(1);
                //             temp_g_p.z() = kf.m_x(2);
                //             pubodom_raw[t_global] = temp_g_p;
                //         }
                //     }
                // }
                // if(conv(0,0)>=1.5)
                // {
                //     we_loss_it = true;
                // }

                if(pubodom_raw.size()>6)
                {
                    double i = -1.0;
                    double total = 0;
                    Eigen::Vector3d sum_temp(0,0,0);
                    std::map<double, Eigen::Vector3d>::iterator iter = pubodom_raw.begin(); 
                    while(iter!=pubodom_raw.end())
                    {
                        // cout<<"final_P22222222222222222222222222222222222="<<iter->second<<endl;
                        // cout<<"final_P33333333333333333333333333333333333="<<i<<endl;
                        // cout<<"final_P44444444444444444444444444444444444="<<sum_temp<<endl;
                        sum_temp.x() += i*iter->second.x();
                        sum_temp.y() += i*iter->second.y();
                        sum_temp.z() += i*iter->second.z();
                        total +=i;
                        // cout<<"final_P55555555555555555555555555555555555="<<total<<endl;
                        i += 1.0;
                        iter++;
                    }
                    // cout<<"total_="<<total<<endl;

                    final_P = sum_temp/total;

                    final_V.x() = V_n.x() -state(3);
                    final_V.y() = V_n.y() -state(4);
                    final_V.z() = V_n.z() -state(5);
                    final_Q = delta_q_to_fix*Q_n_kf;
                    get_Q = true;
                    get_P = true;
                    while(pubodom_raw.size()>6)
                    {
                        pubodom_raw.erase(pubodom_raw.begin());
                        // cout<<"length_of_pubodom_raw"<<pubodom_raw.size()<<endl;
                    }
                    // double elapsed_time9 = timer9.toc(); // 
                    // cout << "Elapsed time: " << elapsed_time9 << " seconds一切正常1" << endl;
                }   
        }
    }
}

void Estimator::inputFeature(double t, const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> &featureFrame)
{
    mBuf.lock();
    featureBuf.push(make_pair(t, featureFrame));//
    mBuf.unlock();

    if(!MULTIPLE_THREAD)//MULTIPLE_THREAD==1
        processMeasurements();//
}


bool Estimator::getIMUInterval(double t0, double t1, vector<pair<double, Eigen::Vector3d>> &accVector, 
                                vector<pair<double, Eigen::Vector3d>> &gyrVector)
{
    if(accBuf.empty())
    {
        printf("not receive imu\n");
        return false;
    }
    //printf("get imu from %f %f\n", t0, t1);
    //printf("imu fornt time %f   imu end time %f\n", accBuf.front().first, accBuf.back().first);
    if(t1 <= accBuf.back().first)
    {
        while (accBuf.front().first <= t0)
        {
            accBuf.pop();
            gyrBuf.pop();
        }
        while (accBuf.front().first < t1)
        {
            accVector.push_back(accBuf.front());
            accBuf.pop();
            gyrVector.push_back(gyrBuf.front());
            gyrBuf.pop();
        }
        accVector.push_back(accBuf.front());
        gyrVector.push_back(gyrBuf.front());
    }
    else
    {
        printf("wait for imu\n");
        return false;
    }
    return true;
}

bool Estimator::IMUAvailable(double t)
{
    if(!accBuf.empty() && t <= accBuf.back().first)
        return true;
    else
        return false;
}

void Estimator::processMeasurements()
{
    while (1)
    {
        //printf("process measurments\n");
        pair<double, map<int, vector<pair<int, Eigen::Matrix<double, 7, 1> > > > > feature;
        vector<pair<double, Eigen::Vector3d>> accVector, gyrVector;
        if(!featureBuf.empty())
        {
            feature = featureBuf.front();
            curTime = feature.first + td;
            while(1)
            {
                if ((!USE_IMU  || IMUAvailable(feature.first + td)))
                    break;
                else
                {
                    printf("wait for imu ... \n");
                    if (! MULTIPLE_THREAD)
                        return;
                    std::chrono::milliseconds dura(5);
                    std::this_thread::sleep_for(dura);
                }
            }
            mBuf.lock();
            if(USE_IMU)
                getIMUInterval(prevTime, curTime, accVector, gyrVector);//

            featureBuf.pop();
            mBuf.unlock();

            if(USE_IMU)
            {
                if(!initFirstPoseFlag)//
                    initFirstIMUPose(accVector);//
                for(size_t i = 0; i < accVector.size(); i++)
                {
                    double dt;
                    if(i == 0)
                        dt = accVector[i].first - prevTime;
                    else if (i == accVector.size() - 1)
                        dt = curTime - accVector[i - 1].first;
                    else
                        dt = accVector[i].first - accVector[i - 1].first;
                    processIMU(accVector[i].first, dt, accVector[i].second, gyrVector[i].second);
                }
            }
            mProcess.lock();
            processImage(feature.second, feature.first);
            prevTime = curTime;
            //-------------------------------------------------------------------------
            for(size_t i = 0; i < accVector.size(); i++)
            {
                if(kf.update_time >0)
                {
                    if(predice_or_not<CONV_THRESHOLD)
                    {
                        kf.predict(ros::Time::now().toSec());//ublox
                        // kf.predict(accVector[i].first);
                    }
                    else
                    {
                        // kf.keep(accVector[i].first);
                        kf.keep(ros::Time::now().toSec());//ublox
                    }
                }
            }
            //-------------------------------------------------------------------------
            printStatistics(*this, 0);

            std_msgs::Header header;
            header.frame_id = "world";
            header.stamp = ros::Time(feature.first);
          
            pubTF(*this, header);//
            pubOdometry(*this, header);//
            if(get_P&&get_Q)
            {
                global_pubOdometry(*this, header);
                gps_pubOdometry(*this, header);
            }
            
            pubKeyPoses(*this, header);
            pubCameraPose(*this, header);//
            pubPointCloud(*this, header);//
            pubKeyframe(*this);
            //pubTF(*this, header);//
            mProcess.unlock();
        }

        if (! MULTIPLE_THREAD)
            break;

        std::chrono::milliseconds dura(2);
        std::this_thread::sleep_for(dura);
    }//while结束
}


void Estimator::initFirstIMUPose(vector<pair<double, Eigen::Vector3d>> &accVector)
{
    printf("init first imu pose\n");
    initFirstPoseFlag = true;
    //return;
    Eigen::Vector3d averAcc(0, 0, 0);
    int n = (int)accVector.size();
    for(size_t i = 0; i < accVector.size(); i++)
    {
        averAcc = averAcc + accVector[i].second;
    }
    averAcc = averAcc / n;
    printf("averge acc %f %f %f\n", averAcc.x(), averAcc.y(), averAcc.z());
    Matrix3d R0 = Utility::g2R(averAcc);
    double yaw = Utility::R2ypr(R0).x();
    R0 = Utility::ypr2R(Eigen::Vector3d{-yaw, 0, 0}) * R0;
    Rs[0] = R0;
    cout << "init R0 " << endl << Rs[0] << endl;
    //Vs[0] = Vector3d(5, 0, 0);
}

void Estimator::initFirstPose(Eigen::Vector3d p, Eigen::Matrix3d r)
{
    Ps[0] = p;
    Rs[0] = r;
    initP = p;
    initR = r;
}


void Estimator::processIMU(double t, double dt, const Vector3d &linear_acceleration, const Vector3d &angular_velocity)
{
    if (!first_imu)//初始为0
    {
        first_imu = true;
        acc_0 = linear_acceleration;
        gyr_0 = angular_velocity;
    }

    if (!pre_integrations[frame_count])
    {
        pre_integrations[frame_count] = new IntegrationBase{acc_0, gyr_0, Bas[frame_count], Bgs[frame_count]};
    }
    if (frame_count != 0)
    {
        pre_integrations[frame_count]->push_back(dt, linear_acceleration, angular_velocity);
        //if(solver_flag != NON_LINEAR)
            tmp_pre_integration->push_back(dt, linear_acceleration, angular_velocity);

        dt_buf[frame_count].push_back(dt);
        linear_acceleration_buf[frame_count].push_back(linear_acceleration);
        angular_velocity_buf[frame_count].push_back(angular_velocity);

        int j = frame_count;         
        Vector3d un_acc_0 = Rs[j] * (acc_0 - Bas[j]) - g;
        Vector3d un_gyr = 0.5 * (gyr_0 + angular_velocity) - Bgs[j];
        Rs[j] *= Utility::deltaQ(un_gyr * dt).toRotationMatrix();
        Vector3d un_acc_1 = Rs[j] * (linear_acceleration - Bas[j]) - g;
        Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
        Ps[j] += dt * Vs[j] + 0.5 * dt * dt * un_acc;
        Vs[j] += dt * un_acc;
    }
    acc_0 = linear_acceleration;
    gyr_0 = angular_velocity; 
}

void Estimator::processImage(const map<int, vector<pair<int, Eigen::Matrix<double, 7, 1>>>> &image, const double header)
{
    ROS_DEBUG("new image coming ------------------------------------------");
    ROS_DEBUG("Adding feature points %lu", image.size());
    if (f_manager.addFeatureCheckParallax(frame_count, image, td))
    {
        marginalization_flag = MARGIN_OLD;
        //printf("keyframe\n");
    }
    else
    {
        marginalization_flag = MARGIN_SECOND_NEW;
        //printf("non-keyframe\n");
    }

    ROS_DEBUG("%s", marginalization_flag ? "Non-keyframe" : "Keyframe");
    ROS_DEBUG("Solving %d", frame_count);
    ROS_DEBUG("number of feature: %d", f_manager.getFeatureCount());
    Headers[frame_count] = header;

    ImageFrame imageframe(image, header);
    imageframe.pre_integration = tmp_pre_integration;
    all_image_frame.insert(make_pair(header, imageframe));
    tmp_pre_integration = new IntegrationBase{acc_0, gyr_0, Bas[frame_count], Bgs[frame_count]};

    if(ESTIMATE_EXTRINSIC == 2)
    {
        ROS_INFO("calibrating extrinsic param, rotation movement is needed");
        if (frame_count != 0)
        {
            vector<pair<Vector3d, Vector3d>> corres = f_manager.getCorresponding(frame_count - 1, frame_count);
            Matrix3d calib_ric;
            if (initial_ex_rotation.CalibrationExRotation(corres, pre_integrations[frame_count]->delta_q, calib_ric))//
            {
                ROS_WARN("initial extrinsic rotation calib success");
                ROS_WARN_STREAM("initial extrinsic rotation: " << endl << calib_ric);
                ric[0] = calib_ric;
                RIC[0] = calib_ric;
                ESTIMATE_EXTRINSIC = 1;
            }
        }
    }

    if (solver_flag == INITIAL)
    {
        // monocular + IMU initilization
        if (!STEREO && USE_IMU)
        {
            if (frame_count == WINDOW_SIZE)
            {
                bool result = false;
                if(ESTIMATE_EXTRINSIC != 2 && (header - initial_timestamp) > 0.1)
                {
                    result = initialStructure();
                    initial_timestamp = header;   
                }
                if(result)
                {
                    optimization();
                    updateLatestStates();
                    solver_flag = NON_LINEAR;
                    slideWindow();
                    ROS_INFO("Initialization finish!");
                }
                else
                    slideWindow();
            }
        }

        // stereo + IMU initilization
        if(STEREO && USE_IMU)
        {
            f_manager.initFramePoseByPnP(frame_count, Ps, Rs, tic, ric);
            f_manager.triangulate(frame_count, Ps, Rs, tic, ric);
            if (frame_count == WINDOW_SIZE)
            {
                map<double, ImageFrame>::iterator frame_it;
                int i = 0;
                for (frame_it = all_image_frame.begin(); frame_it != all_image_frame.end(); frame_it++)
                {
                    frame_it->second.R = Rs[i];
                    frame_it->second.T = Ps[i];
                    i++;
                }
                solveGyroscopeBias(all_image_frame, Bgs);
                for (int i = 0; i <= WINDOW_SIZE; i++)
                {
                    pre_integrations[i]->repropagate(Vector3d::Zero(), Bgs[i]);
                }
                optimization();
                updateLatestStates();
                solver_flag = NON_LINEAR;
                slideWindow();
                ROS_INFO("Initialization finish!");
            }
        }

        // stereo only initilization
        if(STEREO && !USE_IMU)
        {
            f_manager.initFramePoseByPnP(frame_count, Ps, Rs, tic, ric);
            f_manager.triangulate(frame_count, Ps, Rs, tic, ric);
            optimization();

            if(frame_count == WINDOW_SIZE)
            {
                optimization();
                updateLatestStates();
                solver_flag = NON_LINEAR;
                slideWindow();
                ROS_INFO("Initialization finish!");
            }
        }

        if(frame_count < WINDOW_SIZE)
        {
            frame_count++;
            int prev_frame = frame_count - 1;
            Ps[frame_count] = Ps[prev_frame];
            Vs[frame_count] = Vs[prev_frame];
            Rs[frame_count] = Rs[prev_frame];
            Bas[frame_count] = Bas[prev_frame];
            Bgs[frame_count] = Bgs[prev_frame];
        }

    }
    else
    {
        TicToc t_solve;
        if(!USE_IMU)
            f_manager.initFramePoseByPnP(frame_count, Ps, Rs, tic, ric);
        f_manager.triangulate(frame_count, Ps, Rs, tic, ric);
        optimization();
        set<int> removeIndex;
        outliersRejection(removeIndex);
        f_manager.removeOutlier(removeIndex);
        if (! MULTIPLE_THREAD)
        {
            featureTracker.removeOutliers(removeIndex);
            predictPtsInNextFrame();
        }
            
        ROS_DEBUG("solver costs: %fms", t_solve.toc());

        if (failureDetection())
        {
            ROS_WARN("failure detection!");
            failure_occur = 1;
            clearState();
            setParameter();
            ROS_WARN("system reboot!");
            return;
        }

        slideWindow();
        f_manager.removeFailures();
        // prepare output of VINS
        key_poses.clear();
        for (int i = 0; i <= WINDOW_SIZE; i++)
            key_poses.push_back(Ps[i]);

        last_R = Rs[WINDOW_SIZE];
        last_P = Ps[WINDOW_SIZE];
        last_R0 = Rs[0];
        last_P0 = Ps[0];
        updateLatestStates();
    }  
}

bool Estimator::initialStructure()
{
    TicToc t_sfm;
    //check imu observibility
    {
        map<double, ImageFrame>::iterator frame_it;
        Vector3d sum_g;
        for (frame_it = all_image_frame.begin(), frame_it++; frame_it != all_image_frame.end(); frame_it++)
        {
            double dt = frame_it->second.pre_integration->sum_dt;
            Vector3d tmp_g = frame_it->second.pre_integration->delta_v / dt;
            sum_g += tmp_g;
        }
        Vector3d aver_g;
        aver_g = sum_g * 1.0 / ((int)all_image_frame.size() - 1);
        double var = 0;
        for (frame_it = all_image_frame.begin(), frame_it++; frame_it != all_image_frame.end(); frame_it++)
        {
            double dt = frame_it->second.pre_integration->sum_dt;
            Vector3d tmp_g = frame_it->second.pre_integration->delta_v / dt;
            var += (tmp_g - aver_g).transpose() * (tmp_g - aver_g);
            //cout << "frame g " << tmp_g.transpose() << endl;
        }
        var = sqrt(var / ((int)all_image_frame.size() - 1));
        //ROS_WARN("IMU variation %f!", var);
        if(var < 0.25)
        {
            ROS_INFO("IMU excitation not enouth!");
            //return false;
        }
    }
    // global sfm
    Quaterniond Q[frame_count + 1];
    Vector3d T[frame_count + 1];
    map<int, Vector3d> sfm_tracked_points;
    vector<SFMFeature> sfm_f;
    for (auto &it_per_id : f_manager.feature)
    {
        int imu_j = it_per_id.start_frame - 1;
        SFMFeature tmp_feature;
        tmp_feature.state = false;
        tmp_feature.id = it_per_id.feature_id;
        for (auto &it_per_frame : it_per_id.feature_per_frame)
        {
            imu_j++;
            Vector3d pts_j = it_per_frame.point;
            tmp_feature.observation.push_back(make_pair(imu_j, Eigen::Vector2d{pts_j.x(), pts_j.y()}));
        }
        sfm_f.push_back(tmp_feature);
    } 
    Matrix3d relative_R;
    Vector3d relative_T;
    int l;
    if (!relativePose(relative_R, relative_T, l))
    {
        ROS_INFO("Not enough features or parallax; Move device around");
        return false;
    }
    GlobalSFM sfm;
    if(!sfm.construct(frame_count + 1, Q, T, l,
              relative_R, relative_T,
              sfm_f, sfm_tracked_points))
    {
        ROS_DEBUG("global SFM failed!");
        marginalization_flag = MARGIN_OLD;
        return false;
    }

    //solve pnp for all frame
    map<double, ImageFrame>::iterator frame_it;
    map<int, Vector3d>::iterator it;
    frame_it = all_image_frame.begin( );
    for (int i = 0; frame_it != all_image_frame.end( ); frame_it++)
    {
        // provide initial guess
        cv::Mat r, rvec, t, D, tmp_r;
        if((frame_it->first) == Headers[i])
        {
            frame_it->second.is_key_frame = true;
            frame_it->second.R = Q[i].toRotationMatrix() * RIC[0].transpose();
            frame_it->second.T = T[i];
            i++;
            continue;
        }
        if((frame_it->first) > Headers[i])
        {
            i++;
        }
        Matrix3d R_inital = (Q[i].inverse()).toRotationMatrix();
        Vector3d P_inital = - R_inital * T[i];
        cv::eigen2cv(R_inital, tmp_r);
        cv::Rodrigues(tmp_r, rvec);
        cv::eigen2cv(P_inital, t);

        frame_it->second.is_key_frame = false;
        vector<cv::Point3f> pts_3_vector;
        vector<cv::Point2f> pts_2_vector;
        for (auto &id_pts : frame_it->second.points)
        {
            int feature_id = id_pts.first;
            for (auto &i_p : id_pts.second)
            {
                it = sfm_tracked_points.find(feature_id);
                if(it != sfm_tracked_points.end())
                {
                    Vector3d world_pts = it->second;
                    cv::Point3f pts_3(world_pts(0), world_pts(1), world_pts(2));
                    pts_3_vector.push_back(pts_3);
                    Vector2d img_pts = i_p.second.head<2>();
                    cv::Point2f pts_2(img_pts(0), img_pts(1));
                    pts_2_vector.push_back(pts_2);
                }
            }
        }
        cv::Mat K = (cv::Mat_<double>(3, 3) << 1, 0, 0, 0, 1, 0, 0, 0, 1);     
        if(pts_3_vector.size() < 6)
        {
            cout << "pts_3_vector size " << pts_3_vector.size() << endl;
            ROS_DEBUG("Not enough points for solve pnp !");
            return false;
        }
        if (! cv::solvePnP(pts_3_vector, pts_2_vector, K, D, rvec, t, 1))
        {
            ROS_DEBUG("solve pnp fail!");
            return false;
        }
        cv::Rodrigues(rvec, r);
        MatrixXd R_pnp,tmp_R_pnp;
        cv::cv2eigen(r, tmp_R_pnp);
        R_pnp = tmp_R_pnp.transpose();
        MatrixXd T_pnp;
        cv::cv2eigen(t, T_pnp);
        T_pnp = R_pnp * (-T_pnp);
        frame_it->second.R = R_pnp * RIC[0].transpose();
        frame_it->second.T = T_pnp;
    }
    if (visualInitialAlign())
        return true;
    else
    {
        ROS_INFO("misalign visual structure with IMU");
        return false;
    }

}

bool Estimator::visualInitialAlign()
{
    TicToc t_g;
    VectorXd x;
    //solve scale
    bool result = VisualIMUAlignment(all_image_frame, Bgs, g, x);
    if(!result)
    {
        ROS_DEBUG("solve g failed!");
        return false;
    }

    // change state
    for (int i = 0; i <= frame_count; i++)
    {
        Matrix3d Ri = all_image_frame[Headers[i]].R;
        Vector3d Pi = all_image_frame[Headers[i]].T;
        Ps[i] = Pi;
        Rs[i] = Ri;
        all_image_frame[Headers[i]].is_key_frame = true;
    }

    double s = (x.tail<1>())(0);
    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        pre_integrations[i]->repropagate(Vector3d::Zero(), Bgs[i]);
    }
    for (int i = frame_count; i >= 0; i--)
        Ps[i] = s * Ps[i] - Rs[i] * TIC[0] - (s * Ps[0] - Rs[0] * TIC[0]);
    int kv = -1;
    map<double, ImageFrame>::iterator frame_i;
    for (frame_i = all_image_frame.begin(); frame_i != all_image_frame.end(); frame_i++)
    {
        if(frame_i->second.is_key_frame)
        {
            kv++;
            Vs[kv] = frame_i->second.R * x.segment<3>(kv * 3);
        }
    }

    Matrix3d R0 = Utility::g2R(g);
    double yaw = Utility::R2ypr(R0 * Rs[0]).x();
    R0 = Utility::ypr2R(Eigen::Vector3d{-yaw, 0, 0}) * R0;
    g = R0 * g;
    //Matrix3d rot_diff = R0 * Rs[0].transpose();
    Matrix3d rot_diff = R0;
    for (int i = 0; i <= frame_count; i++)
    {
        Ps[i] = rot_diff * Ps[i];
        Rs[i] = rot_diff * Rs[i];
        Vs[i] = rot_diff * Vs[i];
    }
    ROS_DEBUG_STREAM("g0     " << g.transpose());
    ROS_DEBUG_STREAM("my R0  " << Utility::R2ypr(Rs[0]).transpose()); 

    f_manager.clearDepth();
    f_manager.triangulate(frame_count, Ps, Rs, tic, ric);

    return true;
}

bool Estimator::relativePose(Matrix3d &relative_R, Vector3d &relative_T, int &l)
{
    // find previous frame which contians enough correspondance and parallex with newest frame
    for (int i = 0; i < WINDOW_SIZE; i++)
    {
        vector<pair<Vector3d, Vector3d>> corres;
        corres = f_manager.getCorresponding(i, WINDOW_SIZE);
        if (corres.size() > 20)
        {
            double sum_parallax = 0;
            double average_parallax;
            for (int j = 0; j < int(corres.size()); j++)
            {
                Vector2d pts_0(corres[j].first(0), corres[j].first(1));
                Vector2d pts_1(corres[j].second(0), corres[j].second(1));
                double parallax = (pts_0 - pts_1).norm();
                sum_parallax = sum_parallax + parallax;

            }
            average_parallax = 1.0 * sum_parallax / int(corres.size());
            if(average_parallax * 460 > 30 && m_estimator.solveRelativeRT(corres, relative_R, relative_T))
            {
                l = i;
                ROS_DEBUG("average_parallax %f choose l %d and newest frame to triangulate the whole structure", average_parallax * 460, l);
                return true;
            }
        }
    }
    return false;
}

void Estimator::vector2double()
{
    for (int i = 0; i <= WINDOW_SIZE; i++)
    {
        para_Pose[i][0] = Ps[i].x();
        para_Pose[i][1] = Ps[i].y();
        para_Pose[i][2] = Ps[i].z();
        Quaterniond q{Rs[i]};
        para_Pose[i][3] = q.x();
        para_Pose[i][4] = q.y();
        para_Pose[i][5] = q.z();
        para_Pose[i][6] = q.w();

        if(USE_IMU)
        {
            para_SpeedBias[i][0] = Vs[i].x();
            para_SpeedBias[i][1] = Vs[i].y();
            para_SpeedBias[i][2] = Vs[i].z();

            para_SpeedBias[i][3] = Bas[i].x();
            para_SpeedBias[i][4] = Bas[i].y();
            para_SpeedBias[i][5] = Bas[i].z();

            para_SpeedBias[i][6] = Bgs[i].x();
            para_SpeedBias[i][7] = Bgs[i].y();
            para_SpeedBias[i][8] = Bgs[i].z();
        }
    }

    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        para_Ex_Pose[i][0] = tic[i].x();
        para_Ex_Pose[i][1] = tic[i].y();
        para_Ex_Pose[i][2] = tic[i].z();
        Quaterniond q{ric[i]};
        para_Ex_Pose[i][3] = q.x();
        para_Ex_Pose[i][4] = q.y();
        para_Ex_Pose[i][5] = q.z();
        para_Ex_Pose[i][6] = q.w();
    }


    VectorXd dep = f_manager.getDepthVector();//
    for (int i = 0; i < f_manager.getFeatureCount(); i++)
        para_Feature[i][0] = dep(i);

    para_Td[0][0] = td;
}

void Estimator::double2vector()
{
    Vector3d origin_R0 = Utility::R2ypr(Rs[0]);
    Vector3d origin_P0 = Ps[0];

    if (failure_occur)
    {
        origin_R0 = Utility::R2ypr(last_R0);
        origin_P0 = last_P0;
        failure_occur = 0;
    }

    if(USE_IMU)
    {
        Vector3d origin_R00 = Utility::R2ypr(Quaterniond(para_Pose[0][6],
                                                          para_Pose[0][3],
                                                          para_Pose[0][4],
                                                          para_Pose[0][5]).toRotationMatrix());
        double y_diff = origin_R0.x() - origin_R00.x();
        //TODO
        Matrix3d rot_diff = Utility::ypr2R(Vector3d(y_diff, 0, 0));
        if (abs(abs(origin_R0.y()) - 90) < 1.0 || abs(abs(origin_R00.y()) - 90) < 1.0)
        {
            ROS_DEBUG("euler singular point!");
            rot_diff = Rs[0] * Quaterniond(para_Pose[0][6],
                                           para_Pose[0][3],
                                           para_Pose[0][4],
                                           para_Pose[0][5]).toRotationMatrix().transpose();
        }

        for (int i = 0; i <= WINDOW_SIZE; i++)
        {

            Rs[i] = rot_diff * Quaterniond(para_Pose[i][6], para_Pose[i][3], para_Pose[i][4], para_Pose[i][5]).normalized().toRotationMatrix();
            
            Ps[i] = rot_diff * Vector3d(para_Pose[i][0] - para_Pose[0][0],
                                    para_Pose[i][1] - para_Pose[0][1],
                                    para_Pose[i][2] - para_Pose[0][2]) + origin_P0;


                Vs[i] = rot_diff * Vector3d(para_SpeedBias[i][0],
                                            para_SpeedBias[i][1],
                                            para_SpeedBias[i][2]);

                Bas[i] = Vector3d(para_SpeedBias[i][3],
                                  para_SpeedBias[i][4],
                                  para_SpeedBias[i][5]);

                Bgs[i] = Vector3d(para_SpeedBias[i][6],
                                  para_SpeedBias[i][7],
                                  para_SpeedBias[i][8]);
            
        }
    }
    else
    {
        for (int i = 0; i <= WINDOW_SIZE; i++)
        {
            Rs[i] = Quaterniond(para_Pose[i][6], para_Pose[i][3], para_Pose[i][4], para_Pose[i][5]).normalized().toRotationMatrix();
            
            Ps[i] = Vector3d(para_Pose[i][0], para_Pose[i][1], para_Pose[i][2]);
        }
    }

    if(USE_IMU)
    {
        for (int i = 0; i < NUM_OF_CAM; i++)
        {
            tic[i] = Vector3d(para_Ex_Pose[i][0],
                              para_Ex_Pose[i][1],
                              para_Ex_Pose[i][2]);
            ric[i] = Quaterniond(para_Ex_Pose[i][6],
                                 para_Ex_Pose[i][3],
                                 para_Ex_Pose[i][4],
                                 para_Ex_Pose[i][5]).normalized().toRotationMatrix();
        }
    }

    VectorXd dep = f_manager.getDepthVector();
    for (int i = 0; i < f_manager.getFeatureCount(); i++)
        dep(i) = para_Feature[i][0];
    f_manager.setDepth(dep);

    if(USE_IMU)
        td = para_Td[0][0];

}

bool Estimator::failureDetection()
{
    return false;
    if (f_manager.last_track_num < 2)
    {
        ROS_INFO(" little feature %d", f_manager.last_track_num);
        //return true;
    }
    if (Bas[WINDOW_SIZE].norm() > 2.5)
    {
        ROS_INFO(" big IMU acc bias estimation %f", Bas[WINDOW_SIZE].norm());
        return true;
    }
    if (Bgs[WINDOW_SIZE].norm() > 1.0)
    {
        ROS_INFO(" big IMU gyr bias estimation %f", Bgs[WINDOW_SIZE].norm());
        return true;
    }
    /*
    if (tic(0) > 1)
    {
        ROS_INFO(" big extri param estimation %d", tic(0) > 1);
        return true;
    }
    */
    Vector3d tmp_P = Ps[WINDOW_SIZE];
    if ((tmp_P - last_P).norm() > 5)
    {
        //ROS_INFO(" big translation");
        //return true;
    }
    if (abs(tmp_P.z() - last_P.z()) > 1)
    {
        //ROS_INFO(" big z translation");
        //return true; 
    }
    Matrix3d tmp_R = Rs[WINDOW_SIZE];
    Matrix3d delta_R = tmp_R.transpose() * last_R;
    Quaterniond delta_Q(delta_R);
    double delta_angle;
    delta_angle = acos(delta_Q.w()) * 2.0 / 3.14 * 180.0;
    if (delta_angle > 50)
    {
        ROS_INFO(" big delta_angle ");
        //return true;
    }
    return false;
}

void Estimator::optimization()
{
    // TicToc t_whole, t_prepare;
    vector2double();

    ceres::Problem problem;
    ceres::LossFunction *loss_function;
    //loss_function = NULL;
    loss_function = new ceres::HuberLoss(1.0);
    //loss_function = new ceres::CauchyLoss(1.0 / FOCAL_LENGTH);
    //ceres::LossFunction* loss_function = new ceres::HuberLoss(1.0);
    for (int i = 0; i < frame_count + 1; i++)
    {
        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();//???????
        problem.AddParameterBlock(para_Pose[i], SIZE_POSE, local_parameterization);
        if(USE_IMU)
            problem.AddParameterBlock(para_SpeedBias[i], SIZE_SPEEDBIAS);
    }
    if(!USE_IMU)
        problem.SetParameterBlockConstant(para_Pose[0]);

    for (int i = 0; i < NUM_OF_CAM; i++)
    {
        ceres::LocalParameterization *local_parameterization = new PoseLocalParameterization();
        problem.AddParameterBlock(para_Ex_Pose[i], SIZE_POSE, local_parameterization);
        if ((ESTIMATE_EXTRINSIC && frame_count == WINDOW_SIZE && Vs[0].norm() > 0.2) || openExEstimation)
        {
            //ROS_INFO("estimate extinsic param");
            openExEstimation = 1;
        }
        else
        {
            //ROS_INFO("fix extinsic param");
            problem.SetParameterBlockConstant(para_Ex_Pose[i]);//
        }
    }
    problem.AddParameterBlock(para_Td[0], 1);

    if (!ESTIMATE_TD || Vs[0].norm() < 0.2)
        problem.SetParameterBlockConstant(para_Td[0]);

    if (last_marginalization_info && last_marginalization_info->valid)
    {
        // construct new marginlization_factor
        MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);//
        problem.AddResidualBlock(marginalization_factor, NULL,
                                 last_marginalization_parameter_blocks);
    }
    if(USE_IMU)
    {
        for (int i = 0; i < frame_count; i++)
        {
            int j = i + 1;
            if (pre_integrations[j]->sum_dt > 10.0)
                continue;
            IMUFactor* imu_factor = new IMUFactor(pre_integrations[j]);
            problem.AddResidualBlock(imu_factor, NULL, para_Pose[i], para_SpeedBias[i], para_Pose[j], para_SpeedBias[j]);
        }
    }

    int f_m_cnt = 0;
    int feature_index = -1;
    for (auto &it_per_id : f_manager.feature)
    {
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (it_per_id.used_num < 4)
            continue;
 
        ++feature_index;

        int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;
        
        Vector3d pts_i = it_per_id.feature_per_frame[0].point;

        for (auto &it_per_frame : it_per_id.feature_per_frame)
        {
            imu_j++;
            if (imu_i != imu_j)
            {
                Vector3d pts_j = it_per_frame.point;
                ProjectionTwoFrameOneCamFactor *f_td = new ProjectionTwoFrameOneCamFactor(pts_i, pts_j, it_per_id.feature_per_frame[0].velocity, it_per_frame.velocity,
                                                                 it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td);
                problem.AddResidualBlock(f_td, loss_function, para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index], para_Td[0]);
            }

            if(STEREO && it_per_frame.is_stereo)
            {                
                Vector3d pts_j_right = it_per_frame.pointRight;
                if(imu_i != imu_j)
                {
                    ProjectionTwoFrameTwoCamFactor *f = new ProjectionTwoFrameTwoCamFactor(pts_i, pts_j_right, it_per_id.feature_per_frame[0].velocity, it_per_frame.velocityRight,
                                                                 it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td);
                    problem.AddResidualBlock(f, loss_function, para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Ex_Pose[1], para_Feature[feature_index], para_Td[0]);
                }
                else
                {
                    ProjectionOneFrameTwoCamFactor *f = new ProjectionOneFrameTwoCamFactor(pts_i, pts_j_right, it_per_id.feature_per_frame[0].velocity, it_per_frame.velocityRight,
                                                                 it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td);
                    problem.AddResidualBlock(f, loss_function, para_Ex_Pose[0], para_Ex_Pose[1], para_Feature[feature_index], para_Td[0]);
                }
               
            }
            f_m_cnt++;
        }
    }

    ROS_DEBUG("visual measurement count: %d", f_m_cnt);
    //printf("prepare for ceres: %f \n", t_prepare.toc());

    ceres::Solver::Options options;

    options.linear_solver_type = ceres::DENSE_SCHUR;
    //options.num_threads = 2;
    options.trust_region_strategy_type = ceres::DOGLEG;
    options.max_num_iterations = NUM_ITERATIONS;
    //options.use_explicit_schur_complement = true;
    //options.minimizer_progress_to_stdout = true;
    //options.use_nonmonotonic_steps = true;
    if (marginalization_flag == MARGIN_OLD)
        options.max_solver_time_in_seconds = SOLVER_TIME * 4.0 / 5.0;
    else
        options.max_solver_time_in_seconds = SOLVER_TIME;
    TicToc t_solver;
    ceres::Solver::Summary summary;
    ceres::Solve(options, &problem, &summary);
    //cout << summary.BriefReport() << endl;
    ROS_DEBUG("Iterations : %d", static_cast<int>(summary.iterations.size()));
    //printf("solver costs: %f \n", t_solver.toc());

    double2vector();
    //printf("frame_count: %d \n", frame_count);

    if(frame_count < WINDOW_SIZE)
        return;
    
    TicToc t_whole_marginalization;
    if (marginalization_flag == MARGIN_OLD)
    {
        MarginalizationInfo *marginalization_info = new MarginalizationInfo();
        vector2double();

        if (last_marginalization_info && last_marginalization_info->valid)
        {
            vector<int> drop_set;
            for (int i = 0; i < static_cast<int>(last_marginalization_parameter_blocks.size()); i++)
            {
                if (last_marginalization_parameter_blocks[i] == para_Pose[0] ||
                    last_marginalization_parameter_blocks[i] == para_SpeedBias[0])
                    drop_set.push_back(i);
            }
            // construct new marginlization_factor
            MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(marginalization_factor, NULL,
                                                                           last_marginalization_parameter_blocks,
                                                                           drop_set);
            marginalization_info->addResidualBlockInfo(residual_block_info);
        }

        if(USE_IMU)
        {
            if (pre_integrations[1]->sum_dt < 10.0)
            {
                IMUFactor* imu_factor = new IMUFactor(pre_integrations[1]);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(imu_factor, NULL,
                                                                           vector<double *>{para_Pose[0], para_SpeedBias[0], para_Pose[1], para_SpeedBias[1]},
                                                                           vector<int>{0, 1});
                marginalization_info->addResidualBlockInfo(residual_block_info);
            }
        }

        {
            int feature_index = -1;
            for (auto &it_per_id : f_manager.feature)
            {
                it_per_id.used_num = it_per_id.feature_per_frame.size();
                if (it_per_id.used_num < 4)
                    continue;

                ++feature_index;

                int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;
                if (imu_i != 0)
                    continue;

                Vector3d pts_i = it_per_id.feature_per_frame[0].point;

                for (auto &it_per_frame : it_per_id.feature_per_frame)
                {
                    imu_j++;
                    if(imu_i != imu_j)
                    {
                        Vector3d pts_j = it_per_frame.point;
                        ProjectionTwoFrameOneCamFactor *f_td = new ProjectionTwoFrameOneCamFactor(pts_i, pts_j, it_per_id.feature_per_frame[0].velocity, it_per_frame.velocity,
                                                                          it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td);
                        ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f_td, loss_function,
                                                                                        vector<double *>{para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Feature[feature_index], para_Td[0]},
                                                                                        vector<int>{0, 3});
                        marginalization_info->addResidualBlockInfo(residual_block_info);
                    }
                    if(STEREO && it_per_frame.is_stereo)
                    {
                        Vector3d pts_j_right = it_per_frame.pointRight;
                        if(imu_i != imu_j)
                        {
                            ProjectionTwoFrameTwoCamFactor *f = new ProjectionTwoFrameTwoCamFactor(pts_i, pts_j_right, it_per_id.feature_per_frame[0].velocity, it_per_frame.velocityRight,
                                                                          it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td);
                            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, loss_function,
                                                                                           vector<double *>{para_Pose[imu_i], para_Pose[imu_j], para_Ex_Pose[0], para_Ex_Pose[1], para_Feature[feature_index], para_Td[0]},
                                                                                           vector<int>{0, 4});
                            marginalization_info->addResidualBlockInfo(residual_block_info);
                        }
                        else
                        {
                            ProjectionOneFrameTwoCamFactor *f = new ProjectionOneFrameTwoCamFactor(pts_i, pts_j_right, it_per_id.feature_per_frame[0].velocity, it_per_frame.velocityRight,
                                                                          it_per_id.feature_per_frame[0].cur_td, it_per_frame.cur_td);
                            ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(f, loss_function,
                                                                                           vector<double *>{para_Ex_Pose[0], para_Ex_Pose[1], para_Feature[feature_index], para_Td[0]},
                                                                                           vector<int>{2});
                            marginalization_info->addResidualBlockInfo(residual_block_info);
                        }
                    }
                }
            }
        }

        TicToc t_pre_margin;
        marginalization_info->preMarginalize();
        ROS_DEBUG("pre marginalization %f ms", t_pre_margin.toc());
        
        TicToc t_margin;
        marginalization_info->marginalize();
        ROS_DEBUG("marginalization %f ms", t_margin.toc());

        std::unordered_map<long, double *> addr_shift;
        for (int i = 1; i <= WINDOW_SIZE; i++)
        {
            addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i - 1];//只是作了一系列地址映射
            if(USE_IMU)
                addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i - 1];
        }
        for (int i = 0; i < NUM_OF_CAM; i++)
            addr_shift[reinterpret_cast<long>(para_Ex_Pose[i])] = para_Ex_Pose[i];

        addr_shift[reinterpret_cast<long>(para_Td[0])] = para_Td[0];
//111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
        // 检查 para_Pose[1] 是否映射到了 para_Pose[0]
        // long key_to_check = reinterpret_cast<long>(para_Pose[1]);
        // auto it = addr_shift.find(key_to_check);
        // if (it != addr_shift.end()) {
        //     std::cout << "Mapping found for para_Pose[1]: " << it->first << " -> " << it->second << std::endl;
        //     // 进一步，可以验证 it->second 是否等于 para_Pose[0]
        //     if (it->second == para_Pose[0]) {
        //         std::cout << "Mapping para_Pose[1] -> para_Pose[0] is CORRECT." << std::endl;
        //     } else {
        //         std::cout << "ERROR: Expected para_Pose[0] (" << para_Pose[0] 
        //                 << "), but got " << it->second << std::endl;
        //     }
        // } else {
        //     std::cout << "ERROR: No mapping found for para_Pose[1]!" << std::endl;
        // }
// 111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111111
        vector<double *> parameter_blocks = marginalization_info->getParameterBlocks(addr_shift);

        if (last_marginalization_info)
            delete last_marginalization_info;
        last_marginalization_info = marginalization_info;
        last_marginalization_parameter_blocks = parameter_blocks;
        
    }
    else
    {
        if (last_marginalization_info &&
            std::count(std::begin(last_marginalization_parameter_blocks), std::end(last_marginalization_parameter_blocks), para_Pose[WINDOW_SIZE - 1]))
        {

            MarginalizationInfo *marginalization_info = new MarginalizationInfo();
            vector2double();
            if (last_marginalization_info && last_marginalization_info->valid)
            {
                vector<int> drop_set;
                for (int i = 0; i < static_cast<int>(last_marginalization_parameter_blocks.size()); i++)
                {
                    ROS_ASSERT(last_marginalization_parameter_blocks[i] != para_SpeedBias[WINDOW_SIZE - 1]);
                    if (last_marginalization_parameter_blocks[i] == para_Pose[WINDOW_SIZE - 1])
                        drop_set.push_back(i);
                }
                // construct new marginlization_factor
                MarginalizationFactor *marginalization_factor = new MarginalizationFactor(last_marginalization_info);
                ResidualBlockInfo *residual_block_info = new ResidualBlockInfo(marginalization_factor, NULL,
                                                                               last_marginalization_parameter_blocks,
                                                                               drop_set);

                marginalization_info->addResidualBlockInfo(residual_block_info);
            }

            TicToc t_pre_margin;
            ROS_DEBUG("begin marginalization");
            marginalization_info->preMarginalize();
            ROS_DEBUG("end pre marginalization, %f ms", t_pre_margin.toc());

            TicToc t_margin;
            ROS_DEBUG("begin marginalization");
            marginalization_info->marginalize();
            ROS_DEBUG("end marginalization, %f ms", t_margin.toc());
            
            std::unordered_map<long, double *> addr_shift;
            for (int i = 0; i <= WINDOW_SIZE; i++)
            {
                if (i == WINDOW_SIZE - 1)
                    continue;
                else if (i == WINDOW_SIZE)
                {
                    addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i - 1];
                    if(USE_IMU)
                        addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i - 1];
                }
                else
                {
                    addr_shift[reinterpret_cast<long>(para_Pose[i])] = para_Pose[i];
                    if(USE_IMU)
                        addr_shift[reinterpret_cast<long>(para_SpeedBias[i])] = para_SpeedBias[i];
                }
            }
            for (int i = 0; i < NUM_OF_CAM; i++)
                addr_shift[reinterpret_cast<long>(para_Ex_Pose[i])] = para_Ex_Pose[i];

            addr_shift[reinterpret_cast<long>(para_Td[0])] = para_Td[0];

            
            vector<double *> parameter_blocks = marginalization_info->getParameterBlocks(addr_shift);
            if (last_marginalization_info)
                delete last_marginalization_info;
            last_marginalization_info = marginalization_info;
            last_marginalization_parameter_blocks = parameter_blocks;
            
        }
    }
    // std::ofstream foutC("/home/zk/小论文数据支持/cheku/all_opt_tlist.csv", ios::app);
    // foutC.setf(ios::fixed, ios::floatfield);
    // foutC.precision(6);
    // foutC << t_whole.toc() << endl;
    // foutC.close();
    //printf("whole marginalization costs: %f \n", t_whole_marginalization.toc());
    // printf("whole time for ceres: %f \n", t_whole.toc());
}

void Estimator::slideWindow()
{
    TicToc t_margin;
    if (marginalization_flag == MARGIN_OLD)
    {
        double t_0 = Headers[0];
        back_R0 = Rs[0];
        back_P0 = Ps[0];
        if (frame_count == WINDOW_SIZE)
        {
            for (int i = 0; i < WINDOW_SIZE; i++)
            {
                Headers[i] = Headers[i + 1];
                Rs[i].swap(Rs[i + 1]);
                Ps[i].swap(Ps[i + 1]);
                if(USE_IMU)
                {
                    std::swap(pre_integrations[i], pre_integrations[i + 1]);

                    dt_buf[i].swap(dt_buf[i + 1]);
                    linear_acceleration_buf[i].swap(linear_acceleration_buf[i + 1]);
                    angular_velocity_buf[i].swap(angular_velocity_buf[i + 1]);

                    Vs[i].swap(Vs[i + 1]);
                    Bas[i].swap(Bas[i + 1]);
                    Bgs[i].swap(Bgs[i + 1]);
                }
            }
            Headers[WINDOW_SIZE] = Headers[WINDOW_SIZE - 1];
            Ps[WINDOW_SIZE] = Ps[WINDOW_SIZE - 1];
            Rs[WINDOW_SIZE] = Rs[WINDOW_SIZE - 1];

            if(USE_IMU)
            {
                Vs[WINDOW_SIZE] = Vs[WINDOW_SIZE - 1];
                Bas[WINDOW_SIZE] = Bas[WINDOW_SIZE - 1];
                Bgs[WINDOW_SIZE] = Bgs[WINDOW_SIZE - 1];

                delete pre_integrations[WINDOW_SIZE];
                pre_integrations[WINDOW_SIZE] = new IntegrationBase{acc_0, gyr_0, Bas[WINDOW_SIZE], Bgs[WINDOW_SIZE]};

                dt_buf[WINDOW_SIZE].clear();
                linear_acceleration_buf[WINDOW_SIZE].clear();
                angular_velocity_buf[WINDOW_SIZE].clear();
            }

            if (true || solver_flag == INITIAL)
            {
                map<double, ImageFrame>::iterator it_0;
                it_0 = all_image_frame.find(t_0);
                delete it_0->second.pre_integration;
                all_image_frame.erase(all_image_frame.begin(), it_0);
            }
            slideWindowOld();
        }
    }
    else
    {
        if (frame_count == WINDOW_SIZE)
        {
            Headers[frame_count - 1] = Headers[frame_count];
            Ps[frame_count - 1] = Ps[frame_count];
            Rs[frame_count - 1] = Rs[frame_count];

            if(USE_IMU)
            {
                for (unsigned int i = 0; i < dt_buf[frame_count].size(); i++)
                {
                    double tmp_dt = dt_buf[frame_count][i];
                    Vector3d tmp_linear_acceleration = linear_acceleration_buf[frame_count][i];
                    Vector3d tmp_angular_velocity = angular_velocity_buf[frame_count][i];

                    pre_integrations[frame_count - 1]->push_back(tmp_dt, tmp_linear_acceleration, tmp_angular_velocity);

                    dt_buf[frame_count - 1].push_back(tmp_dt);
                    linear_acceleration_buf[frame_count - 1].push_back(tmp_linear_acceleration);
                    angular_velocity_buf[frame_count - 1].push_back(tmp_angular_velocity);
                }

                Vs[frame_count - 1] = Vs[frame_count];
                Bas[frame_count - 1] = Bas[frame_count];
                Bgs[frame_count - 1] = Bgs[frame_count];

                delete pre_integrations[WINDOW_SIZE];
                pre_integrations[WINDOW_SIZE] = new IntegrationBase{acc_0, gyr_0, Bas[WINDOW_SIZE], Bgs[WINDOW_SIZE]};

                dt_buf[WINDOW_SIZE].clear();
                linear_acceleration_buf[WINDOW_SIZE].clear();
                angular_velocity_buf[WINDOW_SIZE].clear();
            }
            slideWindowNew();
        }
    }
}

void Estimator::slideWindowNew()
{
    sum_of_front++;
    f_manager.removeFront(frame_count);
}

void Estimator::slideWindowOld()
{
    sum_of_back++;

    bool shift_depth = solver_flag == NON_LINEAR ? true : false;
    if (shift_depth)
    {
        Matrix3d R0, R1;
        Vector3d P0, P1;
        R0 = back_R0 * ric[0];
        R1 = Rs[0] * ric[0];
        P0 = back_P0 + back_R0 * tic[0];
        P1 = Ps[0] + Rs[0] * tic[0];
        f_manager.removeBackShiftDepth(R0, P0, R1, P1);
    }
    else
        f_manager.removeBack();
}


void Estimator::getPoseInWorldFrame(Eigen::Matrix4d &T)
{
    T = Eigen::Matrix4d::Identity();
    T.block<3, 3>(0, 0) = Rs[frame_count];
    T.block<3, 1>(0, 3) = Ps[frame_count];
}

void Estimator::getPoseInWorldFrame(int index, Eigen::Matrix4d &T)
{
    T = Eigen::Matrix4d::Identity();
    T.block<3, 3>(0, 0) = Rs[index];
    T.block<3, 1>(0, 3) = Ps[index];
}

void Estimator::predictPtsInNextFrame()
{
    //printf("predict pts in next frame\n");
    if(frame_count < 2)
        return;
    // predict next pose. Assume constant velocity motion
    Eigen::Matrix4d curT, prevT, nextT;
    getPoseInWorldFrame(curT);
    getPoseInWorldFrame(frame_count - 1, prevT);
    nextT = curT * (prevT.inverse() * curT);
    map<int, Eigen::Vector3d> predictPts;

    for (auto &it_per_id : f_manager.feature)
    {
        if(it_per_id.estimated_depth > 0)
        {
            int firstIndex = it_per_id.start_frame;
            int lastIndex = it_per_id.start_frame + it_per_id.feature_per_frame.size() - 1;
            //printf("cur frame index  %d last frame index %d\n", frame_count, lastIndex);
            if((int)it_per_id.feature_per_frame.size() >= 2 && lastIndex == frame_count)
            {
                double depth = it_per_id.estimated_depth;
                Vector3d pts_j = ric[0] * (depth * it_per_id.feature_per_frame[0].point) + tic[0];
                Vector3d pts_w = Rs[firstIndex] * pts_j + Ps[firstIndex];
                Vector3d pts_local = nextT.block<3, 3>(0, 0).transpose() * (pts_w - nextT.block<3, 1>(0, 3));
                Vector3d pts_cam = ric[0].transpose() * (pts_local - tic[0]);
                int ptsIndex = it_per_id.feature_id;
                predictPts[ptsIndex] = pts_cam;
            }
        }
    }
    featureTracker.setPrediction(predictPts);
    //printf("estimator output %d predict pts\n",(int)predictPts.size());
}

double Estimator::reprojectionError(Matrix3d &Ri, Vector3d &Pi, Matrix3d &rici, Vector3d &tici,
                                 Matrix3d &Rj, Vector3d &Pj, Matrix3d &ricj, Vector3d &ticj, 
                                 double depth, Vector3d &uvi, Vector3d &uvj)
{
    Vector3d pts_w = Ri * (rici * (depth * uvi) + tici) + Pi;
    Vector3d pts_cj = ricj.transpose() * (Rj.transpose() * (pts_w - Pj) - ticj);
    Vector2d residual = (pts_cj / pts_cj.z()).head<2>() - uvj.head<2>();
    double rx = residual.x();
    double ry = residual.y();
    return sqrt(rx * rx + ry * ry);
}

void Estimator::outliersRejection(set<int> &removeIndex)
{
    //return;
    int feature_index = -1;
    for (auto &it_per_id : f_manager.feature)
    {
        double err = 0;
        int errCnt = 0;
        it_per_id.used_num = it_per_id.feature_per_frame.size();
        if (it_per_id.used_num < 4)
            continue;
        feature_index ++;
        int imu_i = it_per_id.start_frame, imu_j = imu_i - 1;
        Vector3d pts_i = it_per_id.feature_per_frame[0].point;
        double depth = it_per_id.estimated_depth;
        for (auto &it_per_frame : it_per_id.feature_per_frame)
        {
            imu_j++;
            if (imu_i != imu_j)
            {
                Vector3d pts_j = it_per_frame.point;             
                double tmp_error = reprojectionError(Rs[imu_i], Ps[imu_i], ric[0], tic[0], 
                                                    Rs[imu_j], Ps[imu_j], ric[0], tic[0],
                                                    depth, pts_i, pts_j);
                err += tmp_error;
                errCnt++;
                //printf("tmp_error %f\n", FOCAL_LENGTH / 1.5 * tmp_error);
            }
            // need to rewrite projecton factor.........
            if(STEREO && it_per_frame.is_stereo)
            {
                
                Vector3d pts_j_right = it_per_frame.pointRight;
                if(imu_i != imu_j)
                {            
                    double tmp_error = reprojectionError(Rs[imu_i], Ps[imu_i], ric[0], tic[0], 
                                                        Rs[imu_j], Ps[imu_j], ric[1], tic[1],
                                                        depth, pts_i, pts_j_right);
                    err += tmp_error;
                    errCnt++;
                    //printf("tmp_error %f\n", FOCAL_LENGTH / 1.5 * tmp_error);
                }
                else
                {
                    double tmp_error = reprojectionError(Rs[imu_i], Ps[imu_i], ric[0], tic[0], 
                                                        Rs[imu_j], Ps[imu_j], ric[1], tic[1],
                                                        depth, pts_i, pts_j_right);
                    err += tmp_error;
                    errCnt++;
                    //printf("tmp_error %f\n", FOCAL_LENGTH / 1.5 * tmp_error);
                }       
            }
        }
        double ave_err = err / errCnt;
        if(ave_err * FOCAL_LENGTH > 3)
            removeIndex.insert(it_per_id.feature_id);

    }
}


void Estimator::fastPredictIMU(double t, Eigen::Vector3d linear_acceleration, Eigen::Vector3d angular_velocity)
{
    double dt = t - latest_time;
    latest_time = t;
    Eigen::Vector3d un_acc_0 = latest_Q * (latest_acc_0 - latest_Ba) - g;
    Eigen::Vector3d un_gyr = 0.5 * (latest_gyr_0 + angular_velocity) - latest_Bg;
    latest_Q = latest_Q * Utility::deltaQ(un_gyr * dt);
    latest_R = latest_Q.normalized().toRotationMatrix();
    Eigen::Vector3d un_acc_1 = latest_Q * (linear_acceleration - latest_Ba) - g;
    Eigen::Vector3d un_acc = 0.5 * (un_acc_0 + un_acc_1);
    latest_P = latest_P + dt * latest_V + 0.5 * dt * dt * un_acc;
    latest_V = latest_V + dt * un_acc;
    latest_acc_0 = linear_acceleration;
    latest_gyr_0 = angular_velocity;

    if(optimized)
    {
        // cout<<"step_one_P_after_esr_R111111111111111111111111111111111111111111111111111111imu_tehua"<<l_to_g_p<<endl;
        l_to_g_p = const_est_R*latest_P + const_est_t;
        // double t_zitai = ros::Time::now().toSec();
        // temp_t_l_to_g_p = make_pair(t_zitai,l_to_g_p);
        V_n = const_est_R * latest_V;
        // cout<<"latest_v = "<<V_n<<endl;
        Q_n_kf = const_est_R*latest_Q*RIC[0]*R_s_c_enu.transpose();//检查
    }
    if(get_P&&get_Q)
    {
        Eigen::Vector3d global_un_acc_0 = final_Q * (latest_acc_0 - latest_Ba) - const_est_R*g;
        // Eigen::Vector3d global_un_gyr = 0.5 * (latest_gyr_0 + angular_velocity) - latest_Bg;//使用同一个
        final_Q = final_Q * Utility::deltaQ(un_gyr * dt);
        Eigen::Vector3d global_un_acc_1 = final_Q * (linear_acceleration - latest_Ba) - const_est_R*g;
        Eigen::Vector3d global_un_acc = 0.5 * (global_un_acc_0 + global_un_acc_1);
        final_P = final_P + dt * final_V + 0.5 * dt * dt * global_un_acc;
        // cout<<"predict_imu_dt = "<<dt<<endl;
        // cout<<"final_ppppppp = "<<final_P<<endl;
        final_V = final_V + dt * global_un_acc;
        // cout<<"final_vvvvvvv = "<<final_V<<endl;
    }

}

void Estimator::updateLatestStates()
{
    R_s_c_enu << 0, 0, 1,
                -1, 0, 0,
                0, -1, 0;
    mPropagate.lock();
    latest_time = Headers[frame_count] + td;
    double t_local = ros::Time::now().toSec();
    // cout<<"update_latested_td================================== "<<td<<endl;
    latest_P = Ps[frame_count];
    latest_Q = Rs[frame_count];
    latest_V = Vs[frame_count];
    latest_Ba = Bas[frame_count];
    latest_Bg = Bgs[frame_count];
    latest_acc_0 = acc_0;
    latest_gyr_0 = gyr_0;

    if(!optimized)
    {

        // stimulate = get_stimulate(accBuf);
        if((!accBuf.empty())&&(!gyrBuf.empty()))
        {
            stimulate = get_stimulate(accBuf , gyrBuf);
            if(stimulate>0.6)
            {
                stimulating  = true;
            }
            else stimulating = false;
            if(stimulating)
            {
                local_odom_buf[t_local] = latest_P;
            }
        }
        
    }
    else if(optimized)
    {
        // double t_zitai = ros::Time::now().toSec();
        // temp_t_l_to_g_p = make_pair(t_zitai,l_to_g_p);
        
        // time_of_vio_used = ros::Time::now().toSec();

        if(latest_V.x()*latest_V.x()+latest_V.y()*latest_V.y() > 3)
        {
            local_odom_buf[t_local] = latest_P;
        }
        if (local_odom_buf.size() > 400) 
        {
            std::map<double, Eigen::Vector3d>::iterator it = local_odom_buf.begin();
            while (local_odom_buf.size() > 400) 
            {
                local_odom_buf.erase(it); // erase会将迭代器失效，所以不需要手动增加迭代器
            }
        }

        l_to_g_p = const_est_R*latest_P + const_est_t;
        V_n = const_est_R * latest_V;
        // Q_n_kf = const_est_R*latest_Q*RIC[0]*R_s_c_enu.transpose();//检查
        Q_n_kf = const_est_R*latest_Q;
        // cout<<"step_one_P_after_esr_R111111111111111111111111111111111111111111111111111111"<<l_to_g_p<<endl;
        Eigen::Vector3d yk = l_to_g_p - l_to_g_p;

        if(kf.update_time == -1)
        {
            double t_init = ros::Time::now().toSec();
            kf.initState(yk,const_est_R,t_init);//
            // kf.initState(yk,const_est_R,latest_time);//
        }
            
        // else
        // {
        //     if(predice_or_not<=1.0 && ready_to_predict)
        //     {
        //         kf.predict(t_local);
        //         // cout<<"kf.x(1)=====================================loacl================================"<<kf.x(1)<<endl;
        //     }
        //     else if(predice_or_not>1.0||(!ready_to_predict))//一开始先keep
        //     {
        //         kf.keep(t_local);
        //     }
        //     //  cout<<"state_in_localaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa = "<<kf.x<<endl;
        // }
    }
    
    

    mBuf.lock();
    queue<pair<double, Eigen::Vector3d>> tmp_accBuf = accBuf;
    queue<pair<double, Eigen::Vector3d>> tmp_gyrBuf = gyrBuf;
    mBuf.unlock();
    while(!tmp_accBuf.empty())
    {
        double t = tmp_accBuf.front().first;
        Eigen::Vector3d acc = tmp_accBuf.front().second;
        Eigen::Vector3d gyr = tmp_gyrBuf.front().second;
        fastPredictIMU(t, acc, gyr);
        tmp_accBuf.pop();
        tmp_gyrBuf.pop();
    }
    mPropagate.unlock();
}
