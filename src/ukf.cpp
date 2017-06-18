#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <math.h>

#define CLOSE_ENOUGH 0.001

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // initializatoin state
  is_initialized_ = false;
  
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // State dimension 
  n_x_ = 5;

  // Augmented state dimension
  n_aug_ = 7;

  // Sigma point measurement parameter
  lambda_ = 3.0 - n_aug_;

  // Current NIS (radar)
  NIS_radar_ = 0;

  // Current NIS (lidar)
  NIS_laser_ = 0; 

  weights_ = VectorXd(2 * n_aug_ + 1);
  
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);

  // initial time
  time_us_ = 0.0;
  
  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
  
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.5;

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;
  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  Q_ = MatrixXd(2,2);  

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  Tools tools;
  
  if (not is_initialized_) {

    x_.fill(0.0);
    P_ << 1.0, 0, 0, 0, 0,
          0, 1.0, 0, 0, 0,
          0, 0, 1.0, 0, 0,
          0, 0, 0, 1.0, 0,
          0, 0, 0, 0, 1.0;

    Q_ << pow(std_a_, 2), 0,
          0, pow(std_yawdd_, 2);

    
    if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      float px = meas_package.raw_measurements_[0];
      float py = meas_package.raw_measurements_[1];
      
      x_ << px, py, 0, 0, 0;
    }
  
      if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {

	// radar gives us (rho, phi, phi_dot) 
	float rho = meas_package.raw_measurements_[0];
	float phi = meas_package.raw_measurements_[1];

	// state space (px, py, v, psi, psi_dot) 
	float px = rho * cos(phi);
	float py = rho * sin(phi);

	x_ << px, py, 0, 0, 0; 

      }
      
      double weight_0 = lambda_ / (lambda_ + n_aug_);
      weights_(0) = weight_0;

      for (int i=1; i<2*n_aug_+1; i++) {  
	double weight = 0.5/(n_aug_ + lambda_);
	weights_(i) = weight;
      }
      is_initialized_ = true; 
  }


  float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;

  // restrict nonsense (for the first timestamp in particular)
  if (delta_t > 100 * CLOSE_ENOUGH) {
    delta_t = 0.05;
  }
  
  Prediction(delta_t);
  
  if (meas_package.sensor_type_ == MeasurementPackage::LASER && use_laser_) {
    UpdateLidar(meas_package);
  }
  
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_) {
    UpdateRadar(meas_package);
  }

  cout << "x_: " << x_ << endl;
  cout << "P_: " << P_ << endl;
  cout << "NIS_radar_: " << NIS_radar_ << endl;
  cout << "NIS_laser_: " << NIS_laser_ << endl; 
  cout << "_________________________________" << endl;
  
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  Tools tools; 
  
  // augmented state 
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  // augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_ , n_x_) = P_;
  P_aug.bottomRightCorner(2, 2)    = Q_; 
  
  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  
  //augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug.col(0)   = x_aug;  
  for (int i = 0; i< n_aug_; i++) {
    Xsig_aug.col(i+1)        = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  // predict sigma points
  Xsig_pred_.fill(0.0);
  for (int i = 0; i< 2*n_aug_+1; i++) {
    // extract values for better readability
    float p_x      = Xsig_aug(0,i);
    float p_y      = Xsig_aug(1,i);
    float v        = Xsig_aug(2,i);
    float yaw      = Xsig_aug(3,i);
    float yawd     = Xsig_aug(4,i);
    float nu_a     = Xsig_aug(5,i);
    float nu_yawdd = Xsig_aug(6,i);

    // predicted state values
    float px_p, py_p;
    
    // avoid division by zero
    if (fabs(yawd) >= CLOSE_ENOUGH) {
      px_p = p_x + v / yawd * (sin(yaw  + yawd * delta_t) - sin(yaw));
      py_p = p_y + v / yawd * (cos(yaw) - cos(yaw + yawd * delta_t));
    } else {
      px_p = p_x + v * delta_t * cos(yaw);	
      py_p = p_y + v * delta_t * sin(yaw);
    }
    
    float v_p    = v;
    float yaw_p  = yaw + yawd * delta_t;
    float yawd_p = yawd;

    // add noise
    px_p = px_p + 0.5  * nu_a * delta_t * delta_t * cos(yaw);
    py_p = py_p + 0.5  * nu_a * delta_t * delta_t * sin(yaw);

    v_p  = v_p  + nu_a * delta_t;

    yaw_p  = yaw_p  + 0.5 * nu_yawdd * delta_t * delta_t;
    yawd_p = yawd_p + nu_yawdd * delta_t;

    // write predicted sigma point into right column
    Xsig_pred_.col(i) << px_p, py_p, v_p, yaw_p, yawd_p;
  }

  x_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    x_ = x_ + weights_(i) * Xsig_pred_.col(i);
  }
  
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    //angle normalization
    x_diff(3) = tools.WrapAnglePi(x_diff(3));
    
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
    
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  Tools tools;
  
  int n_z  = 2;

  // transform sigma points into measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);
  Zsig.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    float p_x = Xsig_pred_(0,i);
    float p_y = Xsig_pred_(1,i);
    Zsig.col(i) << p_x, p_y;
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R << pow(std_laspx_, 2), 0,
       0, pow(std_laspy_, 2);

  Update(n_z, Zsig, R, meas_package);

}


/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  Tools tools;
  
  // radar give us a measurement vector of length 3 
  int n_z = 3;
  MatrixXd Zsig = MatrixXd(n_z, 2*n_aug_+1);
  
  // transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v   = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    float r     = sqrt(p_x*p_x + p_y*p_y);                        
    float phi   = atan2(p_y,p_x);                                 
    float r_dot = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);

    Zsig.col(i) << r, phi, r_dot; 
    
  }

  // add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z,n_z);
  R << pow(std_radr_, 2), 0, 0,
       0, pow(std_radphi_, 2), 0,
       0, 0, pow(std_radrd_, 2);

  Update(n_z, Zsig, R, meas_package);
  
}

void UKF::Update(int n_z, MatrixXd Zsig, MatrixXd R, MeasurementPackage meas_package){

  Tools tools; 
  
  // mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i=0; i < 2*n_aug_+1; i++) {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  // measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {
    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = tools.WrapAnglePi(z_diff(1)); 
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  S = S + R;

   // calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++) {

    // residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    // angle normalization
    z_diff(1) = tools.WrapAnglePi(z_diff(1)); 

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    // angle normalization
    x_diff(3) = tools.WrapAnglePi(x_diff(3));

    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  // residual
  VectorXd z_diff = meas_package.raw_measurements_ - z_pred;
  z_diff(1) = tools.WrapAnglePi(z_diff(1));

  // NIS
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    NIS_radar_ = z_diff.transpose() * S.inverse() * z_diff;
  }
  if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
    NIS_laser_ = z_diff.transpose() * S.inverse() * z_diff;
  }

  // update state mean and covariance matrix
  MatrixXd K = Tc * S.inverse();
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  
}
