#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;
  is_initialized_ = false;
  
  /////////////Parameter dimensions
  n_x_ = 5;
  n_aug_ = n_x_ + 2;
  n_z_laser = 2;
  n_z_radar = 3;
   
  
  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
 
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 30;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 30;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
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
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.
  
  ///Process noise matrix for prediction
  Q = MatrixXd(2,2);
  Q << std_a_ * std_a_, 0,
        0, std_yawdd_ * std_yawdd_;
  
  ///Measurement noise matrix for laser
  R_laser = MatrixXd(n_z_laser,n_z_laser);
  R_laser << std_laspx_*std_laspx_, 0,
            0, std_laspy_*std_laspy_;
  
  ///Measurement noise matrix for radar
  R_radar = MatrixXd(n_z_radar,n_z_radar);
  R_radar << std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0, std_radrd_*std_radrd_;
  
  
  // set weights
  weights = VectorXd(2*n_aug_+1);
  double weight = (1.0*lambda_)/(lambda_+n_aug_);
  weights(0) = weight;
  weight = 0.5/(n_aug_+lambda_);
  for (int i=1; i<2*n_aug_+1; i++) {  //2n+1 weights
    weights(i) = weight;
  }
  
  
  P_ = MatrixXd(4,4);
/* 	P_ <<  230.7548681 , 0, 0, 0,
			0, 82.98746015, 0, 0,
			0, 0, 13.93748925, 0,
			0, 0, 0, 10.82690745; */
	
	P_ <<  230.5 , 0, 0, 0,
			0, 82.9, 0, 0,
			0, 0, 13.9, 0,
			0, 0, 0, 10.8;
			
/* 	P_ <<  1 , 0, 0, 0,
			0, 1, 0, 0,
			0, 0, 1, 0,
			0, 0, 0, 1; */
  
  
  /**
  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    x_= VectorXd(n_x_);
    previous_timestamp_ = meas_package.timestamp_;
        
        
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      //float rho = meas_package.raw_measurements_[0];
      //float phi = meas_package.raw_measurements_[1];
      
      //x_ << rho*cos(phi), rho*sin(phi), 0, 0;
      // velocity = sqrt(5.199747^2 + 0.001796856^2);  //
      x_ << 0.8599968, 0.6000449, 5.19974731, 0.0003455661, 0.01382155;
      //cout << "initial x_ in radar= " << x_ << endl;
      is_initialized_ = true;
      return;	  
	  
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      
      /* Initialize state. */
		//set the state with the initial location and zero velocity
		
      //x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0;
      x_ << 0.6, 0.6, 5.199937, 0, 0.006911322;
      //cout << "initial x_ in laser= " << x_ << endl;
      is_initialized_ = true;
      return;
    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  
  }
  
  double delta_t = (meas_package.timestamp_ - previous_timestamp_)/1000000.0;	//delta_t - expressed in seconds
	previous_timestamp_ = meas_package.timestamp_;
  //////Calling UKF predict function for each iteration.
  Prediction(delta_t);
  
  /////Calling measurement
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR and use_radar_ == true) {
    UpdateRadar(meas_package);
    
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER  and use_laser_==true) {  
    UpdateLidar(meas_package);
  
  }

  return;
  
}



/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
    //create augmented mean state
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug << x_, 0, 0;

  //create augmented covariance matrix
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug.bottomRightCorner(2,2) = Q;
  
    //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  
  //create augmented sigma points
  lambda_ = 3-n_x_;
  float lamda_n_x_root = sqrt(lambda_+n_x_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  Xsig_aug.col(0) << x_aug;
  for (unsigned int i=0; i < n_aug_ ; i++){
    Xsig_aug.col(i+1) << x_aug + lamda_n_x_root * A.col(i);
    Xsig_aug.col(i+1+n_aug_) << x_aug - lamda_n_x_root * A.col(i);
      
  }
  
  
  // apply prediction function on sigma points 
  for (unsigned int i = 0; i< 2*n_aug_+1; i++)
  {
    //extract values for better readability
    double p_x = Xsig_aug(0,i);
    double p_y = Xsig_aug(1,i);
    double v = Xsig_aug(2,i);
    double yaw = Xsig_aug(3,i);
    double yawd = Xsig_aug(4,i);
    double nu_a = Xsig_aug(5,i);
    double nu_yawdd = Xsig_aug(6,i);

    //predicted state values
    double px_p, py_p;

    //avoid division by zero
    if (fabs(yawd) > 0.001) {
        px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
        py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
    }
    else {
        px_p = p_x + v*delta_t*cos(yaw);
        py_p = p_y + v*delta_t*sin(yaw);
    }

    double v_p = v;
    double yaw_p = yaw + yawd*delta_t;
    double yawd_p = yawd;

    //add noise
    px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
    py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
    v_p = v_p + nu_a*delta_t;

    yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
    yawd_p = yawd_p + nu_yawdd*delta_t;

    //write predicted sigma point into right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
  
  ///////////////predict covariance function
  
  //predicted state mean
  VectorXd x_mean = VectorXd(n_aug_);
  x_mean.fill(0.0);
  for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_mean = x_mean + weights(i) * Xsig_pred_.col(i);
  }

  //predicted state covariance matrix
  MatrixXd P = MatrixXd(n_x_,n_x_);
  P.fill(0.0);
  for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_mean;
    //angle normalization
/*     while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI; */
    
    x_diff(3) = atan2(cos(x_diff(3)), sin(x_diff(3)) );
    

    P = P + weights(i) * x_diff * x_diff.transpose() ;
  }
  
  x_ = x_mean;
  P_ = P;
  //return from prediction step
   return;
  
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  ///updating number of measurement points i.e 2 for lider/laser
  n_z_ = n_z_laser;
  R = R_laser;
 
//transform sigma points into measurement space would be direct mapping of x,y sigma points 
  Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  Zsig.row(0) = Xsig_pred_.row(0); /// x measurement
  Zsig.row(1) = Xsig_pred_.row(1); /// y measurement   
    
  Update(meas_package);

  return;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  ///updating number of measurement points i.e 3 for radar
  n_z_ = n_z_radar;
  
  R = R_radar;
  //////////////////////////////////Predict Radar Measurement
    //create matrix for sigma points in measurement space
  Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  
    //transform sigma points into measurement space
  for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // extract values for better readibility
    double p_x = Xsig_pred_(0,i);
    double p_y = Xsig_pred_(1,i);
    double v  = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    double v1 = cos(yaw)*v;
    double v2 = sin(yaw)*v;

    // measurement model
    Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                        //r
    Zsig(1,i) = atan2(p_y,p_x);                                 //phi
    Zsig(2,i) = (p_x*v1 + p_y*v2 ) / sqrt(p_x*p_x + p_y*p_y);   //r_dot
  }
  
  Update(meas_package);

  return;
}


//////comman update function to be called for both sensor update functions: UpdateLidar, UpdateRadar after setting the respective functions
void UKF::Update(MeasurementPackage meas_package) {
  //mean predicted measurement
  
  VectorXd z_pred = VectorXd(n_z_);
  z_pred.fill(0.0);
  for (unsigned int i=0; i < 2*n_aug_+1; i++) {
      z_pred = z_pred + weights(i) * Zsig.col(i);
  }

  //innovation covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  S.fill(0.0);
  for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
/*     while (z_diff(1)>  M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)< -M_PI) z_diff(1)+=2.*M_PI; */
    z_diff(1) = atan2(cos(z_diff(1)), sin(z_diff(1)) );
    
    
    S = S + weights(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S = S + R;

  
  //////////////////////////////////////////
    //create  vector for incoming radar/lazer measurement
  VectorXd z = VectorXd(n_z_);
  z = meas_package.raw_measurements_;
   /*    5.9214,
      0.2187,
      2.0062; */

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
/*     while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI; */
    z_diff(1) = atan2(cos(z_diff(1)), sin(z_diff(1)) );
    
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
/*     while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI; */
    x_diff(3) = atan2(cos(x_diff(3)), sin(x_diff(3)) );    
    
    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
/*   while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI; */
  z_diff(1) = atan2(cos(z_diff(1)), sin(z_diff(1)) );
  
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
   
  return;
}
