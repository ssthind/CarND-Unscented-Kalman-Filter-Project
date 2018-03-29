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
  //for running initializing code in ProcessMeasurement fuction
  is_initialized_ = false;
  
  //// Parameter setting
  n_x_ = 5;
  n_aug_ = n_x_ + 2;
  n_z_laser = 2;
  n_z_radar = 3;
  lambda_ = 3-n_aug_; 
  
  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);
 
  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1; // 1; // 3;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.9; //1; //0.3;
  
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
  
  
  P_ = MatrixXd(n_x_,n_x_);
			
	P_ <<  0.1 , 0, 0, 0, 0,
			0, 0.1, 0, 0, 0,
			0, 0, 1, 0, 0,
			0, 0, 0, 1, 0,
			0, 0, 0, 0, 1 ;
   
  ///NIS check initialize
  NIS_calc = 0.0;
  NIS_chk_lazer = 0;
  NIS_chk_radar = 0;
  counter_reading = 0;
}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  Complete this function! Switches between lidar and radar
  measurements.
  */
  if (!is_initialized_) {
    x_= VectorXd(n_x_);
    previous_timestamp_ = meas_package.timestamp_;
        
    x_ << 1,1,1,1,0.1;    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float rho = meas_package.raw_measurements_[0];
      float phi = meas_package.raw_measurements_[1];

	  x_(0) = rho*cos(phi);
	  x_(1) = rho*sin(phi);
      //cout << "initial x_ in radar= " << x_ << endl;
      is_initialized_ = true;
      return;	  
	  
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      
      /* Initialize state. */
	  //set the state with the initial location and zero velocity
      x_(0) = meas_package.raw_measurements_[0];
	  x_(1) = meas_package.raw_measurements_[1];
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
  //Calling UKF predict function for each iteration.
  Prediction(delta_t);

  //Setting up log files for NIS
  ofstream logfile;
  logfile.open ("NIS_log.csv", ios::out | ios::app );
  
  //Calling measurement
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR and use_radar_ == true) {
    UpdateRadar(meas_package);
	//cout << "NIS Radar: "<< NIS_calc << endl;
	logfile << "NIS Radar , "  << NIS_calc  << endl;
	if ( NIS_calc > 7.815){   ///  value below 95% for degree of freedom 3
		NIS_chk_radar +=1; 
	}

	
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER  and use_laser_==true) {  
    UpdateLidar(meas_package);
	//cout << "NIS Laser: "<< NIS_calc << endl;
    logfile << "NIS Laser , "  << NIS_calc  << endl;
	if ( NIS_calc > 5.991){    ///  value below 95% for degree of freedom 2
		NIS_chk_lazer +=1; 
	}
  }
  //closing log file for NIS
  logfile.close();
  // updates NIS related calculation
  counter_reading +=1;
  if (counter_reading ==498){
	cout << "Percentage NIS:  " << 100*(1- (NIS_chk_lazer + NIS_chk_radar)/counter_reading) << "%" << endl;
	cout << "Percentage NIS Radar:  " << 100*(1- (NIS_chk_radar)*2/counter_reading) << "%" << endl;
	cout << "Percentage NIS Laser:  " << 100*(1- (NIS_chk_lazer)*2/counter_reading) << "%" << endl;
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
  Estimated the object's location. Modify the state
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

  //create square root matrix of P augmented matrix
  MatrixXd A = P_aug.llt().matrixL();
  
  //create augmented sigma points
  float lamda_n_x_root = sqrt(lambda_+n_aug_);
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2*n_aug_+1);
  Xsig_aug.col(0) << x_aug;
  for (unsigned int i=0; i < n_aug_ ; i++){
    Xsig_aug.col(i+1) << x_aug + lamda_n_x_root * A.col(i);
    Xsig_aug.col(i+1+n_aug_) << x_aug - lamda_n_x_root * A.col(i);
      
  }

  // apply prediction function on sigma points 
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_+1);
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
    if (fabs(yawd) > 0.00001) {
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
 
  //predicted state mean
  VectorXd x_mean = VectorXd(n_x_);
  x_mean.fill(0.0);
  for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points
    x_mean = x_mean + weights(i) * Xsig_pred_.col(i);
  }
  
  //predicted state covariance matrix
  MatrixXd P = MatrixXd(n_x_,n_x_);
  VectorXd x_diff;
  P.fill(0.0);
  for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {  //iterate over sigma points

    // state difference
    x_diff = Xsig_pred_.col(i) - x_mean;
    //angle normalization
    x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3)) );

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
  Used lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  calculate the lidar NIS in common NIS function.
  */
  //updating number of measurement points i.e 2 for lider/laser
  n_z_ = n_z_laser;
  //setting R for laser measurement noise
  R = R_laser;
 
  //transform sigma points into measurement space would be direct mapping of x,y sigma points 
  Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  //Zsig.row(0) = Xsig_pred_.row(0); /// x measurement
  //Zsig.row(1) = Xsig_pred_.row(1); /// y measurement  
  Zsig = Xsig_pred_.block(0,0, n_z_, 2 * n_aug_ + 1);
  Update(meas_package);

  return;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  Used radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.
  calculate the radar NIS in common NIS function.
  */
  
  ///updating number of measurement points i.e 3 for radar
  n_z_ = n_z_radar;
  ////setting R for radar measurement noise
  R = R_radar;

  //create matrix for sigma points in measurement space
  Zsig = MatrixXd(n_z_, 2 * n_aug_ + 1);
  
    //transform sigma points into measurement space
  for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    // converting cartesian to polar coordinates for differencing with Radar measurements  in polar coordinates
   double rho = sqrt(Xsig_pred_.col(i)[0]*Xsig_pred_.col(i)[0] + Xsig_pred_.col(i)[1]*Xsig_pred_.col(i)[1]);
   double psi = atan2(Xsig_pred_.col(i)[1],Xsig_pred_.col(i)[0]);
   double phi = Xsig_pred_.col(i)[3];
   double rhodot = Xsig_pred_.col(i)[2]*(Xsig_pred_.col(i)[0]*cos(phi) + Xsig_pred_.col(i)[1]*sin(phi))/rho;
   // measurement model
   Zsig.col(i) <<  rho, psi, rhodot;

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

  //innovation/Predicted measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z_,n_z_);
  S.fill(0.0);
  for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    //angle normalization
	if (n_z_==3){
		z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)) );
    }
    S = S + weights(i) * z_diff * z_diff.transpose();
  }
  //add measurement noise covariance matrix
  S = S + R;
  
  //create  vector for incoming radar/lazer measurement
  VectorXd z = VectorXd(n_z_);
  z = meas_package.raw_measurements_;

  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (unsigned int i = 0; i < 2 * n_aug_ + 1; i++) {  //2n+1 simga points

    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    if (n_z_==3){
	  z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)) );
    }
    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    //angle normalization
    x_diff(3) = atan2(sin(x_diff(3)), cos(x_diff(3)) );    
    
    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }

  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //residual
  VectorXd z_diff = z - z_pred;

  //angle normalization
  if (n_z_==3){
	z_diff(1) = atan2(sin(z_diff(1)), cos(z_diff(1)) );
  }
  //update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K*S*K.transpose();
  
  
  //NIS calcualtion for laser/radar update step
  NIS_calc = z_diff.transpose() * S.inverse() * z_diff;
  
  return;
}
