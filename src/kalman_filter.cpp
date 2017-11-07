#include "kalman_filter.h"
#include <math.h> 


using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;  // State vector
  P_ = P_in;  // Object Covariance matrix
  F_ = F_in;  // State Transition matrix
  H_ = H_in;  // Measurement matrix
  R_ = R_in;  // Measurement Covariance matrix
  Q_ = Q_in;  // Process Covariance Matrix
}

void KalmanFilter::Predict() {
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  VectorXd y = z - h(x_);
  NormalizePhi(y);
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

VectorXd KalmanFilter::h(const VectorXd &cart){
  float px = cart(0);
  float py = cart(1);
  float vx = cart(2);
  float vy = cart(3);

  VectorXd polar_out(3);

  float rho = sqrt(px * px + py * py);

  float phi;
  if (fabs(px) < 0.0001 )
    phi = 0.0;
  else
    phi = atan2(py, px);

  float rho_dot;
  if (fabs(rho) < 0.0001 )
    rho_dot = 0.0;
  else
    rho_dot = (px * vx + py * vy) / rho;

  polar_out << rho, phi, rho_dot;

  return polar_out;
}

void KalmanFilter::NormalizePhi(VectorXd &polar_coord){
  float phi = polar_coord(1);

  phi = atan2(sin(phi),cos(phi));

  polar_coord(1) = phi;
}