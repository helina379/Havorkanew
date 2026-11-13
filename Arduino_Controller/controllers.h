// controllers.h -- struct definitions used by Arduino_Controller.ino
// Rebuilt from the field names used in the sketch. If mam's original turns up, use that instead.
#ifndef CONTROLLERS_H
#define CONTROLLERS_H

/* ---------------- PID ---------------- */
struct pid_param {
  double Ts;
  double Kp;
  double Ti;
  double Td;
  double beta;
  double alpha;
};

struct pid_calc {
  double setpoint;
  double past_setpoint;
  double procout;
  double past_procout;
  double proportional;
  double integral;
  double next_integral;
  double derivative;
  double past_derivative;
  double der_term;
  double int_term;
  double ctrlsignal_u;
};

/* ---------------- LQR ---------------- */
struct matG  { double G11, G12, G13, G21, G22, G23, G31, G32, G33; };
struct matH  { double H11, H21, H31; };
struct matC  { double C11, C12, C13; };
struct matKe { double Ke_11, Ke_21, Ke_31; };
struct matK1 { double K1_11; };
struct matK2 { double K2_11, K2_12, K2_13; };
struct lqr_states { double x1, x2, x3; };

#endif
