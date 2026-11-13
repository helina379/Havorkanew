// Hovorka plant (I2C slave) for the master/slave PID setup.
// Replaces Arduino_BergmanModel.ino. Controller sketch stays unchanged.
//
// - Model code is copied from hovorka.ino (same RK4, same Srhb quirk), glucose in mg/dL.
// - Real time: 1 simulated minute = 1 real minute. One RK4 step = 0.1 min = 6 s real.
// - No FreeRTOS on this board (Hovorka's stack use is too big for 128-byte tasks on an Uno).
// - Insulin comes ONLY from the controller (scheduled boluses removed). Meals kept.
// - Logs "time_min,glucose_mgdl,insulin_mUmin" once per simulated minute at 115200 baud.

#include <Wire.h>
#include <EasyTransferI2C.h>
#include <math.h>

// ---------------- I2C / EasyTransfer ----------------
EasyTransferI2C ET;

struct RECEIVE_DATA_STRUCTURE {
  // MUST match the controller sketch exactly
  double request;
  double glucose;
  double insulin;
};
RECEIVE_DATA_STRUCTURE mydata;

#define I2C_SLAVE_ADDRESS  9
#define I2C_MASTER_ADDRESS 2

int flag_glucose = 0;

// ---------------- simulation settings ----------------
const float STEP = 0.1f;                  // min per RK4 step
const unsigned long STEP_MS = 6000UL;     // real ms per step (real time)
const int PRINT_EVERY = 10;               // log once per simulated minute
const float U_MAX = 200.0f;               // clamp on controller insulin (mU/min)

unsigned long nextStep;
long k = 0;

// ---------------- patient ----------------
const float BW      = 70.0f;
const float U_BASAL = 12.9127f;           // mU/min

// ---------------- meals (minutes since start) ----------------
struct Event { float start; float dur; float val; };
const Event MEALS[] = { {210,10,124.17f}, {390,20,34.5f}, {600,10,34.5f}, {765,20,76.95f}, {900,10,34.5f} };
const int N_MEALS = sizeof(MEALS) / sizeof(MEALS[0]);

// ---------------- parameters ----------------
const float tau_d = 40.0f, Ag = 0.8f, Mwg = 180.0f;
const float F01 = 0.00097f * BW;
const float Vg  = 0.16f * BW;
const float k12 = 0.066f;
const float Gb  = 90.0f;
const float Vi  = 0.12f * BW;
const float tau_s = 55.0f;
const float ka1 = 0.006f, ka2 = 0.06f, ka3 = 0.03f;
const float ke  = 0.138f;
const float Sit = 51.2f / 10000.0f, Sid = 8.2f / 10000.0f, Sie = 520.0f / 10000.0f;
const float kb1 = Sit * ka1, kb2 = Sid * ka2, kb3 = Sie * ka3;
const float kp2 = 0.0007f, ke1 = 0.007f;
const float Gth = 162.0f, Gth1 = 60.0f;
const float K6gp = 0.034f, Sc = 297.0f;
const float cth = 8.0f / 1000000.0f;
const float tD = 59.9f, zz = 23.24f;
const float Ggng1b = 0.495f, Ggg1b = 0.7425f;
const float nn = 0.01f;
const float cb = 58.0f / 10000000.0f;
const float rho = 0.86f;
const float sigma = 0.01f / 10000000.0f;
const float S_glyc = 0.98f / 10000000.0f;

// ---------------- state ----------------
float G, Gt, x1, x2, x3, S1, S2, I, Dm1, Dm2;
float c, G6p, EGP, Srhs, Srhd, Srhb, Srh;
float u_in = U_BASAL;   // insulin currently applied (held between controller updates)

static inline float max0(float x) { return x > 0 ? x : 0.0f; }

void initState() {
  const float u = U_BASAL;
  G = 90.0f; Gt = 70.0f;
  S1 = tau_s * u; S2 = tau_s * u;
  I  = u / (0.01656f * BW);
  x1 = 0.30898f * u / BW;
  x2 = 0.04951f * u / BW;
  x3 = 3.2206f  * u / BW;
  const float EGPb = 1.23f, g6po = 5.50f;
  EGP = EGPb;
  G6p = EGPb / K6gp + g6po;
  c = cb;
  Srhs = 0; Srhd = 0; Srhb = 0; Srh = 0;
  Dm1 = 0; Dm2 = 0;
}

void stepModel(float t, float d_cho, float u) {
  const float h = STEP;

  float D_meal = 1000.0f * d_cho / Mwg;

  // ---- meal absorption ----
  float rkm11 = h * (Ag * D_meal - Dm1 / tau_d);
  float rkm21 = h * (Dm1 / tau_d - Dm2 / tau_d);
  float rkm12 = h * (Ag * D_meal - (Dm1 + rkm11 / 2) / tau_d);
  float rkm22 = h * ((Dm1 + rkm11 / 2) / tau_d - (Dm2 + rkm21 / 2) / tau_d);
  float rkm13 = h * (Ag * D_meal - (Dm1 + rkm12 / 2) / tau_d);
  float rkm23 = h * ((Dm1 + rkm12 / 2) / tau_d - (Dm2 + rkm22 / 2) / tau_d);
  float rkm14 = h * (Ag * D_meal - (Dm1 + rkm13) / tau_d);
  float rkm24 = h * ((Dm1 + rkm13) / tau_d - (Dm2 + rkm23) / tau_d);

  float Ug = Dm2 / tau_d;  // pre-update Dm2
  Dm1 += (rkm11 + 2 * rkm12 + 2 * rkm13 + rkm14) / 6;
  Dm2 += (rkm21 + 2 * rkm22 + 2 * rkm23 + rkm24) / 6;

  // ---- glucose subsystem terms (fixed across the 4 stages) ----
  float Ugc  = 18 * Ug / Vg;
  float EGPc = 18 * (EGP / Vg) * (1 - x3) - kp2 * (G - Gb);

  float Erc = 0;
  if (G >= Gth) Erc = ke1 * (G - Gth);

  float F01uc = 18 * F01 * G / (Vg * 81);
  if (G >= 81) F01uc = 18 * F01 / Vg;

  float G0 = G, Gt0 = Gt;

  float rkg11 = h * (Ugc - F01uc - Erc + k12 * (Gt0 - Gb) - x1 * (G0 - Gb) + EGPc);
  float rkg21 = h * (x1 * (G0 - Gb) - (k12 + x2) * (Gt0 - Gb));

  float rkg12 = h * (Ugc - F01uc - Erc + k12 * ((Gt0 + rkg21 / 2) - Gb) - x1 * ((G0 + rkg11 / 2) - Gb) + EGPc);
  float rkg22 = h * (x1 * ((G0 + rkg11 / 2) - Gb) - (k12 + x2) * ((Gt0 + rkg21 / 2) - Gb));

  float rkg13 = h * (Ugc - F01uc - Erc + k12 * ((Gt0 + rkg22 / 2) - Gb) - x1 * ((G0 + rkg12 / 2) - Gb) + EGPc);
  float rkg23 = h * (x1 * ((G0 + rkg12 / 2) - Gb) - (k12 + x2) * ((Gt0 + rkg22 / 2) - Gb));

  float rkg14 = h * (Ugc - F01uc - Erc + k12 * ((Gt0 + rkg23) - Gb) - x1 * ((G0 + rkg13) - Gb) + EGPc);
  float rkg24 = h * (x1 * ((G0 + rkg13) - Gb) - (k12 + x2) * ((Gt0 + rkg23) - Gb));

  // ---- glycogenolysis drive ----
  float E   = (1 - tanhf((t - tD) / zz)) / 2;
  float Ggg = Ggg1b + Sc * max0(c - cth) * E;

  float dGt_scaled = (rkg21 + 2 * rkg22 + 2 * rkg23 + rkg24) / 6 * h;
  float new_EGP;
  if (dGt_scaled >= 0) new_EGP = K6gp * G6p - x3 * dGt_scaled - kp2 * (G - Gb);
  else                 new_EGP = K6gp * G6p - kp2 * (G - Gb);

  // ---- Srh (bug-for-bug, same as model.py) ----
  if (G >= Gb) Srhs = rho * (Srhs - Srhb);
  else         Srhb = nn * cb;
  Srhs = rho * (Srhs - max0(sigma * (Gth1 - G) / (I + 1) + Srhb));

  Srhd = S_glyc * max0(-dGt_scaled);
  Srh  = Srhs + Srhd;

  float rkc1 = h * (-nn * c + Srh);
  float rkc2 = h * (-nn * (c + rkc1 / 2) + Srh);
  float rkc3 = h * (-nn * (c + rkc2 / 2) + Srh);
  float rkc4 = h * (-nn * (c + rkc3) + Srh);
  c += (rkc1 + 2 * rkc2 + 2 * rkc3 + rkc4) / 6;

  G  += (rkg11 + 2 * rkg12 + 2 * rkg13 + rkg14) / 6;
  Gt += (rkg21 + 2 * rkg22 + 2 * rkg23 + rkg24) / 6;

  float rkg61 = h * (-K6gp * G6p + Ggg + Ggng1b);
  float rkg62 = h * (-K6gp * (G6p + rkg61 / 2) + Ggg + Ggng1b);
  float rkg63 = h * (-K6gp * (G6p + rkg62 / 2) + Ggg + Ggng1b);
  float rkg64 = h * (-K6gp * (G6p + rkg63) + Ggg + Ggng1b);
  G6p += (rkg61 + 2 * rkg62 + 2 * rkg63 + rkg64) / 6;

  EGP = new_EGP;

  // ---- insulin absorption ----
  float s1 = S1, s2 = S2;
  float rki11 = h * (u - s1 / tau_s);
  float rki21 = h * ((s1 - s2) / tau_s);
  float rki12 = h * (u - (s1 + rki11 / 2) / tau_s);
  float rki22 = h * (((s1 + rki11 / 2) - (s2 + rki21 / 2)) / tau_s);
  float rki13 = h * (u - (s1 + rki12 / 2) / tau_s);
  float rki23 = h * (((s1 + rki12 / 2) - (s2 + rki22 / 2)) / tau_s);
  float rki14 = h * (u - (s1 + rki13) / tau_s);
  float rki24 = h * (((s1 + rki13) - (s2 + rki23)) / tau_s);

  float Ui = S2 / tau_s;
  S1 += (rki11 + 2 * rki12 + 2 * rki13 + rki14) / 6;
  S2 += (rki21 + 2 * rki22 + 2 * rki23 + rki24) / 6;

  // ---- insulin action ----
  float I0 = I;
  float rkx11 = h * (-ka1 * x1 + kb1 * I0), rkx21 = h * (-ka2 * x2 + kb2 * I0), rkx31 = h * (-ka3 * x3 + kb3 * I0);
  float rkx12 = h * (-ka1 * (x1 + rkx11 / 2) + kb1 * I0), rkx22 = h * (-ka2 * (x2 + rkx21 / 2) + kb2 * I0), rkx32 = h * (-ka3 * (x3 + rkx31 / 2) + kb3 * I0);
  float rkx13 = h * (-ka1 * (x1 + rkx12 / 2) + kb1 * I0), rkx23 = h * (-ka2 * (x2 + rkx22 / 2) + kb2 * I0), rkx33 = h * (-ka3 * (x3 + rkx32 / 2) + kb3 * I0);
  float rkx14 = h * (-ka1 * (x1 + rkx13) + kb1 * I0),     rkx24 = h * (-ka2 * (x2 + rkx23) + kb2 * I0),     rkx34 = h * (-ka3 * (x3 + rkx33) + kb3 * I0);

  x1 += (rkx11 + 2 * rkx12 + 2 * rkx13 + rkx14) / 6;
  x2 += (rkx21 + 2 * rkx22 + 2 * rkx23 + rkx24) / 6;
  x3 += (rkx31 + 2 * rkx32 + 2 * rkx33 + rkx34) / 6;

  // ---- plasma insulin ----
  float rkia1 = h * (Ui / Vi - ke * I0);
  float rkia2 = h * (Ui / Vi - ke * (I0 + rkia1 / 2));
  float rkia3 = h * (Ui / Vi - ke * (I0 + rkia2 / 2));
  float rkia4 = h * (Ui / Vi - ke * (I0 + rkia3));
  I += (rkia1 + 2 * rkia2 + 2 * rkia3 + rkia4) / 6;
}

float mealRate(float t) {
  for (int i = 0; i < N_MEALS; i++)
    if (t >= MEALS[i].start && t <= MEALS[i].start + MEALS[i].dur) return MEALS[i].val / MEALS[i].dur;
  return 0.0f;
}

void receive(int numBytes) {}   // empty on purpose, same as the Bergman sketch

void setup() {
  Serial.begin(115200);
  Wire.begin(I2C_SLAVE_ADDRESS);
  ET.begin(details(mydata), &Wire);
  Wire.onReceive(receive);

  mydata.request = 0;
  mydata.glucose = 90.0;
  mydata.insulin = U_BASAL;

  initState();
  Serial.println(F("time_min,glucose_mgdl,insulin_mUmin"));
  nextStep = millis();
}

void loop() {
  // ---- I2C: poll for a request or an insulin value from the controller ----
  if (ET.receiveData()) {
    if (mydata.request == 1) {
      mydata.request = 0;
      flag_glucose = 1;
    } else {
      float u = (float)mydata.insulin;
      if (u < 0) u = 0;               // pump can't deliver negative insulin
      if (u > U_MAX) u = U_MAX;
      u_in = u;
    }
  }

  if (flag_glucose == 1) {
    mydata.glucose = G;
    mydata.insulin = u_in;
    ET.sendData(I2C_MASTER_ADDRESS);
    flag_glucose = 0;
  }

  // ---- model: one RK4 step every STEP_MS of real time ----
  if ((long)(millis() - nextStep) >= 0) {
    nextStep += STEP_MS;
    float t = k * STEP;
    if (k % PRINT_EVERY == 0) {
      Serial.print(t, 1);
      Serial.print(',');
      Serial.print(G, 3);
      Serial.print(',');
      Serial.println(u_in, 3);
    }
    stepModel(t, mealRate(t), u_in);
    k++;
  }
}
