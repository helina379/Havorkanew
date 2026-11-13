#include <Arduino_FreeRTOS.h>
#include <stdio.h>
#include <Wire.h>
#include <EasyTransferI2C.h>
#include "controllers.h"

void TaskOne(void *pvParameters);
void TaskTwo(void *pvParameters);

/**********************************/
/********* Easy Transfer *********/
/*********************************/

//create object
EasyTransferI2C ET; 

struct SEND_DATA_STRUCTURE{
  //put your variable definitions here for the data you want to send
  //THIS MUST BE EXACTLY THE SAME ON THE OTHER ARDUINO
  double request;
  double glucose;
  double insulin;
};

//give a name to the group of data
SEND_DATA_STRUCTURE mydata;

//define slave and master i2c address
#define I2C_SLAVE_ADDRESS 9
#define I2C_MASTER_ADDRESS 2

/*********************************/
/**** Definiciones Generales ****/
/********************************/
int data;
double ctrlsig_u;  
double insulin;
int flag_insulin;

/************************************/
/********* Controlador PID *********/
/***********************************/

// Definicion de estructuras para el controlador
struct pid_param params;
struct pid_calc calc;

// Funcion del controlador PID
double pidController(struct pid_param *params, struct pid_calc *calc)
{
  double ctrlsig_u;

  // Calculo de valores actuales
  calc->proportional = params->Kp * (params->beta * calc->setpoint - calc->procout);
  calc->der_term = (params->alpha * params->Td) / (params->alpha * params->Td + params->Ts);
  calc->derivative = calc->der_term * calc->past_derivative - calc->der_term * (calc->procout - calc->past_procout);
  calc->ctrlsignal_u = 48.2328 + calc->proportional + calc->integral + calc->derivative; // salida del controlador (48.2328 +)
//
//  Serial.print("Proportional: ");
//  Serial.print(calc->proportional);
//  Serial.print(" , Integral: ");
//  Serial.print(calc->integral);
//  Serial.print(" , Derivative: ");
//  Serial.print(calc->derivative);
//  Serial.print(" , Ctrl_sig: ");
//  Serial.print(calc->ctrlsignal_u);
//  Serial.print("\n");
  
    // Calculo de valores anteriores
  calc->past_derivative = calc->derivative;
  calc->past_setpoint = calc->setpoint;
  calc->past_procout = calc->procout;
  calc->integral = calc->next_integral;

  // Calculo de valores siguientes
  calc->int_term = params->Ts * params->Kp / params->Ti;
  calc->next_integral = calc->integral + calc->int_term * (calc->setpoint - calc->procout);

  ctrlsig_u = calc->ctrlsignal_u;
  return ctrlsig_u;
};

/************************************/
/********* Controlador LQR *********/
/***********************************/
  //Definicion de matrices del controlador
  struct matG G;
  struct matH H;
  struct matC C;
  struct matKe Ke;
  struct matK1 K1;
  struct matK2 K2;
  struct lqr_states currentStates;

  //Definicion de otros parametros 
  double r_act;
  double I_past;
  double y_act;

// Funcion del controlador LQR
double lqrController(
  struct matG *G,
  struct matH *H,
  struct matC *C,
  struct matKe *Ke,
  struct matK1 *K1,
  struct matK2 *K2,
  struct lqr_states *x,
  double r_act,
  double y_act,
  double *I_past)
{
  double I_act;
  double ctrlsig_u;
  double x1_next;
  double x2_next;
  double x3_next;

  // Valores actuales
  I_act = *I_past + r_act - y_act;
  ctrlsig_u = K1->K1_11 * I_act - (K2->K2_11 * x->x1 + K2->K2_12 * x->x2 + K2->K2_13 * x->x3) + 48.2328;

  // Guardar valor actual
  *I_past = I_act;

  // Valores siguientes
  x1_next = G->G11 * x->x1 + G->G12 * x->x2 + G->G13 * x->x3 + H->H11 * ctrlsig_u + Ke->Ke_11 * (y_act - C->C11 * x->x1 - C->C12 * x->x2 - C->C13 * x->x3);
  x2_next = G->G21 * x->x1 + G->G22 * x->x2 + G->G23 * x->x3 + H->H21 * ctrlsig_u + Ke->Ke_21 * (y_act - C->C11 * x->x1 - C->C12 * x->x2 - C->C13 * x->x3);
  x3_next = G->G31 * x->x1 + G->G32 * x->x2 + G->G33 * x->x3 + H->H31 * ctrlsig_u + Ke->Ke_31 * (y_act - C->C11 * x->x1 - C->C12 * x->x2 - C->C13 * x->x3);

  return ctrlsig_u;
}

/*----------------------------------------------------------------------------------------*/
/*-------------------------------------Set up---------------------------------------------*/
/*----------------------------------------------------------------------------------------*/
void setup() {
  Serial.begin(9600);
  Wire.begin(I2C_MASTER_ADDRESS);
  //start the library, pass in the data details and the name of the serial port. Can be Serial, Serial1, Serial2, etc.
  ET.begin(details(mydata), &Wire);  
  //define handler function on receiving data
  Wire.onReceive(receive);

  flag_insulin = 0;
  mydata.insulin = 0;
  mydata.request = 0;

  /********************************/
  /**** Setup Controlador PID ****/
  /*******************************/
  
  // Inicializacion de los parametros del controlador PID
  params = { 0 };
  calc = { 0 };
  
  params.Ts = 1;
  params.Kp = -1.228;
  params.Ti = 81.26;
  params.Td = 19.5;
  params.beta = 1;
  params.alpha = 0.01;

  //calc.past_procout = 0.0;
  calc.past_procout = 120.0;
  calc.setpoint = calc.past_setpoint = 120.0;

  /********************************/
  /**** Setup Controlador LQR ****/
  /*******************************/
  
  // Inicializacion de las matrices del controlador
  G.G11 = 0.9523;
  G.G12 = 0.0302;
  G.G13 = 0.0008;
  G.G21 = -0.1906;
  G.G22 = 0.8342;
  G.G23 = 0.0658;
  G.G31 = -0.2211;
  G.G32 = -0.5331;
  G.G33 = 0.3039;

  H.H11 = 0.0001;
  H.H21 = -0.0130;
  H.H31 = -0.1167;

  C.C11 = 9.1569;
  C.C12 = 0.1399;
  C.C13 = 0.0044;

  Ke.Ke_11 = 0.0398;
  Ke.Ke_21 = -0.0135;
  Ke.Ke_31 = -0.0076;

  K1.K1_11 = -0.0182;

  K2.K2_11 = -1.9176;
  K2.K2_12 = -0.2889;
  K2.K2_13 = -0.0314;

  // Inicializacion de los estados iniciales
  currentStates.x1 = 0.0;
  currentStates.x2 = 0.0;
  currentStates.x3 = 0.0;

  // Inicializacion de otras variables
  r_act = 120.0;
  //y_act = 0.0;
  y_act = 120.0;
  I_past = 0.0;
  ctrlsig_u = 0.0;
  ctrlsig_u = 48.238;

      xTaskCreate(
    TaskOne
    , "TaskOne"
    , 128 // This stack size can be checked & adjusted by reading Highwater
    , NULL
    , 1  // priority
    , NULL);

    xTaskCreate(
    TaskTwo
    , "TaskTwo"
    , 128 // This stack size can be checked & adjusted by reading Highwater
    , NULL
    , 1  // priority
    , NULL);

}

void loop() {

}

void TaskOne(void *pvParameters)  // This is a task.
{
  (void)pvParameters;
  TickType_t xLastWakeTime;
  const TickType_t xFrequency = (60000/portTICK_PERIOD_MS);

  // Initialise the xLastWakeTime variable with the current time.
  xLastWakeTime = xTaskGetTickCount();

  for (;;)
  {
    //Serial.println("Ask for glucose");
    mydata.request = 1;
    ET.sendData(I2C_SLAVE_ADDRESS);
    mydata.request = 0;
    vTaskDelayUntil( &xLastWakeTime, xFrequency );
  }
}

void TaskTwo(void *pvParameters)  // This is a task.
{
  (void)pvParameters;
  TickType_t xLastWakeTime;
  const TickType_t xFrequency = (20000/portTICK_PERIOD_MS);

  // Initialise the xLastWakeTime variable with the current time.
  xLastWakeTime = xTaskGetTickCount();

  for (;;)
  {
    if(ET.receiveData()){ 
    //Serial.println("Receiving glucose");
    //Serial.println(mydata.glucose);
    y_act = mydata.glucose;
    calc.procout = mydata.glucose;
    flag_insulin = 1;
  }
    if(flag_insulin == 1){
      flag_insulin = 0;
      //ctrlsig_u = lqrController(&G, &H, &C, &Ke, &K1, &K2, &currentStates, r_act, y_act, &I_past);
      ctrlsig_u = pidController(&params, &calc);
      mydata.insulin = ctrlsig_u;
      //Serial.println("Sending Insulin");
      Serial.println(mydata.insulin);
      mydata.request = 0;
      ET.sendData(I2C_SLAVE_ADDRESS);
      }
    vTaskDelayUntil( &xLastWakeTime, xFrequency );
  }
}
/*----------------------------------------------------------------------------------------*/
/*-----------------------------Recibir dato de glucosa------------------------------------*/
/*----------------------------------------------------------------------------------------*/
//
void receive(int numBytes) {
}
