#include "tobu.hpp"
/* -----------------------------------------------------------------------------------------------------------------------*/
//マルチコア設定
semaphore_t sem;

//for Kalman filter
Matrix<float, 7 ,1> Xp = MatrixXf::Zero(7,1);
Matrix<float, 7 ,1> Xe = MatrixXf::Zero(7,1);
Matrix<float, 7 ,7> P = MatrixXf::Identity(7,7);
Matrix<float, 6 ,1> Z = MatrixXf::Zero(6,1);
Matrix<float, 3, 1> Omega_m = MatrixXf::Zero(3, 1);
Matrix<float, 6, 6> Q = MatrixXf::Identity(6, 6)*1;
Matrix<float, 6, 6> R = MatrixXf::Identity(6, 6)*1;
Matrix<float, 7 ,6> G;
Matrix<float, 3 ,1> Beta;
volatile float Ax,Ay,Az,Wp,Wq,Wr,Mx,My,Mz;
volatile float Dmx,Dmy,Dmz;
volatile float Wqa=0.0, Wpa=0.0,Wra=0.0; 
volatile float Phi=0.0,Theta=0.0,Psi=0.0;
volatile float Kalman_time=0.0;
volatile uint8_t Hzcount=0;
volatile float Psiav=0.0, Phiav=0.0, Thetaav=0.0;
//for Log
volatile uint32_t Logcount=0;
volatile uint32_t Printcount=0;
const uint16_t Logdatanum=45000;
volatile float Logdata[Logdatanum]={0.0};
char sbuf[500];
const uint8_t DATANUM=33;
volatile float Tlog = 0.0;
volatile float Led_timer = 0.0;
volatile uint8_t Led_flag = 1;

//for Control
volatile float Ref_phi, Ref_t, Ref_psi;
volatile float Ref_p,Ref_q,Ref_r;
volatile float Olderr1_p=0.0,Olderr2_p=0.0,Olderr3_p=0.0;
volatile float Olderr1_q=0.0,Olderr2_q=0.0,Olderr3_q=0.0;
volatile float Olderr1_r=0.0,Olderr2_r=0.0,Olderr3_r=0.0;
volatile float Olderr1_phi=0.0, Olderr2_phi=0.0, Olderr3_phi=0.0;
volatile float Olderr1_t=0.0,   Olderr2_t=0.0,   Olderr3_t=0.0;
volatile float Olderr1_psi=0.0;
volatile float Sk_p=0.0, Sk_q=0.0, Sk_r=0.0;
volatile float Sk_phi=0.0, Sk_t=0.0, Sk_psi=0.0;
volatile float Dk_p=0.0, Dk_q=0.0, Dk_r=0.0;
volatile float Dk_phi=0.0, Dk_t=0.0;
volatile float Up=0.0,Uq=0.0,Ur=0.0;

//Filter Time Constance
const float Tc_angl = 0.018;
const float Tc_rate = 0.012;

//PID Gain
const float Kp_phi   = 4.0;//3.5
const float Ti_phi   = 3000.0;//5000
const float Td_phi = 0.002;//0.002;

const float Kp_theta = 3.2;
const float Ti_theta = 5000.0;
const float Td_theta = 0.002;//0.00004;

const float Kp_psi = 0.7;

const float Kp_p = 0.2;//0.1;
const float Ti_p = 50000.0;
const float Td_p = 0.0;//0.0000014;

const float Kp_q = 0.29;//0.1;
const float Ti_q = 40000.0;
const float Td_q = 0.0;//0.0000014;

const float Kp_r = 0.72;//0.1;
const float Ti_r = 50000.0;
const float Td_r = 0.0;//0.0000014;

//for Motor control
volatile float Duty_rr,Duty_rl,Duty_fl,Duty_fr; 
volatile float Com_rr,Com_rl,Com_fl,Com_fr; 

//Rocking wings
float Rocking_time=0.0;
float Rocking_flag = 0;


//for etc
const uint LED_PIN = 25;
volatile uint8_t Safetycount=1;
uint32_t s_time, e_time, d_time;
volatile uint8_t Kalman_flag=1;
float Q0av=0.0;


#if 0
//volatile float PHI,THETA,PSI;
//volatile float Spsi=0.0,Sphi=0.0,St=0.0;

volatile float data2MID,data4MID;
volatile float Time_for_debug = 0.0;
#endif

//ローカル関数宣言
void read_sensor(void);
void kalman(void);
/*--------------------------------------------------------------------------------------------------------------------------------------*/

void kalman(void){

  float err_psi, err_phi, err_t;
  float dk_a,dk_e,dk_r;
  float dk_phi,dk_psi,dk_t;
  float dt=0.01;

  /*------------------------------------------------------------------------------------------------------------------------------------*/
  while(1)
  {
    sem_acquire_blocking(&sem);
    sem_reset(&sem, 0);
    //printf("%f\n",time);

    //s_time=time_us_32();
    
    dt=0.01;
    Omega_m << (float)Wp, (float)Wq, (float)Wr;
    Z << (float)Ax, (float)Ay, (float)Az, (float)Mx, (float)My, (float)Mz;//ここに入れる
    //--Begin Extended Kalman Filter--
    if(Kalman_flag==1)
    {
      //std::cout << "Omega" << Omega_m << std::endl;
      //std::cout << "Z" << Z << std::endl;
      //printf("%8.2f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
      //  Kalman_time, Xe(0,0),Xe(1,0),Xe(2,0),Xe(3,0),Xe(4,0),Xe(5,0),Xe(6,0));
      
      ekf(Xp, Xe, P, Z, Omega_m, Q, R, G*dt, Beta, dt);

      Phi = CalcPhi(Xe);
      Theta = CalcTheta(Xe);
      Psi = CalcPsi(Xe);
    
      Kalman_time = Kalman_time + 0.01;
    }

    //スティック上
    if(Data6<-0.2)
    {
      //printf("upper log sw\n");
      if(Logcount < Logdatanum-DATANUM){
        // LED Blink
        if(Led_timer < 20){
          gpio_put(LED_PIN, Led_flag);
          Led_timer++;
        }
        else
        {
          Led_flag=1-Led_flag;
          Led_timer = 0;
        }
        // Log write
        Logdata[Logcount++]=Xe(0,0);
        Logdata[Logcount++]=Xe(1,0);
        Logdata[Logcount++]=Xe(2,0);
        Logdata[Logcount++]=Xe(3,0);
        Logdata[Logcount++]=Xe(4,0);

        Logdata[Logcount++]=Xe(5,0);
        Logdata[Logcount++]=Xe(6,0);
        Logdata[Logcount++]=Wp;
        Logdata[Logcount++]=Wq;
        Logdata[Logcount++]=Wr;
        
        Logdata[Logcount++]=Ax;
        Logdata[Logcount++]=Ay;
        Logdata[Logcount++]=Az;
        Logdata[Logcount++]=Mx;
        Logdata[Logcount++]=My;
        
        Logdata[Logcount++]=Mz;
        Logdata[Logcount++]=Ref_p;
        Logdata[Logcount++]=Ref_q;
        Logdata[Logcount++]=Ref_r;
        Logdata[Logcount++]=Phi-Phiav;
        
        Logdata[Logcount++]=Theta-Thetaav;
        Logdata[Logcount++]=Psi-Psiav;
        Logdata[Logcount++]=Ref_phi;
        Logdata[Logcount++]=Ref_t;
        Logdata[Logcount++]=Ref_psi;
        
        Logdata[Logcount++]=Up;
        Logdata[Logcount++]=Uq;
        Logdata[Logcount++]=Ur;
        Logdata[Logcount++]=Dk_p;
        Logdata[Logcount++]=Dk_q;
        
        Logdata[Logcount++]=Dk_r;
        Logdata[Logcount++]=Dk_phi;
        Logdata[Logcount++]=Dk_t;
      }
      else{
        gpio_put(LED_PIN,0);
      }
    }

    //スティック下
    else if(Data6>0.2){
      //printf("lower log sw\n");
      if(Data3 < 0.3){
        if(Printcount+DATANUM < Logdatanum ){
          for (uint8_t i=0;i<DATANUM;i++){
            if(i==0){
              printf("%8.2f ", Tlog);
              Tlog = Tlog + 0.01;
            }
            sprintf(sbuf,"%16.5f ",Logdata[Printcount+i]);
            printf("%s",sbuf);
          }

          printf("\n");
          Printcount=Printcount+DATANUM;
        }
        else
        {
          gpio_put(LED_PIN, 0);
        }
      }
    }
    else if(-0.1<Data6 && Data6<0.1)
    {
      //printf("mid log sw\n");
      if(Kalman_time>12.0)gpio_put(LED_PIN, 0);
      Printcount = 0;
      Logcount = 0;
      Tlog = 0.0;
    }
#if 1
    //角度PID制御

    Ref_phi = (Data4)*0.52398775598299*2;
    Ref_t =   (Data2)*0.523598775598299*2;
    Ref_psi =  Data1;

    if(Data3>=0.25){
      //phi
      err_phi = (Ref_phi - (Phi - Phiav) );
      Sk_phi = Sk_phi + err_phi;//修正しました
      if (Sk_phi>30000.0)         //修正しました
      {                         //修正しました
        Sk_phi = 30000.0;         //修正しました
      }                         //修正しました
      else if(Sk_phi<-30000.0)    //修正しました
      {                         //修正しました
        Sk_phi=-30000.0;          //修正しました
      }                         //修正しました
      dk_phi = (err_phi - Olderr3_phi)*100.0;
      //if(dk_phi > 100) dk_phi = 100;
      //else if(dk_phi < -100) dk_phi = -100;
      Dk_phi = Dk_phi*Tc_angl/(Tc_angl + 0.01) + dk_phi*0.01/(Tc_angl + 0.01);
      Ref_p = Kp_phi*(err_phi + Sk_phi*0.01/Ti_phi + Td_phi * Dk_phi);
      if(Ref_p>36.0)Ref_p = 36.0;
      else if(Ref_p<-36.0)Ref_p =-36.0;
      Olderr3_phi = Olderr2_phi;
      Olderr2_phi = Olderr1_phi;
      Olderr1_phi = err_phi;
      
      //theta
      err_t = (Ref_t - (Theta - Thetaav) );
      Sk_t = Sk_t + err_t;//修正しました
      if (Sk_t>30000.0)     //修正しました
      {                   //修正しました
        Sk_t = 30000.0;     //修正しました
      }                   //修正しました
      else if(Sk_t<-30000.0)//修正しました
      {                   //修正しました
        Sk_t =-30000.0;     //修正しました
      }                   //修正しました
      dk_t = (err_t - Olderr3_t)*100.0;
      Dk_t = Dk_t * Tc_angl/(Tc_angl + 0.01) + dk_t * 0.01/(Tc_angl + 0.01);
      Ref_q = Kp_theta * (err_t + Sk_t*0.01/Ti_theta + Td_theta * Dk_t);
      if(Ref_q>36.0)Ref_q = 36.0;
      else if(Ref_q<-36.0)Ref_q =-36.0;
      Olderr3_t = Olderr2_t;
      Olderr2_t = Olderr1_t;
      Olderr1_t = err_t;
      
      //Psi
      Ref_r = Ref_psi;

    }
    else{
      Sk_phi=0.0;
      Sk_t=0.0;
      //Sphi=0.0;
      //St=0.0;
      Ref_p = 0.0;
      Ref_q = 0.0;
      Ref_r = 0.0;
      Olderr1_phi=0.0;
      Olderr2_phi=0.0;
      Olderr3_phi=0.0;
      Olderr1_t=0.0;
      Olderr2_t=0.0;
      Olderr3_t=0.0;
      //現在の角度（水平の姿勢）を記憶
      Phiav=Phi;
      Thetaav=Theta;
      //data2MID=Data2;
      //data4MID=Data4;
      Dk_phi = 0.0;//追加しました
      Dk_t = 0.0;  //追加しました
    }


    //printf("kalman %04f\n", Data6);    

#endif
  
    //e_time=time_us_32();
    //d_time=e_time-s_time;
    //printf("%u\n", d_time);
  }
}


/*--------------------------------------------------------------------------------------------------------------------------------------*/
void MAINLOOP(void)
{
  float mx1,my1,mz1,mag_norm;
  float dk_p, dk_q, dk_r;
  float err_p, err_q, err_r;
  float p_rate, q_rate, r_rate;
  
  //e エレベータ a エルロン r ラダー t スロットル
  pwm_clear_irq(2);
  imu_mag_data_read();
  Ax=   -acceleration_mg[0]*0.001*GRAV;
  Ay=   -acceleration_mg[1]*0.001*GRAV;
  Az=    acceleration_mg[2]*0.001*GRAV;
  Wp=    angular_rate_mdps[0]*0.001*0.017453292;
  Wq=    angular_rate_mdps[1]*0.001*0.017453292;
  Wr=   -angular_rate_mdps[2]*0.001*0.017453292;
  Dmx=  -(magnetic_field_mgauss[0]);
  Dmy=   (magnetic_field_mgauss[1]);
  Dmz=  -(magnetic_field_mgauss[2]);
  /*
     回転行列
     [[-0.78435472 -0.62015392 -0.01402787]
     [ 0.61753358 -0.78277935  0.07686857]
     [-0.05865107  0.05162955  0.99694255]]
     中心座標
     -109.32529343620176 72.76584808916506 759.2285249891385
     W
     0.5498054412471614
     拡大係数
     0.002034773458122364 0.002173892202021849 0.0021819494099235273
    */

  //回転行列
  const float rot[9]={-0.78435472, -0.62015392, -0.01402787,
    0.61753358, -0.78277935,  0.07686857,
    -0.05865107,  0.05162955,  0.99694255};
  //中心座標
  const float center[3]={-109.32529343620176, 72.76584808916506, 759.2285249891385};
  //拡大係数
  const float zoom[3]={0.002034773458122364, 0.002173892202021849, 0.0021819494099235273};


  //回転・平行移動・拡大
  mx1 = zoom[0]*( rot[0]*Dmx +rot[1]*Dmy +rot[2]*Dmz -center[0]);
  my1 = zoom[1]*( rot[3]*Dmx +rot[4]*Dmy +rot[5]*Dmz -center[1]);
  mz1 = zoom[2]*( rot[6]*Dmx +rot[7]*Dmy +rot[8]*Dmz -center[2]);
  //逆回転
  Mx = rot[0]*mx1 +rot[3]*my1 +rot[6]*mz1;
  My = rot[1]*mx1 +rot[4]*my1 +rot[7]*mz1;
  Mz = rot[2]*mx1 +rot[5]*my1 +rot[8]*mz1; 
  mag_norm=sqrt(Mx*Mx +My*My +Mz*Mz);
  Mx/=mag_norm;
  My/=mag_norm;
  Mz/=mag_norm;

  //printf("%9.4f %12.5f %12.5f %12.5f\n",Time_for_debug, mx, my, mz );
  //Time_for_debug = Time_for_debug + 0.0025; 

  Hzcount=Hzcount+1;
  if(Hzcount==4){
    sem_release(&sem);
    Hzcount=0;
  }

  //角速度PID制御
  
  //最大角速度 e,a 6π,r 2π
  //最大角度30°

  //エレベータピッチレートq
  q_rate = Wq - Wqa;
  err_q = (Ref_q - q_rate );
  Sk_q = Sk_q + err_q;       //修正しました
  if (Sk_q > 30000.0){     //修正しました
    Sk_q = 30000.0;        //修正しました
  }                      //修正しました
  else if(Sk_q <-30000.0 ){//修正しました
    Sk_q =-30000.0;        //修正しました
  }                      //修正しました
  dk_q  = (err_q - Olderr3_q)*400;
  Dk_q  = Dk_q * Tc_rate/(Tc_rate + 0.0025) + dk_q*0.0025/(Tc_rate + 0.0025);
  Uq = Kp_q * (err_q + Sk_q*0.0025/Ti_q + Td_q * Dk_q);
  Olderr3_q = Olderr2_q;
  Olderr2_q = Olderr1_q;
  Olderr1_q = err_q;

  //エルロンロールレートp
  p_rate = Wp - Wpa;
  err_p = (Ref_p - p_rate );
  Sk_p = Sk_p + err_p;    //修正しました
  if (Sk_p > 30000.0)     //修正しました
  {                     //修正しました
    Sk_p = 30000.0;       //修正しました
  }                     //修正しました
  else if(Sk_p <-30000.0){//修正しました
    Sk_p =-30000.0;       //修正しました
  }                     //修正しました
  dk_p = (err_p - Olderr3_p)*400;
  Dk_p = Dk_p * Tc_rate/(Tc_rate + 0.0025) + dk_p * 0.0025/(Tc_rate + 0.0025);  
  Up = Kp_p * (err_p + Sk_p * 0.0025/Ti_p + Td_p * Dk_p);
  Olderr3_p = Olderr2_p;
  Olderr2_p = Olderr1_p;
  Olderr1_p = err_p;

  //ラダーr
  r_rate = Wr - Wra;
  err_r = (Ref_r - r_rate);
  Sk_r = Sk_r + err_r;   //修正しました
  if (Sk_r > 30000.0)    //修正しました
  {                    //修正しました
    Sk_r = 30000.0;      //修正しました
  }                    //修正しました
  else if(Sk_r <-30000.0)//修正しました
  {                    //修正しました
    Sk_r =-30000.0;      //修正しました
  }                    //修正しました
  dk_r =(err_r - Olderr3_r)*400;
  Dk_r =Dk_r * Tc_rate/(Tc_rate + 0.0025) +dk_r* 0.0025/(Tc_rate + 0.0025);  
  Ur=Kp_r * (err_r + Sk_r*0.0025/Ti_r * Sk_r + Td_r * Dk_r);
  Olderr3_r = Olderr2_r;
  Olderr2_r = Olderr1_r;
  Olderr1_r = err_r;

  if(Kalman_time>12.0){
    if(Data3<0.1){
      if(Data2>0.9 && Data4>0.9 && Safetycount==0){
        Safetycount=1;
        //      printf("A %d %f %f \n",Safetycount,Data2,Data4);
      }
      else if(Data2<-0.9 && Data4<-0.9 && Safetycount==1){
        Safetycount=0;
        //    printf("B %d %f %f \n",Safetycount,Data2,Data4);
      }
    }
    if(Safetycount==1){
      Com_fr=0.0;
      Com_fl=0.0;
      Com_rr=0.0;         
      Com_rl=0.0;
    }
    else if(Safetycount==0){
      Com_fr = Data3 + ( Uq -Up +Ur) * 0.25;
      Com_fl = Data3 + ( Uq +Up -Ur) * 0.25;
      Com_rr = Data3 + (-Uq -Up -Ur) * 0.25;         
      Com_rl = Data3 + (-Uq +Up +Ur) * 0.25;
    }
  }


#if 1    
  //Duty_rr=(float)(DUTYMAX-DUTYMIN)*Com_rr+DUTYMIN;
  //Duty_fr=(float)(DUTYMAX-DUTYMIN)*Com_fr+DUTYMIN;
  //Duty_rl=(float)(DUTYMAX-DUTYMIN)*Com_rl+DUTYMIN;
  //Duty_fl=(float)(DUTYMAX-DUTYMIN)*Com_fl+DUTYMIN;

  //if (Duty_rr>DUTYMAX-50.0)Duty_rr=DUTYMAX-50.0;
  //if (Duty_rr<DUTYMIN+15.0)Duty_rr=DUTYMIN+15.0;
  //if (Duty_fr>DUTYMAX-50.0)Duty_fr=DUTYMAX-50.0;
  //if (Duty_fr<DUTYMIN+15.0)Duty_fr=DUTYMIN+15.0;
  //if (Duty_rl>DUTYMAX-50.0)Duty_rl=DUTYMAX-50.0;
  //if (Duty_rl<DUTYMIN+15.0)Duty_rl=DUTYMIN+15.0;
  //if (Duty_fl>DUTYMAX-50.0)Duty_fl=DUTYMAX-50.0;
  //if (Duty_fl<DUTYMIN+15.0)Duty_fl=DUTYMIN+15.0;

  if(Data3<0.05 || Safetycount==1)
  {
    set_duty_rr(0.0);
    set_duty_fr(0.0);
    set_duty_rl(0.0);
    set_duty_fl(0.0);
    //pwm_set_chan_level(slice_num[0], PWM_CHAN_A, DUTYMIN);
    //pwm_set_chan_level(slice_num[0], PWM_CHAN_B, DUTYMIN);
    //pwm_set_chan_level(slice_num[1], PWM_CHAN_A, DUTYMIN);
    //pwm_set_chan_level(slice_num[1], PWM_CHAN_B, DUTYMIN);
    Sk_p=0.0;
    Sk_q=0.0;
    Sk_r=0.0;
    Olderr1_p=0.0;
    Olderr1_q=0.0;
    Olderr1_r=0.0;
    Olderr2_p=0.0;
    Olderr2_q=0.0;
    Olderr2_r=0.0;
    Olderr3_p=0.0;
    Olderr3_q=0.0;
    Olderr3_r=0.0;
    Dk_p=0.0;
    Dk_q=0.0;
    Dk_r=0.0;     //追加しました
    Up = 0.0;
    Uq = 0.0;
    Ur = 0.0;
    //oldDk_phi=0.0; //角度制御に持っていきました
    //oldDk_t=0.0;   //角度制御に持っていきました
  }
  else
  {
    set_duty_rr(Com_rr);
    set_duty_fr(Com_fr);
    set_duty_rl(Com_rl);
    set_duty_fl(Com_fl);
    //pwm_set_chan_level(slice_num[0], PWM_CHAN_A, Duty_rr);
    //pwm_set_chan_level(slice_num[0], PWM_CHAN_B, Duty_fr);
    //pwm_set_chan_level(slice_num[1], PWM_CHAN_A, Duty_rl);
    //pwm_set_chan_level(slice_num[1], PWM_CHAN_B, Duty_fl);
  }
#endif
}
/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/

void init_kalman(float mx, float my, float mz)
{
  Xe << 1.00, 0.0, 0.0, 0.0, mx, my, mz;
  Xp =Xe;
/*
  Q <<  0.00872, 0.0   , 0.0     , 0.0   , 0.0   , 0.0,
        0.0   , 0.00198, 0.0     , 0.0   , 0.0   , 0.0,
        0.0   , 0.0   , 0.000236 , 0.0   , 0.0   , 0.0,
        0.0   , 0.0   , 0.0      , 0.5e-5, 0.0   , 0.0,
        0.0   , 0.0   , 0.0      , 0.0   , 0.5e-5, 0.0,
        0.0   , 0.0   , 0.0      , 0.0   , 0.0   , 0.1e0;

  R <<  1e-1   , 0.0    , 0.0    , 0.0     , 0.0     , 0.0,
        0.0    , 1e-1   , 0.0    , 0.0     , 0.0     , 0.0,
        0.0    , 0.0    , 5e-1   , 0.0     , 0.0     , 0.0,
        0.0    , 0.0    , 0.0    , 1e-2    , 0.0     , 0.0,
        0.0    , 0.0    , 0.0    , 0.0     , 1e-2    , 0.0,
        0.0    , 0.0    , 0.0    , 0.0     , 0.0     , 1e-3;
*/
  Q <<  1.0e-5, 0.0   , 0.0     , 0.0   , 0.0   , 0.0,
        0.0   , 1.0e-5, 0.0     , 0.0   , 0.0   , 0.0,
        0.0   , 0.0   , 1.0e-5 , 0.0   , 0.0   , 0.0,
        0.0   , 0.0   , 0.0      , 0.5e-3, 0.0   , 0.0,
        0.0   , 0.0   , 0.0      , 0.0   , 0.2e-3, 0.0,
        0.0   , 0.0   , 0.0      , 0.0   , 0.0   , 0.1e-2;

  R <<  1e-1   , 0.0    , 0.0    , 0.0     , 0.0     , 0.0,
        0.0    , 1e-1   , 0.0    , 0.0     , 0.0     , 0.0,
        0.0    , 0.0    , 5e-1   , 0.0     , 0.0     , 0.0,
        0.0    , 0.0    , 0.0    , 1e-4    , 0.0     , 0.0,
        0.0    , 0.0    , 0.0    , 0.0     , 1e-4    , 0.0,
        0.0    , 0.0    , 0.0    , 0.0     , 0.0     , 1e-2;


  G << -1.0,-1.0,-1.0, 0.0, 0.0, 0.0, 
        1.0,-1.0, 1.0, 0.0, 0.0, 0.0, 
        1.0, 1.0,-1.0, 0.0, 0.0, 0.0, 
       -1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

  Beta << 0.0, 0.0, 0.0;

  P <<  1e0,   0,   0,   0,     0,     0,     0,  
          0, 1e0,   0,   0,     0,     0,     0,
          0,   0, 1e0,   0,     0,     0,     0,
          0,   0,   0, 1e0,     0,     0,     0, 
          0,   0,   0,   0, 1.0e0,     0,     0,  
          0,   0,   0,   0,     0, 1.0e0,     0,  
          0,   0,   0,   0,     0,     0, 1.0e0;
}


int main(void)
{
  uint16_t f;
  uint8_t waittime=5;
  float old_kalman_time;
  //float Thetaplas=0.0,Psiplas=0.0,Phiplas=0.0;
  const uint LED_PIN = 25;        //LED_PIN=0  
  gpio_init(LED_PIN);             //gpioを使えるようにする
  gpio_set_dir(25, GPIO_OUT);
  stdio_init_all();
  imu_mag_init();
  serial_settei();
  gpio_put(LED_PIN, 1);
  sleep_ms(1000);
  gpio_put(LED_PIN, 0);
#if 0
  for (uint8_t i=0;i<waittime;i++)
  {
     printf("#Please wait %d[s] ! \n",waittime-i);
     sleep_ms(1000);
  }
  printf("#Start Kalman Filter\n");
#endif
  
  sem_init(&sem, 0, 1);
  multicore_launch_core1(kalman);
  pwm_settei();
  printf("#Start Kalman filter test !\n");

  Q0av = 0.0;
  while(Q0av<0.995 || isinf(Q0av))
  {
    MN = 0.0;
    ME = 0.0;
    MD = 0.0;
    Wpa = 0.0;
    Wqa = 0.0;
    Wra = 0.0;
    f = 0;
    Kalman_time = 0.0;
    old_kalman_time = Kalman_time;
    while(Kalman_time<2.0)
    {
      if(Kalman_time!=old_kalman_time)
      {
        MN=MN+Mx;
        ME=ME+My;
        MD=MD+Mz;
        Wqa=Wq+Wqa;
        Wpa=Wp+Wpa;
        Wra=Wr+Wra;
        f=f+1;
        old_kalman_time = Kalman_time;
      }
    }
    Wpa=Wpa/f;
    Wqa=Wqa/f;
    Wra=Wra/f;
    MN=MN/f;
    ME=ME/f;
    MD=MD/f;
    
    Kalman_flag = 0;//Kalman filter 一時停止
    init_kalman(Wpa, Wqa, Wra);//カルマンフィルタ初期化
    //Xe(4,0)=Wpa;
    //Xe(5,0)=Wqa;
    //Xe(6,0)=Wra;
    //Xp(4,0)=Wpa;
    //Xp(5,0)=Wqa;
    //Xp(6,0)=Wra;

    printf("#omega ave %f %f %f\n",Wpa, Wqa, Wra);
    printf("#omega now %f %f %f\n",Wp, Wq, Wr);
    printf("#mag ave %f %f %f\n",MN, ME, MD);
    printf("#mag now %f %f %f\n",Mx, My, Mz);
    printf("#Kalman init Xe %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
        Xe(0,0),Xe(1,0),Xe(2,0),Xe(3,0),Xe(4,0),Xe(5,0),Xe(6,0));
    printf("#Kalman init Xp %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
        Xp(0,0),Xp(1,0),Xp(2,0),Xp(3,0),Xp(4,0),Xp(5,0),Xp(6,0));
    printf("#Kalman init  P %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n",
        P(0,0),P(1,1),P(2,2),P(3,3),P(4,4),P(5,5),P(6,6));

    Kalman_time=0.0;
    Kalman_flag = 1;//Kalman filter 再開


    f=0;
    old_kalman_time=Kalman_time;
    while(Kalman_time<6.0){
      if(Kalman_time!=old_kalman_time)
      {
        printf("#Kalman time %5.2f\r", Kalman_time);
        Q0av += Xe(0,0);
        f=f+1;
        old_kalman_time=Kalman_time;
      }

    }
    printf("\n");
    Q0av/=f;
    Q0av=abs(Q0av);
    printf("#Q0av=%f\n",Q0av);
    if(Q0av<0.995 || isinf(Q0av)){
      printf("#Kalman filter Failure !\n");
      for (int i=0;i<100;i++)
      {
        gpio_put(LED_PIN, 1);
        sleep_ms(10);
        gpio_put(LED_PIN, 0);
        sleep_ms(10);
      }
    }
    else{
      printf("#Kalman filter Sucess !\n");
      for (int i=0;i<10;i++)
      {
        gpio_put(LED_PIN, 1);
        sleep_ms(100);
        gpio_put(LED_PIN, 0);
        sleep_ms(100);
      }
    }
  }

  f=0;
  old_kalman_time=Kalman_time;
  while(Kalman_time<10.0){
    if(Kalman_time!=old_kalman_time)
    {
      printf("#Kalman time %5.2f\r", Kalman_time);
      Phiav += Phi;
      Thetaav += Theta;
      Psiav += Psi;
      f=f+1;
      old_kalman_time=Kalman_time;
    }
  }
  Phiav/=f;
  Thetaav/=f;
  Psiav/=f;
  printf("\n");
  printf("#angle ave %f %f %f\n",Phiav, Thetaav, Psiav);

  while(1)
  { 
  }

}
void read_sensor(void)
{
  float mx1,my1,mz1;

  imu_mag_data_read();
  Wp=    angular_rate_mdps[0]*0.001*0.017453292;
  Wq=    angular_rate_mdps[1]*0.001*0.017453292;
  Wr=   -angular_rate_mdps[2]*0.001*0.017453292;
  Dmx=  -(magnetic_field_mgauss[0]);
  Dmy=   (magnetic_field_mgauss[1]);
  Dmz=  -(magnetic_field_mgauss[2]);

  //回転行列
  const float rot[9]={-0.78435472, -0.62015392, -0.01402787,
    0.61753358, -0.78277935,  0.07686857,
    -0.05865107,  0.05162955,  0.99694255};
  //中心座標
  const float center[3]={-109.32529343620176, 72.76584808916506, 759.2285249891385};
  //拡大係数
  const float zoom[3]={0.002034773458122364, 0.002173892202021849, 0.0021819494099235273};

  //回転・平行移動・拡大
  mx1 = zoom[0]*( rot[0]*Dmx +rot[1]*Dmy +rot[2]*Dmz -center[0]);
  my1 = zoom[1]*( rot[3]*Dmx +rot[4]*Dmy +rot[5]*Dmz -center[1]);
  mz1 = zoom[2]*( rot[6]*Dmx +rot[7]*Dmy +rot[8]*Dmz -center[2]);
  //逆回転
  Mx = rot[0]*mx1 +rot[3]*my1 +rot[6]*mz1;
  My = rot[1]*mx1 +rot[4]*my1 +rot[7]*mz1;
  Mz = rot[2]*mx1 +rot[5]*my1 +rot[8]*mz1; 
  float mag_norm=sqrt(Mx*Mx +My*My +Mz*Mz);
  Mx/=mag_norm;
  My/=mag_norm;
  Mz/=mag_norm;
}


float rocking_wings(void)
{
  float freq = 1.5;//Hz
  float amp  = 25;//deg
  
  if(Rocking_flag == 0)
  {
    Rocking_time = Kalman_time;
    Rocking_flag = 1;
  }
  else if(Kalman_time - Rocking_time < 3.0)
  {
    return amp*M_PI/180.0*sin(freq*2*M_PI*Kalman_time);
  }
  else
  {
    return 0.0;
  }
  return 0.0;
}
