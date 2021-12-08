#include "tobu.hpp"
/* -----------------------------------------------------------------------------------------------------------------------*/
//マルチコア設定
semaphore_t sem;

Matrix<float, 7 ,1> xp = MatrixXf::Zero(7,1);
Matrix<float, 7 ,1> xe = MatrixXf::Zero(7,1);
Matrix<float, 7 ,1> x_sim = MatrixXf::Zero(7,1);
Matrix<float, 7 ,7> P = MatrixXf::Identity(7,7);
Matrix<float, 6 ,1> z = MatrixXf::Zero(6,1);
Matrix<float, 3, 1> omega_m = MatrixXf::Zero(3, 1);
Matrix<float, 6, 6> Q = MatrixXf::Identity(6, 6)*1;
Matrix<float, 6, 6> R = MatrixXf::Identity(6, 6)*1;
Matrix<float, 7 ,6> G;
Matrix<float, 3 ,1> beta;
float ax,ay,az,wp,wq,wr,mx,my,mz;
volatile float dmx,dmy,dmz;
float Wqa=0.0, Wpa=0.0,Wra=0.0; 
float phi,theta,psi;
float PHI,THETA,PSI;
float psiav=0.0, phiav=0.0, thetaav=0.0;
volatile float kalman_time=0.0;
float spsi=0.0,sphi=0.0,st=0.0;
float se=0.0,sa=0.0,sr=0.0;
float dt=0.01;
short Hzcount=1;
short safetycount=1;
float Duty_rr,Duty_rl,Duty_fl,Duty_fr; 
float ref_e,ref_r,ref_a;
float err_e=0.0,err_a=0.0,err_r=0.0;
float sk_e=0.0, sk_a=0.0, sk_r=0.0;
float sk_phi=0.0, sk_t=0.0, sk_psi=0.0;
float Dk_a=0.0, Dk_e=0.0, Dk_r=0.0;
float Dk_phi=0.0, Dk_t=0.0;
float olderr_e=0.0, olderr_a=0.0, olderr_r=0.0;
float olderr2_a=0.0,olderr3_a=0.0,olderr2_r=0.0,olderr3_r=0.0,olderr2_e=0.0,olderr3_e=0.0;
float olderr_phi=0.0, olderr2_phi=0.0, olderr3_phi=0.0;
float olderr_t=0.0, olderr2_t=0.0, olderr3_t=0.0;
float olderr_psi=0.0;
char sbuf[256];
const uint LED_PIN = 25;
uint32_t logcount=0;
uint32_t printcount=0;
const uint16_t Logdatanum=48000;
float logdata[Logdatanum]={0.0};
const uint8_t DATANUM=28;
float Tlog = 0.0;
float Led_timer = 0.0;
uint8_t Led_flag = 1;
float Time_for_debug = 0.0;
//Gain
const float Ta=0.018;
const float Tb=0.012;
const float Ti_roll=100;
const float Ti_pitch=100;
const float Ti=100;
const float Td=0.0000014;
const float Td_roll=0.00004;
const float Td_pitch=0.00004;
const float kpe=0.1,kpa=0.1,kpr=0.1;
const float kppsi=1.0,kpphi=1.0,kpt=1.0;

/*--------------------------------------------------------------------------------------------------------------------------------------*/

void kalman(void){

  float ref_t,ref_phi,ref_psi;
  float err_psi, err_phi, err_t;
  float dk_a,dk_e,dk_r;
  float dk_phi,dk_psi,dk_t;
  float data2MID,data4MID;

  /*------------------------------------------------------------------------------------------------------------------------------------*/
  while(1)
  {
    sem_acquire_blocking(&sem);
    sem_reset(&sem, 0);
    //printf("%f\n",time);
    ref_t=(Data2)*0.523598775598299*2;
    ref_phi=(Data4)*0.52398775598299*2;
    ref_psi=Data1;

    dt=0.01;
    omega_m << wp, wq, wr;
    z << ax, ay, az, mx, my, mz;//ここに入れる
    //--Begin Extended Kalman Filter--
    ekf(xp, xe, P, z, omega_m, Q, R, G*dt, beta, dt);
    kalman_time = kalman_time + 0.01;

    phi=Phl(xe);
    theta=Theta(xe);
    psi=Psi(xe);

    PSI=psi-psiav;
    PHI=phi-phiav;
    THETA=theta-thetaav;


    if(Data3>=0.25){

      //phi
      olderr3_phi = olderr2_phi;
      olderr2_phi = olderr_phi;
      olderr_phi = err_phi;
      err_phi = (ref_phi - PHI);
      sk_phi = sk_phi + err_phi;//修正しました
      if (sk_phi>30000)         //修正しました
      {                         //修正しました
        sk_phi = 30000;         //修正しました
      }                         //修正しました
      else if(sk_phi<-30000)    //修正しました
      {                         //修正しました
        sk_phi=-30000;          //修正しました
      }                         //修正しました
      dk_phi = (err_phi-olderr3_phi)*100;
      Dk_phi = Dk_phi*Ta/(Ta+0.01) + 0.01/(Ta+0.01)*dk_phi;
      sphi = kpphi*(err_phi + 1/Ti_roll*sk_phi*0.01 + Td_roll*Dk_phi);
      //theta
      olderr3_t = olderr2_t;
      olderr2_t = olderr_t;
      olderr_t = err_t;
      err_t = (ref_t - THETA);
      sk_t = sk_t + err_t;//修正しました
      if (sk_t>30000)     //修正しました
      {                   //修正しました
        sk_t = 30000;     //修正しました
      }                   //修正しました
      else if(sk_t<-30000)//修正しました
      {                   //修正しました
        sk_t =-30000;     //修正しました
      }                   //修正しました
      dk_t = (err_t - olderr3_t)*100;
      Dk_t = Dk_t*Ta/(Ta+0.01) + 0.01/(Ta+0.01)*dk_t;
      st = kpt*(err_t + 1/Ti_pitch*sk_t*0.01 + Td_pitch*Dk_t);
    }
    else{
      sk_phi=0;
      sk_t=0;
      sphi=0;
      st=0;
      olderr_phi=0;
      olderr2_phi=0;
      olderr3_phi=0;
      olderr_t=0;
      olderr2_t=0;
      olderr3_t=0;
      phiav=phi;
      thetaav=theta;
      data2MID=Data2;
      data4MID=Data4;
      Dk_phi = 0.0;//追加しました
      Dk_t = 0.0;  //追加しました
    }


    //printf("kalman %04f\n", Data6);    


    //スティック上
    if(Data6<-0.2)
    {
      //printf("upper log sw\n");
      if(logcount<47500-DATANUM){
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
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
        logdata[logcount++]=0;
#if 0
        logdata[logcount++]=xe(0,0);
        logdata[logcount++]=xe(1,0);
        logdata[logcount++]=xe(2,0);
        logdata[logcount++]=xe(3,0);
        logdata[logcount++]=xe(4,0);
        logdata[logcount++]=xe(5,0);
        logdata[logcount++]=xe(6,0);
        logdata[logcount++]=wp;
        logdata[logcount++]=wq;
        logdata[logcount++]=wr;
        logdata[logcount++]=ax;
        logdata[logcount++]=ay;
        logdata[logcount++]=az;
        logdata[logcount++]=mx;
        logdata[logcount++]=my;
        logdata[logcount++]=mz;
        logdata[logcount++]=ref_a;
        logdata[logcount++]=ref_e;
        logdata[logcount++]=ref_r;
        logdata[logcount++]=PHI;
        logdata[logcount++]=THETA;
        logdata[logcount++]=PSI;
        logdata[logcount++]=ref_phi;
        logdata[logcount++]=ref_t;
        logdata[logcount++]=ref_psi;
        logdata[logcount++]=sa;
        logdata[logcount++]=se;
        logdata[logcount++]=sr;
#endif
      }
      else{
        gpio_put(LED_PIN,0);
      }
    }

    //スティック下
    else if(Data6>0.2){
      //printf("lower log sw\n");
      if(Data3 < 0.3){
        if(printcount+DATANUM < 47500 ){
          for (uint8_t i=0;i<DATANUM;i++){
            if(i==0){
              printf("%8.2f ", Tlog);
              Tlog = Tlog + 0.01;
            }
            sprintf(sbuf,"%12.5f",logdata[printcount+i]);
            printf("%s",sbuf);
          }

          printf("\n");
          printcount=printcount+DATANUM;
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
      if(kalman_time>22.0)gpio_put(LED_PIN, 0);
      printcount = 0;
      logcount = 0;
      Tlog = 0.0;
    }
  }
}


/*--------------------------------------------------------------------------------------------------------------------------------------*/
void MAINLOOP(void)
{
  float mx1,my1,mz1;
  float dk_e, dk_a, dk_r;
  float duty_rr,duty_rl,duty_fl,duty_fr; 
  
  //e エレベータ a エルロン r ラダー t スロットル
  pwm_clear_irq(2);
  imu_mag_data_read();
  ax=   -acceleration_mg[0]*0.001*GRAV;
  ay=   -acceleration_mg[1]*0.001*GRAV;
  az=    acceleration_mg[2]*0.001*GRAV;
  wp=    angular_rate_mdps[0]*0.001*0.017453292;
  wq=    angular_rate_mdps[1]*0.001*0.017453292;
  wr=   -angular_rate_mdps[2]*0.001*0.017453292;
  dmx=  -(magnetic_field_mgauss[0]);
  dmy=   (magnetic_field_mgauss[1]);
  dmz=  -(magnetic_field_mgauss[2]);
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
  mx1 = zoom[0]*( rot[0]*dmx +rot[1]*dmy +rot[2]*dmz -center[0]);
  my1 = zoom[1]*( rot[3]*dmx +rot[4]*dmy +rot[5]*dmz -center[1]);
  mz1 = zoom[2]*( rot[6]*dmx +rot[7]*dmy +rot[8]*dmz -center[2]);
  //逆回転
  mx = rot[0]*mx1 +rot[3]*my1 +rot[6]*mz1;
  my = rot[1]*mx1 +rot[4]*my1 +rot[7]*mz1;
  mz = rot[2]*mx1 +rot[5]*my1 +rot[8]*mz1; 
  float mag_norm=sqrt(mx*mx +my*my +mz*mz);
  mx/=mag_norm;
  my/=mag_norm;
  mz/=mag_norm;

  //printf("%9.4f %12.5f %12.5f %12.5f\n",Time_for_debug, mx, my, mz );
  //Time_for_debug = Time_for_debug + 0.0025; 

  if(Hzcount<=3){
    Hzcount=Hzcount+1;
  }
  else{
    Hzcount=1;
    sem_release(&sem);
  }

#if 1    
  //姿勢安定化
  //最大角速度 e,a 6π,r 2π
  //最大角度30°
  ref_e=st;
  ref_a=sphi;
  ref_r=Data1*6.283185307;


  //エレベータq
  olderr3_e=olderr2_e;
  olderr2_e=olderr_e;
  olderr_e=err_e;
  err_e=(ref_e - (wq-Wqa));
  sk_e=sk_e+err_e;       //修正しました
  if (sk_e > 30000){     //修正しました
    sk_e = 30000;        //修正しました
  }                      //修正しました
  else if(sk_e <-30000 ){//修正しました
    sk_e =-30000;        //修正しました
  }                      //修正しました
  dk_e=(err_e-olderr3_e)*400;
  Dk_e=Dk_e*Tb/(Tb+0.0025)+0.0025/(Tb+0.0025)*dk_e;
  se=kpe*(err_e+1/Ti*sk_e*0.0025+Td*Dk_e);
     

  //エルロンp
  olderr3_a=olderr2_a;
  olderr2_a=olderr_a;
  olderr_a=err_a;
  err_a=(ref_a - (wp-Wpa));
  sk_a = sk_a+err_a;    //修正しました
  if (sk_a > 30000)     //修正しました
  {                     //修正しました
    sk_a = 30000;       //修正しました
  }                     //修正しました
  else if(sk_a <-30000){//修正しました
    sk_a =-30000;       //修正しました
  }                     //修正しました
  dk_a=(err_a-olderr3_a)*400;
  Dk_a=Dk_a*Tb/(Tb+0.0025)+0.0025/(Tb+0.0025)*dk_a;  
  sa=kpa*(err_a+1/Ti*sk_a*0.0025+Td*Dk_a);
  

  //ラダーr
  olderr3_r=olderr2_r;
  olderr2_r=olderr_r;
  olderr_r=err_r;
  err_r=(ref_r - (wr-Wra));
  sk_r = sk_r+err_r;   //修正しました
  if (sk_r > 30000)    //修正しました
  {                    //修正しました
    sk_r = 30000;      //修正しました
  }                    //修正しました
  else if(sk_r <-30000)//修正しました
  {                    //修正しました
    sk_r =-30000;      //修正しました
  }                    //修正しました
  dk_r=(err_r-olderr3_r)*400;
  Dk_r=Dk_r*Tb/(Tb+0.0025)+0.0025/(Tb+0.0025)*dk_r;  
  sr=kpr*(err_r+1/Ti*sk_r*0.0025+Td*Dk_r);

  if(kalman_time>=25.0){
    if(Data3<0.1){
#if 1
      if(Data2>0.9 && Data4>0.9 && safetycount==0){
        safetycount=1;
        //      printf("A %d %f %f \n",safetycount,Data2,Data4);
      }
      else if(Data2<-0.9 && Data4<-0.9 && safetycount==1){
        safetycount=0;
        //    printf("B %d %f %f \n",safetycount,Data2,Data4);
      }
    }


    if(safetycount==1){
      Duty_fr=0.0;
      Duty_fl=0.0;
      Duty_rr=0.0;         
      Duty_rl=0.0;

    }
    else if(safetycount==0){
      Duty_fr = Data3 + ( se -sa +sr) * 0.25;
      Duty_fl = Data3 + ( se +sa -sr) * 0.25;
      Duty_rr = Data3 + (-se -sa -sr) * 0.25;         
      Duty_rl = Data3 + (-se +sa +sr) * 0.25;
    }
  }
#endif
#endif
#if 0   
  Duty_fr=Data3;
  Duty_fl=Data3;
  Duty_rr=Data3;         
  Duty_rl=Data3;
  else{
    Duty_fr=Data3*0;
    Duty_fl=Data3*0;
    Duty_rr=Data3*0;         
    Duty_rl=Data3*0;
  }
#endif
  tight_loop_contents();

  duty_rr=(float)(DUTYMAX-DUTYMIN)*Duty_rr+DUTYMIN;
  duty_fr=(float)(DUTYMAX-DUTYMIN)*Duty_fr+DUTYMIN;
  duty_rl=(float)(DUTYMAX-DUTYMIN)*Duty_rl+DUTYMIN;
  duty_fl=(float)(DUTYMAX-DUTYMIN)*Duty_fl+DUTYMIN;

  if (duty_rr>DUTYMAX-50.0)duty_rr=DUTYMAX-50.0;
  if (duty_rr<DUTYMIN+15.0)duty_rr=DUTYMIN+15.0;
  if (duty_fr>DUTYMAX-50.0)duty_fr=DUTYMAX-50.0;
  if (duty_fr<DUTYMIN+15.0)duty_fr=DUTYMIN+15.0;
  if (duty_rl>DUTYMAX-50.0)duty_rl=DUTYMAX-50.0;
  if (duty_rl<DUTYMIN+15.0)duty_rl=DUTYMIN+15.0;
  if (duty_fl>DUTYMAX-50.0)duty_fl=DUTYMAX-50.0;
  if (duty_fl<DUTYMIN+15.0)duty_fl=DUTYMIN+15.0;

  if(Data3<0.05 || safetycount==1)
  {
    pwm_set_chan_level(slice_num[0], PWM_CHAN_A, DUTYMIN);
    pwm_set_chan_level(slice_num[0], PWM_CHAN_B, DUTYMIN);
    pwm_set_chan_level(slice_num[1], PWM_CHAN_A, DUTYMIN);
    pwm_set_chan_level(slice_num[1], PWM_CHAN_B, DUTYMIN);
    sk_r=0.0;
    sk_a=0.0;
    sk_e=0.0;
    err_r=0.0;
    err_a=0.0;
    err_e=0.0;
    Dk_r=0.0;
    Dk_a=0.0;
    Dk_e=0.0;     //追加しました
    //oldDk_phi=0.0; //角度制御に持っていきました
    //oldDk_t=0.0;   //角度制御に持っていきました


  }
  else
  {
    pwm_set_chan_level(slice_num[0], PWM_CHAN_A, duty_rr);
    pwm_set_chan_level(slice_num[0], PWM_CHAN_B, duty_fr);
    pwm_set_chan_level(slice_num[1], PWM_CHAN_A, duty_rl);
    pwm_set_chan_level(slice_num[1], PWM_CHAN_B, duty_fl);
  }
  //printf("%04f %04f %04f %04f %04f %04f \n",Olddata[0],Olddata[1],Olddata[2],Olddata[3],Olddata[4],Olddata[5]);    
  //e_time=time_us_64();
  //d_time=e_time-s_time;
  //printf("右前%04f   左前%04f   右後ろ%04f   左後ろ%04f\n",duty_fr,duty_fl,duty_rr,duty_rl);
  //printf("%04f \n",t);
  //t=t+1;

}
/*-----------------------------------------------------------------------------------------------------------------------------------------------------------*/
int main(void)
{
  float w,f=1;
  uint8_t waittime=5;
  float Thetaplas=0.0,Psiplas=0.0,Phiplas=0.0;
  const uint LED_PIN = 25;        //LED_PIN=0  
  gpio_init(LED_PIN);             //gpioを使えるようにする
  gpio_set_dir(25, GPIO_OUT);

  xe << 1.00, 0.0, 0.0, 0.0,0.0,0.0, 0.0;
  xp =xe;

  Q <<  0.0363, 0.0   , 0.0     , 0.0   , 0.0   , 0.0,
        0.0   , 0.0298, 0.0     , 0.0   , 0.0   , 0.0,
        0.0   , 0.0   , 0.000513, 0.0   , 0.0   , 0.0,
        0.0   , 0.0   , 0.0     , 0.5e-5, 0.0   , 0.0,
        0.0   , 0.0   , 0.0     , 0.0   , 0.5e-5, 0.0,
        0.0   , 0.0   , 0.0     , 0.0   , 0.0   , 1.0;

  R <<  29.1   , 0.0    , 0.0    , 0.0     , 0.0     , 0.0,
        0.0    , 29.9   , 0.0    , 0.0     , 0.0     , 0.0,
        0.0    , 0.0    , 32.2   , 0.0     , 0.0     , 0.0,
        0.0    , 0.0    , 0.0    , 0.000416, 0.0     , 0.0,
        0.0    , 0.0    , 0.0    , 0.0     , 0.000372, 0.0,
        0.0    , 0.0    , 0.0    , 0.0     , 0.0     , 0.000467;

  G <<   1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 
        -1.0, 1.0,-1.0, 0.0, 0.0, 0.0, 
        -1.0,-1.0, 1.0, 0.0, 0.0, 0.0, 
        1.0,-1.0,-1.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0;

  beta << 0.0, 0.0, 0.0;

  P <<  1e0,0,0,0,0,0,0,  
          0,1e0,0,0,0,0,0,
          0,  0,1e0,0,0,0,
          0,0,0,1e0,0,0,0, 
          0,0,0,0,1e0,0,0,  
          0,0,0,0,0,1e0,0,  
          0,0,0,0,0,0,1e0;


  stdio_init_all();
  imu_mag_init();
  serial_settei();

  gpio_put(LED_PIN, 1);
  sleep_ms(1000);
  gpio_put(LED_PIN, 0);

  for (uint8_t i=0;i<waittime;i++)
  {
     printf("#Please wait %d[s] ! \n",waittime-i);
     sleep_ms(1000);
  }
  printf("#Start Kalman Filter\n");
  

  while(f<=400){
    float mx1,my1,mz1;
    imu_mag_data_read();
    wp=    angular_rate_mdps[0]*0.001*0.017453292;
    wq=    angular_rate_mdps[1]*0.001*0.017453292;
    wr=   -angular_rate_mdps[2]*0.001*0.017453292;
    dmx=  -(magnetic_field_mgauss[0]);
    dmy=   (magnetic_field_mgauss[1]);
    dmz=  -(magnetic_field_mgauss[2]);

    //回転行列
    const float rot[9]={-0.78435472, -0.62015392, -0.01402787,
      0.61753358, -0.78277935,  0.07686857,
      -0.05865107,  0.05162955,  0.99694255};
    //中心座標
    const float center[3]={-109.32529343620176, 72.76584808916506, 759.2285249891385};
    //拡大係数
    const float zoom[3]={0.002034773458122364, 0.002173892202021849, 0.0021819494099235273};

    //回転・平行移動・拡大
    mx1 = zoom[0]*( rot[0]*dmx +rot[1]*dmy +rot[2]*dmz -center[0]);
    my1 = zoom[1]*( rot[3]*dmx +rot[4]*dmy +rot[5]*dmz -center[1]);
    mz1 = zoom[2]*( rot[6]*dmx +rot[7]*dmy +rot[8]*dmz -center[2]);
    //逆回転
    mx = rot[0]*mx1 +rot[3]*my1 +rot[6]*mz1;
    my = rot[1]*mx1 +rot[4]*my1 +rot[7]*mz1;
    mz = rot[2]*mx1 +rot[5]*my1 +rot[8]*mz1; 
    float mag_norm=sqrt(mx*mx +my*my +mz*mz);
    mx/=mag_norm;
    my/=mag_norm;
    mz/=mag_norm;

    MN=MN+mx;
    ME=ME+my;
    MD=MD+mz;
    Wqa=wq+Wqa;
    Wpa=wp+Wpa;
    Wra=wr+Wra;
    f=f+1;
    sleep_us(2500);
  }
  Wqa=Wqa/400;
  Wpa=Wpa/400;
  Wra=Wra/400;
  MN=MN/400;
  MD=MN/400;
  ME=MN/400;
  
  sem_init(&sem, 0, 1);
  multicore_launch_core1(kalman);

  pwm_settei();

  while(kalman_time<20.0);

  f=0.0;
  while(f<4000){
    Phiplas=Phiplas+phi;
    Thetaplas=Thetaplas+theta;
    Psiplas=Psiplas+psi;
    f=f+1;
  }
  phiav=Phiplas/4000;
  thetaav=Thetaplas/4000;
  psiav=Psiplas/4000;
  
  for (int i=0;i<10;i++)
  {
    gpio_put(LED_PIN, 1);
    sleep_ms(100);
    gpio_put(LED_PIN, 0);
    sleep_ms(100);
  }

  while(1);
  

}

