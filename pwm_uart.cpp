#include "pwm_uart.hpp"

float shigma_e,shiguma_a,shiguma_r;
static int chars_rxed = 0;
static int data_num=0;
uint8_t sbus_data[25];
uint8_t ch =0;
uint slice_num[2];
uint16_t Olddata[6];
uint16_t Chdata[6];
float Data1=0.0,Data2=0.0,Data3=0.0,Data4=0.0,Data5=0.0,Data6=0.0;

//関数の宣言
//uint8_t serial_settei(void);
//uint8_t pwm_settei();
void on_uart_rx();
/*
 *   WARNING:
 *   Functions declare in this section are defined at the end of this file
 *   and are strictly related to the hardware platform used.
 *
 */

//シリアル設定
uint8_t serial_settei(void){

  // Set up our UART with a basic baud rate.
  uart_init(UART_ID, 2400);

  // Set the TX and RX pins by using the function select on the GPIO
  // Set dat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asheet for more information on function select
  gpio_set_function(UART_TX_PIN, GPIO_FUNC_UART);
  gpio_set_function(UART_RX_PIN, GPIO_FUNC_UART);

  // Actually, we want a different speed
  // The call will return the actual baud rate selected, which will be as close as
  // possible to that requested
  int actual = uart_set_baudrate(UART_ID, BAUD_RATE);

  // Set UART flow control CTS/RTS, we don't want these, so turn them off
  uart_set_hw_flow(UART_ID, false, false);

  // Set our data format
  uart_set_format(UART_ID, DATA_BITS, STOP_BITS, PARITY);

  // Turn off FIFO's - we want to do this character by character
  uart_set_fifo_enabled(UART_ID, true);

  // Set up a RX interrupt
  // We need to set up the handler first
  // Select correct interrupt for the UART we are using
  int UART_IRQ = UART_ID == uart0 ? UART0_IRQ : UART1_IRQ;

  // And set up and enable the interrupt handlers
  irq_set_exclusive_handler(UART_IRQ, on_uart_rx);
  irq_set_enabled(UART_IRQ, true);

  // Now enable the UART to send interrupts - RX only
  uart_set_irq_enables(UART_ID, true, false);
  return 0;
}

uint8_t pwm_settei(void){
  // PWMの設定
  //
  // GPIO 2 RR Motor 
  // GPOI 3 FR Motor
  // GPIO 4 RL Motor
  // GPIO 5 FL Motor
  // GPIO 6 Servo
  //
  // Set period T
  // T=(wrap+1)*clkdiv/sysclock
  // T = (3124+1)*100/125e6 = 3125e2/125e6=25e-4=0.0025s(=400Hz)
  //
  // Set Duty
  // Duty=clkdiv*PWM_CHAN_level/sysclock
  // Duty = 100 * 2500 /125e6 = 25.0e4/125e6=0.2e-2=0.002s=2ms
  // Duty = 100 * 1250 /125e6 = 12.5e4/125e6=0.1e-2=0.001s=1ms
  //
  // Tell GPIO 2-5 they are allocated to the PWM
  gpio_set_function(2, GPIO_FUNC_PWM);
  gpio_set_function(3, GPIO_FUNC_PWM);
  gpio_set_function(4, GPIO_FUNC_PWM);
  gpio_set_function(5, GPIO_FUNC_PWM);
  // Find out which PWM slice is connected to GPIO 3 and 4
  slice_num[0] = pwm_gpio_to_slice_num(3);
  slice_num[1] = pwm_gpio_to_slice_num(4);

  // Set period
  pwm_set_wrap(slice_num[0], 3124);
  pwm_set_clkdiv(slice_num[0], 100.0);
  pwm_set_wrap(slice_num[1], 3124);
  pwm_set_clkdiv(slice_num[1], 100.0);

//以下の#ifで1にするとESCキャリブレーションが動作する
#if 0
  //Start ESC calibration
  // Set the PWM Duty Maximum
  pwm_set_chan_level(slice_num[0], PWM_CHAN_A, DUTYMAX);
  pwm_set_chan_level(slice_num[0], PWM_CHAN_B, DUTYMAX);
  pwm_set_chan_level(slice_num[1], PWM_CHAN_A, DUTYMAX);
  pwm_set_chan_level(slice_num[1], PWM_CHAN_B, DUTYMAX);
  //Enable PWM Signal
  pwm_set_enabled(slice_num[0], true);
  pwm_set_enabled(slice_num[1], true);
  // Wait 4s
  sleep_ms(4000);
  // Set the PWM Duty minimum
  pwm_set_chan_level(slice_num[0], PWM_CHAN_A, DUTYMIN);
  pwm_set_chan_level(slice_num[0], PWM_CHAN_B, DUTYMIN);
  pwm_set_chan_level(slice_num[1], PWM_CHAN_A, DUTYMIN);
  pwm_set_chan_level(slice_num[1], PWM_CHAN_B, DUTYMIN);
#else
  //Start ESC control without ESC calibration 
  // Set the PWM Duty minimum
  pwm_set_chan_level(slice_num[0], PWM_CHAN_A, DUTYMIN);
  pwm_set_chan_level(slice_num[0], PWM_CHAN_B, DUTYMIN);
  pwm_set_chan_level(slice_num[1], PWM_CHAN_A, DUTYMIN);
  pwm_set_chan_level(slice_num[1], PWM_CHAN_B, DUTYMIN);
  //Enable PWM Signal
  pwm_set_enabled(slice_num[0], true);
  pwm_set_enabled(slice_num[1], true);
#endif

  pwm_clear_irq(slice_num[1]);
  pwm_set_irq_enabled(slice_num[1], true);
  irq_set_exclusive_handler(PWM_IRQ_WRAP,MAINLOOP);
  irq_set_enabled(PWM_IRQ_WRAP, true);
  return 0;
}

void set_duty_fr(float duty)
{
    duty=(float)(DUTYMAX-DUTYMIN)*duty+DUTYMIN;
    if (duty>DUTYMAX-50)duty=DUTYMAX-50;
    if (duty<DUTYMIN+15)duty=DUTYMIN+15;
    pwm_set_chan_level(slice_num[0], PWM_CHAN_B, duty);
}

void set_duty_fl(float duty)
{
    duty=(float)(DUTYMAX-DUTYMIN)*duty+DUTYMIN;
    if (duty>DUTYMAX-50)duty=DUTYMAX-50;
    if (duty<DUTYMIN+15)duty=DUTYMIN+15;
    pwm_set_chan_level(slice_num[1], PWM_CHAN_B, duty);
}

void set_duty_rr(float duty)
{
    duty=(float)(DUTYMAX-DUTYMIN)*duty+DUTYMIN;
    if (duty>DUTYMAX-50)duty=DUTYMAX-50;
    if (duty<DUTYMIN+15)duty=DUTYMIN+15;
    pwm_set_chan_level(slice_num[0], PWM_CHAN_A, duty);
}

void set_duty_rl(float duty)
{
    duty=(float)(DUTYMAX-DUTYMIN)*duty+DUTYMIN;
    if (duty>DUTYMAX-50)duty=DUTYMAX-50;
    if (duty<DUTYMIN+15)duty=DUTYMIN+15;
    pwm_set_chan_level(slice_num[1], PWM_CHAN_A, duty);
}


void on_uart_rx() {
  short data;
  while (uart_is_readable(UART_ID)) {
    ch = uart_getc(UART_ID);
    if(ch==0x0f&&chars_rxed==00){
      sbus_data[chars_rxed]=ch;
      //printf("%02X ",ch);
      chars_rxed++;
    }
    else if(chars_rxed>0){
      sbus_data[chars_rxed]=ch;
      //printf("%02X ",ch);
      chars_rxed++;            
    }
    //if (uart_is_writable(UART_ID)) {
    //    // Change it slightly first!
    //ch++;
    //    uart_putc(UART_ID, ch);
    //}

    switch(chars_rxed){
      case 3:
        Olddata[0]=(sbus_data[1]|(sbus_data[2]<<8)&0x07ff);
        Data1=(float)(Olddata[0]-CH1MID)*2/((CH1MAX-CH1MIN));
        //printf("%04d ",Olddata[0]);
        //printf("%04f ",Data1);
        break;
      case 4:
        Olddata[1]=(sbus_data[3]<<5|sbus_data[2]>>3)&0x07ff;	
        Data2=(float)(Olddata[1]-CH2MID)*2/((CH2MAX-CH2MIN));
        //printf("%04d ",Olddata[1]);
        //printf("%04f ",Data2);
        break;
      case 6:
        Olddata[2]=(sbus_data[3]>>6|sbus_data[4]<<2|sbus_data[5]<<10)&0x07ff;
        Data3=(float)(Olddata[2]-CH3MIN)/(CH3MAX-CH3MIN);
        //printf("%04d ",Olddata[2]);
        //printf("%04f ",Data3);
        break;
      case 7:
        Olddata[3]=(sbus_data[6]<<7|sbus_data[5]>>1)&0x07ff;
        Data4=(float)(Olddata[3]-CH4MID)*2/((CH4MAX-CH4MIN));
        //printf("%04d ",Olddata[3]);
        //printf("%04f ",Data4);
        break;
      case 8:
        Olddata[4]=(sbus_data[7]<<4|sbus_data[6]>>4)&0x07ff;
        Data5=(float)(Olddata[4]-CH5MID)*2/((CH5MAX-CH5MIN));
        //printf("%04d ",Olddata[4]);
        //printf("%04f ",Data5);
        break;
      case 10:
        Olddata[5]=(sbus_data[7]>>7|sbus_data[8]<<1|sbus_data[9]<<9)&0x07ff;
        Data6=(float)(Olddata[5]-CH6MID)*2/((CH6MAX-CH6MIN));
        //printf("%04d ",Olddata[5]);
        //printf("%04f ",Data6);
        break;
    }

    if(chars_rxed==25)
    {
      //printf("\n");
      chars_rxed=0;
    }
  }
}

