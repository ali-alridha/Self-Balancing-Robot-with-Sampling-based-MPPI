#include "MPU6050.h"
#include <AFMotor.h>
MPU6050 mpu;
AF_DCMotor motor1(1);
AF_DCMotor motor2(2);

int prev_time = millis();
int curr_time = millis();
float x = 0.0;
float x_dot;

void setup() {
  Wire.begin();
  Serial.begin(9600);
  mpu.initialize();     // запускаем датчик
  motor1.run(FORWARD);  // задаем движение вперед
  motor2.run(FORWARD);  // задаем движение вперед
}

void loop() {
  int16_t ax = mpu.getAccelerationX();  // ускорение по оси Х
  int16_t gyrY = mpu.getRotationY();
  // Serial.println(gyrY);
  float gyrY_f = gyrY * 250 / 32768.0;

  // стандартный диапазон: +-2g
  ax = constrain(ax, -16384, 16384);    // ограничиваем +-1g
  float angle = ax / 16384.0;           // переводим в +-1.0

  // и в угол через арксинус
  if (angle < 0) angle = 90 - degrees(acos(angle));
  else angle = degrees(acos(-angle)) - 90;

  int u = (int) ((angle-2) * 8.5 - gyrY_f * 4);
  u = constrain(u, -255, 255);

  x_dot = (u / 255.0) * 0.086; // under 5.5 volts

  if (u > 0)
  {
     motor1.run(FORWARD);  // задаем движение вперед
     motor2.run(FORWARD);  // задаем движение вперед
     motor1.setSpeed(u);   // задаем скорость движения
     motor2.setSpeed(u);   // задаем скорость движения 
  }
  else
  {
     motor1.run(BACKWARD);  // задаем движение вперед
     motor2.run(BACKWARD);  // задаем движение вперед
     motor1.setSpeed(-u);   // задаем скорость движения
     motor2.setSpeed(-u);   // задаем скорость движения 
  }

  curr_time = millis();
  x += x_dot * (curr_time - prev_time) / 1000;
  prev_time = curr_time;

  // Serial.print(x);
  // Serial.print(" ");
  // Serial.println(x_dot);

  Serial.print(angle);
  Serial.print(" ");
  Serial.print(gyrY_f);
  Serial.print(" ");
  Serial.println(u);
  delay(5);
}