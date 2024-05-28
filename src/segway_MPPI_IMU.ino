#include <Arduino.h>
#include <math.h>
#include "MPU6050.h"
#include <AFMotor.h>


// Constants
const float g = 9.81;
const float mass_of_cart = 0.074;
const float mass_of_pole = 0.38;
const float length_of_pole = 0.15;
const float max_force_abs = 100.0;
const float delta_t = 0.02;
const int horizon_step_T = 20;
const int number_of_samples_K = 5;
const float param_lambda = 1.2;
const float Sigma = 10.8;

// constraints for state[0] (linear displacement) and state[2] (linear velocity)
const float x_min = -0.3;
const float x_max = 0.3;
const float x_dot_min = -0.3;
const float x_dot_max = 0.3;

// Weights for cost functions
float stage_cost_weight[] = {20.0, 20.0, 0.1, 0.1};
float terminal_cost_weight[] = {20.0, 20.0, 0.1, 0.1};

// State and control arrays
float state[] = {0.0, M_PI/2, 0.0, 0.0};
float u_prev[100] = {0};

int prev_time = millis();
int curr_time = millis();
float x = 0.0;
float x_dot = 0.0;

MPU6050 mpu;
AF_DCMotor motor1(1);
AF_DCMotor motor2(2);

// Utility functions
float clip(float value, float min_value, float max_value)
{
    return max(min_value, min(value, max_value));
}

float wrap_angle(float theta)
{
    return fmod(theta + M_PI, 2 * M_PI) - M_PI;
}

// CartPole dynamics
void update_cartpole(float force)
{
    float x = state[0];
    float theta = state[1];
    float x_dot = state[2];
    float theta_dot = state[3];

    float temp = (force + mass_of_pole * length_of_pole * theta_dot * theta_dot * sin(theta)) / (mass_of_cart + mass_of_pole);
    float new_theta_ddot = (g * sin(theta) - cos(theta) * temp) /
                            (length_of_pole * (4.0 / 3.0 - mass_of_pole * cos(theta) * cos(theta) / (mass_of_cart + mass_of_pole)));
    float new_x_ddot = temp - mass_of_pole * length_of_pole * new_theta_ddot * cos(theta) / (mass_of_cart + mass_of_pole);

    state[0] += x_dot * delta_t;
    state[1] += theta_dot * delta_t;
    state[1] = wrap_angle(state[1]);
    state[2] += new_x_ddot * delta_t;
    state[3] += new_theta_ddot * delta_t;
    // bound the displacement and linear velocity to upper and lower limits
    state[0] = constrain(state[0], x_min, x_max);
    state[2] = constrain(state[2], x_dot_min, x_dot_max);
}

// Cost functions
float stage_cost(float *x_t)
{
    float x = x_t[0];
    float theta = wrap_angle(x_t[1]);
    float x_dot = x_t[2];
    float theta_dot = x_t[3];

    return stage_cost_weight[0] * x * x + stage_cost_weight[1] * theta * theta + stage_cost_weight[2] * x_dot * x_dot + stage_cost_weight[3] * theta_dot * theta_dot;
}

float terminal_cost(float *x_T)
{
    float x = x_T[0];
    float theta = wrap_angle(x_T[1]);
    float x_dot = x_T[2];
    float theta_dot = x_T[3];

    return terminal_cost_weight[0] * x * x + terminal_cost_weight[1] * theta * theta + terminal_cost_weight[2] * x_dot * x_dot + terminal_cost_weight[3] * theta_dot * theta_dot;
}

// Main MPPI function
void mppi_control(float *observed_x, float &input_force)
{
    float S[number_of_samples_K] = {0};
    float epsilon[number_of_samples_K][horizon_step_T] = {0};
    float u[horizon_step_T] = {0};
    float x0[4];
    memcpy(x0, observed_x, 4 * sizeof(float));
    float param_gamma = param_lambda * (1.0 - 1.0);

    // Generate epsilon
    for (int k = 0; k < number_of_samples_K; ++k)
    {
        for (int t = 0; t < horizon_step_T; ++t)
        {
            epsilon[k][t] = random(-Sigma * 1000, Sigma * 1000) / 1000.0; // Random normal distribution
        }
    }

    // Cost calculation
    for (int k = 0; k < number_of_samples_K; ++k)
    {
        float v[horizon_step_T] = {0};
        float x[4];
        memcpy(x, x0, 4 * sizeof(float));

        for (int t = 0; t < horizon_step_T; ++t)
        {
            if (k < (1.0 - 0.0) * number_of_samples_K)
            {
                v[t] = u_prev[t] + epsilon[k][t];
            }
            else
            {
                v[t] = epsilon[k][t];
            }

            update_cartpole(v[t]);
            S[k] += stage_cost(state) + param_gamma * u_prev[t] * (1.0 / Sigma) * v[t];
        }

        S[k] += terminal_cost(state);
    }

    // Compute weights
    float rho = S[0];
    for (int k = 1; k < number_of_samples_K; ++k)
    {
        if (S[k] < rho)
        {
            rho = S[k];
        }
    }
    float eta = 0.0;
    for (int k = 0; k < number_of_samples_K; ++k)
    {
        eta += exp(-1.0 / param_lambda * (S[k] - rho));
    }

    float w[number_of_samples_K];
    for (int k = 0; k < number_of_samples_K; ++k)
    {
        w[k] = exp(-1.0 / param_lambda * (S[k] - rho)) / eta;
    }

    // Compute control input
    float w_epsilon[horizon_step_T] = {0};
    for (int t = 0; t < horizon_step_T; ++t)
    {
        for (int k = 0; k < number_of_samples_K; ++k)
        {
            w_epsilon[t] += w[k] * epsilon[k][t];
        }
    }

    for (int t = 0; t < horizon_step_T; ++t)
    {
        u[t] += w_epsilon[t];
    }
 
    for (int i = 0; i < horizon_step_T - 1; ++i)
    {
        u_prev[i] = u_prev[i + 1];
    }
    u_prev[horizon_step_T - 1] = u[horizon_step_T - 1];

    input_force = u[0];
}

// Simulation loop
void setup()
{
    Serial.begin(9600);
    Wire.begin();
    mpu.initialize();     // запускаем датчик
    motor1.run(FORWARD);  // задаем движение вперед
    motor2.run(FORWARD);  // задаем движение вперед
    Serial.println("[INFO] Starting MPPI Controller for CartPole");
}

void loop()
{
    int16_t ax = mpu.getAccelerationX();  // ускорение по оси Х
    int16_t gyrY = mpu.getRotationY();
    float gyrY_f = gyrY / 32768.0 * 250 * 3.14 / 180;

    // стандартный диапазон: +-2g
    ax = constrain(ax, -16384, 16384);    // ограничиваем +-1g
    float angle = ax / 16384.0;           // переводим в +-1.0
  
    angle = acos(angle);

    state[0] = x;
    state[1] = angle;
    state[2] = x_dot;
    state[3] = gyrY_f;

    float input_force;

    mppi_control(state, input_force);
    input_force *= (255 / 10.0 * 1.5);

    input_force = constrain(input_force, -255, 255);

  if (input_force > 0)
  {
     motor1.run(FORWARD);  // задаем движение вперед
     motor2.run(FORWARD);  // задаем движение вперед
     motor1.setSpeed(input_force);   // задаем скорость движения
     motor2.setSpeed(input_force);   // задаем скорость движения 
  }
  else
  {
     motor1.run(BACKWARD);  // задаем движение вперед
     motor2.run(BACKWARD);  // задаем движение вперед
     motor1.setSpeed(-input_force);   // задаем скорость движения
     motor2.setSpeed(-input_force);   // задаем скорость движения 
  }

    x_dot = (input_force / 255.0) * 0.086; // under 5.5 volts

    curr_time = millis();
    x += x_dot * (curr_time - prev_time) / 1000;
    prev_time = curr_time;

    Serial.print("x: ");
    Serial.print(state[0]);
    Serial.print(", theta: ");
    Serial.print(state[1]);
    Serial.print(", x_dot: ");
    Serial.print(state[2]);
    Serial.print(", theta_dot: ");
    Serial.print(state[3]);
    Serial.print(", input force: ");
    Serial.println(input_force);

    update_cartpole(input_force);

    delay(20); // Delay to simulate delta_t
}
