#include <Arduino.h>
#include <math.h>

// Constants
const double g = 9.81;
const double mass_of_cart = 1.0;
const double mass_of_pole = 0.01;
const double length_of_pole = 2.0;
const double max_force_abs = 100.0;
const double delta_t = 0.02;
const int horizon_step_T = 100;
const int number_of_samples_K = 1000;
const double param_lambda = 50.0;
const double Sigma = 10.0;

// Weights for cost functions
double stage_cost_weight[] = {5.0, 10.0, 0.1, 0.1};
double terminal_cost_weight[] = {5.0, 10.0, 0.1, 0.1};

// State and control arrays
double state[] = {0.0, M_PI, 0.0, 0.0};
double u_prev[100] = {0};

// Utility functions
double clip(double value, double min_value, double max_value)
{
    return max(min_value, min(value, max_value));
}

double wrap_angle(double theta)
{
    return fmod(theta + M_PI, 2 * M_PI) - M_PI;
}

// Sigway dynamics
void update(double force)
{
    double x = state[0];
    double theta = state[1];
    double x_dot = state[2];
    double theta_dot = state[3];

    double temp = (force + mass_of_pole * length_of_pole * theta_dot * theta_dot * sin(theta)) / (mass_of_cart + mass_of_pole);
    double new_theta_ddot = (g * sin(theta) - cos(theta) * temp) /
                            (length_of_pole * (4.0 / 3.0 - mass_of_pole * cos(theta) * cos(theta) / (mass_of_cart + mass_of_pole)));
    double new_x_ddot = temp - mass_of_pole * length_of_pole * new_theta_ddot * cos(theta) / (mass_of_cart + mass_of_pole);

    state[0] += x_dot * delta_t;
    state[1] += theta_dot * delta_t;
    state[1] = wrap_angle(state[1]);
    state[2] += new_x_ddot * delta_t;
    state[3] += new_theta_ddot * delta_t;
}

// Cost functions
double stage_cost(double *x_t)
{
    double x = x_t[0];
    double theta = wrap_angle(x_t[1]);
    double x_dot = x_t[2];
    double theta_dot = x_t[3];

    return stage_cost_weight[0] * x * x + stage_cost_weight[1] * theta * theta + stage_cost_weight[2] * x_dot * x_dot + stage_cost_weight[3] * theta_dot * theta_dot;
}

double terminal_cost(double *x_T)
{
    double x = x_T[0];
    double theta = wrap_angle(x_T[1]);
    double x_dot = x_T[2];
    double theta_dot = x_T[3];

    return terminal_cost_weight[0] * x * x + terminal_cost_weight[1] * theta * theta + terminal_cost_weight[2] * x_dot * x_dot + terminal_cost_weight[3] * theta_dot * theta_dot;
}

// Main MPPI function
void mppi_control(double *observed_x, double &input_force)
{
    double S[number_of_samples_K] = {0};
    double epsilon[number_of_samples_K][horizon_step_T] = {0};
    double u[horizon_step_T] = {0};
    double x0[4];
    memcpy(x0, observed_x, 4 * sizeof(double));
    double param_gamma = param_lambda * (1.0 - 1.0);

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
        double v[horizon_step_T] = {0};
        double x[4];
        memcpy(x, x0, 4 * sizeof(double));

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

            update(v[t]);
            S[k] += stage_cost(state) + param_gamma * u_prev[t] * (1.0 / Sigma) * v[t];
        }

        S[k] += terminal_cost(state);
    }

    // Compute weights
    double rho = S[0];
    for (int k = 1; k < number_of_samples_K; ++k)
    {
        if (S[k] < rho)
        {
            rho = S[k];
        }
    }
    double eta = 0.0;
    for (int k = 0; k < number_of_samples_K; ++k)
    {
        eta += exp(-1.0 / param_lambda * (S[k] - rho));
    }

    double w[number_of_samples_K];
    for (int k = 0; k < number_of_samples_K; ++k)
    {
        w[k] = exp(-1.0 / param_lambda * (S[k] - rho)) / eta;
    }

    // Compute control input
    double w_epsilon[horizon_step_T] = {0};
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
    Serial.begin(115200);
    Serial.println("[INFO] Starting MPPI Controller for Sigway");
}

void loop()
{
    double input_force;
    mppi_control(state, input_force);

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

    update(input_force);

    delay(20); // Delay to simulate delta_t
}
