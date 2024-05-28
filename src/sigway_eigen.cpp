#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <algorithm>
#include <numeric>
#include <eigen3/Eigen/Dense>

using namespace Eigen;


double clip(double value, double min_value, double max_value)
{
    return std::max(min_value, std::min(value, max_value));
}


class Sigway
{
public:
    Sigway(double mass_of_cart = 1.0, double mass_of_pole = 0.01, double length_of_pole = 2.0, double max_force_abs = 100.0, double delta_t = 0.02, bool visualize = true): mass_of_cart(mass_of_cart), mass_of_pole(mass_of_pole), length_of_pole(length_of_pole), max_force_abs(max_force_abs), delta_t(delta_t), visualize_flag(visualize)
    {
        g = 9.81;
        reset();
    }

    void reset(const VectorXd &init_state = (VectorXd(4) << 0.0, M_PI, 0.0, 0.0).finished())
    {
        state = init_state;
    }

    void update(const VectorXd &u, double delta_t = 0.0)
    {
        double x = state[0];
        double theta = state[1];
        double x_dot = state[2];
        double theta_dot = state[3];

        double dt = (delta_t == 0.0) ? this->delta_t : delta_t;

        double force = clip(u[0], -max_force_abs, max_force_abs);

        double temp = (force + mass_of_pole * length_of_pole * theta_dot * theta_dot * sin(theta)) / (mass_of_cart + mass_of_pole);
        double new_theta_ddot = (g * sin(theta) - cos(theta) * temp) /
                                (length_of_pole * (4.0 / 3.0 - mass_of_pole * cos(theta) * cos(theta) / (mass_of_cart + mass_of_pole)));
        double new_x_ddot = temp - mass_of_pole * length_of_pole * new_theta_ddot * cos(theta) / (mass_of_cart + mass_of_pole);

        double new_x = x + x_dot * dt;
        double new_theta = theta + theta_dot * dt;
        new_theta = fmod(new_theta + M_PI, 2 * M_PI) - M_PI;

        double new_theta_dot = theta_dot + new_theta_ddot * dt;
        double new_x_dot = x_dot + new_x_ddot * dt;

        state << new_x, new_theta, new_x_dot, new_theta_dot;
    }

    VectorXd get_state() const
    {
        return state;
    }

private:
    double g;
    double mass_of_cart;
    double mass_of_pole;
    double length_of_pole;
    double max_force_abs;
    double delta_t;
    bool visualize_flag;


    VectorXd state;
};


class MPPIControllerForCartPole
{
public:
    MPPIControllerForCartPole(double delta_t = 0.02, double mass_of_cart = 1.0, double mass_of_pole = 0.01,
                              double length_of_pole = 2.0, double max_force_abs = 100.0, int horizon_step_T = 100,
                              int number_of_samples_K = 1000, double param_exploration = 0.0, double param_lambda = 50.0,
                              double param_alpha = 1.0, double sigma = 10.0,
                              const VectorXd &stage_cost_weight = (VectorXd(4) << 5.0, 10.0, 0.1, 0.1).finished(),
                              const VectorXd &terminal_cost_weight = (VectorXd(4) << 5.0, 10.0, 0.1, 0.1).finished())
        : delta_t(delta_t), mass_of_cart(mass_of_cart), mass_of_pole(mass_of_pole), length_of_pole(length_of_pole),
          max_force_abs(max_force_abs), T(horizon_step_T), K(number_of_samples_K), param_exploration(param_exploration),
          param_lambda(param_lambda), param_alpha(param_alpha), param_gamma(param_lambda * (1.0 - param_alpha)),
          Sigma(sigma), stage_cost_weight(stage_cost_weight), terminal_cost_weight(terminal_cost_weight)
    {
        u_prev = VectorXd::Zero(T);
    }

    std::pair<double, VectorXd> calc_control_input(const VectorXd &observed_x)
    {
        VectorXd u = u_prev;
        VectorXd x0 = observed_x;
        VectorXd S = VectorXd::Zero(K);
        MatrixXd epsilon = calc_epsilon(Sigma, K, T);

        for (int k = 0; k < K; ++k)
        {
            VectorXd v = VectorXd::Zero(T);
            VectorXd x = x0;

            for (int t = 1; t <= T; ++t)
            {
                if (k < (1.0 - param_exploration) * K)
                {
                    v[t - 1] = u[t - 1] + epsilon(k, t - 1);
                }
                else
                {
                    v[t - 1] = epsilon(k, t - 1);
                }

                x = F(x, g(v[t - 1]));
                S[k] += c(x) + param_gamma * u[t - 1] * (1.0 / Sigma) * v[t - 1];
            }

            S[k] += phi(x);
        }

        VectorXd w = compute_weights(S);

        VectorXd w_epsilon = VectorXd::Zero(T);
        for (int t = 0; t < T; ++t)
        {
            for (int k = 0; k < K; ++k)
            {
                w_epsilon[t] += w[k] * epsilon(k, t);
            }
        }

        w_epsilon = moving_average_filter(w_epsilon, 10);

        u += w_epsilon;

        u_prev.head(T - 1) = u.tail(T - 1);
        u_prev[T - 1] = u[T - 1];

        return {u[0], u};
    }

private:
    double delta_t;
    double mass_of_cart;
    double mass_of_pole;
    double length_of_pole;
    double max_force_abs;
    int T;
    int K;
    double param_exploration;
    double param_lambda;
    double param_alpha;
    double param_gamma;
    double Sigma;
    VectorXd stage_cost_weight;
    VectorXd terminal_cost_weight;

    VectorXd u_prev;

    MatrixXd calc_epsilon(double sigma, int size_sample, int size_time_step)
    {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::normal_distribution<> d(0, sigma);

        MatrixXd epsilon(size_sample, size_time_step);
        for (int i = 0; i < size_sample; ++i)
        {
            for (int j = 0; j < size_time_step; ++j)
            {
                epsilon(i, j) = d(gen);
            }
        }
        return epsilon;
    }

    double g(double v)
    {
        return clip(v, -max_force_abs, max_force_abs);
    }

    double c(const VectorXd &x_t)
    {
        double x = x_t[0];
        double x_dot = x_t[2];
        double theta = x_t[1];
        double theta_dot = x_t[3];
        theta = fmod(theta + M_PI, 2 * M_PI) - M_PI;

        return stage_cost_weight[0] * x * x + stage_cost_weight[1] * theta * theta +
               stage_cost_weight[2] * x_dot * x_dot + stage_cost_weight[3] * theta_dot * theta_dot;
    }

    double phi(const VectorXd &x_T)
    {
        double x = x_T[0];
        double x_dot = x_T[2];
        double theta = x_T[1];
        double theta_dot = x_T[3];
        theta = fmod(theta + M_PI, 2 * M_PI) - M_PI;

        return terminal_cost_weight[0] * x * x + terminal_cost_weight[1] * theta * theta +
               terminal_cost_weight[2] * x_dot * x_dot + terminal_cost_weight[3] * theta_dot * theta_dot;
    }

    VectorXd F(const VectorXd &x_t, double v_t)
    {
        double x = x_t[0];
        double theta = x_t[1];
        double x_dot = x_t[2];
        double theta_dot = x_t[3];

        double temp = (v_t + mass_of_pole * length_of_pole * theta_dot * theta_dot * sin(theta)) / (mass_of_cart + mass_of_pole);
        double new_theta_ddot = (g * sin(theta) - cos(theta) * temp) /
                                (length_of_pole * (4.0 / 3.0 - mass_of_pole * cos(theta) * cos(theta) / (mass_of_cart + mass_of_pole)));
        double new_x_ddot = temp - mass_of_pole * length_of_pole * new_theta_ddot * cos(theta) / (mass_of_cart + mass_of_pole);

        theta += theta_dot * delta_t;
        x += x_dot * delta_t;
        theta_dot += new_theta_ddot * delta_t;
        x_dot += new_x_ddot * delta_t;

        return (VectorXd(4) << x, theta, x_dot, theta_dot).finished();
    }

    VectorXd compute_weights(const VectorXd &S)
    {
        VectorXd w(K);
        double rho = S.minCoeff();
        double eta = 0.0;

        for (int k = 0; k < K; ++k)
        {
            eta += exp((-1.0 / param_lambda) * (S[k] - rho));
        }

        for (int k = 0; k < K; ++k)
        {
            w[k] = (1.0 / eta) * exp((-1.0 / param_lambda) * (S[k] - rho));
        }

        return w;
    }

    VectorXd moving_average_filter(const VectorXd &xx, int window_size)
    {
        VectorXd b = VectorXd::Ones(window_size) / window_size;
        VectorXd xx_mean = VectorXd::Zero(xx.size());

        for (int i = 0; i < xx.size(); ++i)
        {
            int start = std::max(0, i - window_size / 2);
            int end = std::min(static_cast<int>(xx.size()), i + window_size / 2 + 1);
            xx_mean[i] = xx.segment(start, end - start).mean();
        }

        return xx_mean;
    }
};


int main()
{
    double delta_t = 0.02;
    int sim_steps = 200;

    std::cout << "[INFO] delta_t: " << delta_t << "[s], sim_steps: " << sim_steps
              << "[steps], total_sim_time: " << delta_t * sim_steps << "[s]" << std::endl;

    
    Sigway cartpole(1.0, 0.01, 2.0, 100.0);
    cartpole.reset((VectorXd(4) << 0.0, M_PI, 0.0, 0.0).finished());

    
    MPPIControllerForCartPole mppi(delta_t, 1.0, 0.01, 2.0, 100.0, 100, 1000, 0.0, 50.0, 1.0, 10.0,
                                   (VectorXd(4) << 5.0, 10.0, 0.1, 0.1).finished(),
                                   (VectorXd(4) << 5.0, 10.0, 0.1, 0.1).finished());

    
    for (int i = 0; i < sim_steps; ++i)
    {
        VectorXd current_state = cartpole.get_state();

        auto [input_force, input_force_sequence] = mppi.calc_control_input(current_state);

        std::cout << "Time: " << i * delta_t << "[s], x=" << current_state[0] << "[m], theta=" << current_state[1]
                  << "[rad], input force=" << input_force << "[N]" << std::endl;

        cartpole.update(VectorXd::Constant(1, input_force), delta_t);
    }

    return 0;
}
