#ifndef PATH_SMOOTHER_HPP
#define PATH_SMOOTHER_HPP

#include "cubic_spline.hpp"
#include "lbfgs.hpp"

#include <Eigen/Eigen>

#include <cmath>
#include <cfloat>
#include <iostream>
#include <vector>
#include <algorithm>

#define MAX(a, b) (a > b)?a:b

namespace path_smoother
{

    class PathSmoother
    {
    private:
        cubic_spline::CubicSpline cubSpline;

        int pieceN;
        Eigen::Matrix3Xd diskObstacles;
        double penaltyWeight;
        Eigen::Vector2d headP;
        Eigen::Vector2d tailP;
        Eigen::Matrix2Xd points;
        Eigen::Matrix2Xd gradByPoints;

        Eigen::MatrixXd inv_A;
        Eigen::MatrixXd invA_dB;
        Eigen::MatrixXd dx;

        unsigned int _m;
        std::vector<Eigen::VectorXd> x_buff;
        std::vector<Eigen::VectorXd> g_buff;

        lbfgs::lbfgs_parameter_t lbfgs_params;

    private:
        static inline double costFunction(void *ptr,
                                          const Eigen::VectorXd &x,
                                          Eigen::VectorXd &g)
        {
           //TODO
            return 0.0f;
        }

        int line_search_lewisoverton_easy(Eigen::VectorXd &x,
                                        double &f,
                                        Eigen::VectorXd &g,
                                        double &stp,
                                        const Eigen::VectorXd &s,
                                        const Eigen::VectorXd &xp,
                                        const Eigen::VectorXd &gp,
                                        const double stpmin,
                                        const double stpmax,
                                        const lbfgs::lbfgs_parameter_t &param,
                                        CubicCurve &curve)
    {
        // x is the decision variable vector
        // f is function value at x
        // g is the gradient value at x
        // stp is the initial stepsize for line search 
        // s is the search direction vector
        // xp is the decision variable vector at the current iteration
        // gp is the gradient vector at the current iteration
        // stpmin is the minimum allowable stepsize
        // stpmax is the maximum allowable stepsize
        // the struct param contains all necessary parameters 
        // the cd contains all necessary callback function

        // eg.             x = xp; f = cd.proc_evaluate(cd.instance, x, g);
        // the above line assigns x with xp and computes the function and grad at x
        
        // note the output x, f and g which satisfy the weak wolfe condition when the function returns

        //////////////////////////// HOMEWORK 1 START ////////////////////////////

        // PUT YOUR CODE FOR Lewis-Overton line search here
        int cnt = 0;
        if(stp <= 0.0f)
        {
            return lbfgs::LBFGSERR_INVALID_FUNCVAL;
        }
        double l = stpmin;
        double u = stpmax + 1e-5f;
        stp = (stp > stpmin)? stp : stpmin;
        stp = (stp < stpmax)? stp : stpmax;
        double fk = f;
        double cost_loss_k = s.dot(gp);
        if(cost_loss_k > 0.0f)
        {
            return lbfgs::LBFGSERR_INCREASEGRADIENT;
        }
        
        while (true)
        {
            cnt++;
            x = xp + stp * s;
            f = getGrad(x, g, curve);
            if (std::isinf(f) || std::isnan(f))
            {
                return lbfgs::LBFGSERR_INVALID_FUNCVAL;
            }
            // sufficient decrease condition
            if(fk - f < -param.f_dec_coeff * stp * cost_loss_k)
            {
                u = stp;
            }
            // curvature condition
            else if(s.dot(g) < param.s_curv_coeff * cost_loss_k)
                l = stp;
            else
                return cnt;

            if(cnt > param.max_linesearch)
            {
                return lbfgs::LBFGSERR_MAXIMUMLINESEARCH;
            }

            if((u - l) < 2.0f * param.machine_prec)
            {
                return lbfgs::LBFGSERR_WIDTHTOOSMALL;
            }
            if(u < stpmax)
                stp = 0.5f * (l + u);
            else
            {
                stp = 2 * l;
                if(stp > stpmax)
                {
                    return lbfgs::LBFGSERR_MAXIMUMSTEP;
                }
            }
        }
    }

    public:
        inline bool setup(const Eigen::Vector2d &initialP,
                          const Eigen::Vector2d &terminalP,
                          const int &pieceNum,
                          const Eigen::Matrix3Xd &diskObs,
                          const double penaWeight)
        {
            _m = 8;
            pieceN = pieceNum;
            diskObstacles = diskObs;
            penaltyWeight = penaWeight;
            headP = initialP;
            tailP = terminalP;

            cubSpline.setConditions(headP, tailP, pieceN);

            points.resize(2, pieceN - 1);
            gradByPoints.resize(2, pieceN - 1);

            inv_A = 4.0f * inv_A.setIdentity(pieceN-1, pieceN-1);
            invA_dB = Eigen::MatrixXd::Zero(pieceN, pieceN - 1);
            Eigen::MatrixXd dB = Eigen::MatrixXd::Zero(pieceN - 1, pieceN - 1);
            dx = Eigen::MatrixXd::Zero(pieceN, pieceN - 1);
            for(int i = 0; i < pieceN-2; i++)
            {
                inv_A(i, i+1) = 1.0f;
                inv_A(i+1, i) = 1.0f;
                dB(i, i+1) = 3.0f;
                dB(i+1, i) = -3.0f;
                dx(i, i) = -1.0f;
                dx(i+1, i) = 1.0f;
            }
            dx(pieceN-2, pieceN-2) = -1.0f;
            dx(pieceN-1, pieceN-2) = 1.0f;
            inv_A = inv_A.inverse();
            invA_dB.block(0, 0, pieceN - 1, pieceN - 1) = inv_A * dB;

            return true;
        }

        void matrix2vec(Eigen::VectorXd &output, const Eigen::Matrix2Xd &input)
        {
            output.resize(input.rows() * input.cols());
            for(int i = 0; i < input.cols(); i++)
            {
                output.block(i * input.rows(), 0, input.rows(), 1) = input.col(i);
            }
        }
        void vec2matrix(Eigen::Matrix2Xd &output, const Eigen::VectorXd &input)
        {
            output.resize(2, int(input.size()/2));
            for(int i = 0; i < output.cols(); i++)
            {
                output.col(i) = input.block(2*i, 0, 2, 1);
            }
        }

        double getGrad(const Eigen::VectorXd &in, Eigen::VectorXd &gout, CubicCurve &curve)
        {
            Eigen::Matrix2Xd g_points;
            double cost = 0.0f;
            vec2matrix(g_points, in);
            gradByPoints = Eigen::MatrixXd::Zero(2, pieceN - 1);
            
            Eigen::MatrixXd tmp_x;
            tmp_x.resize(2, pieceN+1);
            tmp_x.col(0) = headP;
            tmp_x.col(pieceN) = tailP;
            tmp_x.block(0, 1, 2, pieceN-1) = g_points;

            Eigen::MatrixXd diff_x = 3.0f * (tmp_x.block(0, 2, 2, pieceN-1) - tmp_x.block(0, 0, 2, pieceN-1)).transpose();
            Eigen::MatrixXd D = Eigen::MatrixXd::Zero(pieceN, 2);
            D.block(0, 0, pieceN-1, 2) = inv_A * diff_x;

            Eigen::Matrix<double, 2, 4> cur;
            curve.clear();
            for (int i = 0; i < pieceN-1; i++)
            {
                Eigen::MatrixXd dd = 2.0f * dx.row(i) + invA_dB.row(i) + invA_dB.row(i+1);
                Eigen::MatrixXd dc = -3.0f * dx.row(i) - 2.0f * invA_dB.row(i) - invA_dB.row(i+1);
                Eigen::Vector2d c = 3.0f * (tmp_x.col(i+1) - tmp_x.col(i)) - 2.0f * D.row(i).transpose() - D.row(i+1).transpose();
                Eigen::Vector2d d = 2.0f * (tmp_x.col(i) - tmp_x.col(i+1)) + D.row(i).transpose()+ D.row(i+1).transpose();
                cur.col(3) = tmp_x.col(i);
                cur.col(2) = D.row(i).transpose();
                cur.col(1) = c;
                cur.col(0) = d;
                double tmp = 1.0f;
                curve.emplace_back(tmp, cur);
                cost = cost + 12.0f * d.dot(d) + 6.0f * d.dot(c) + 4.0f * c.dot(c);
                gradByPoints = gradByPoints + 24.0f * d * dd + 12.0f * d * dc + 12.0f * c * dd + 8.0f * c * dc;
            }

            for (int i = 0; i < pieceN-1; i++)
            {
                for (int j = 0; j < diskObstacles.cols(); j++)
                {
                    double dist_2 = (g_points(0, i) - diskObstacles(0, j)) * (g_points(0, i) - diskObstacles(0, j)) + 
                                    (g_points(1, i) - diskObstacles(1, j)) * (g_points(1, i) - diskObstacles(1, j));
                    if(dist_2 < diskObstacles(2, j) * diskObstacles(2, j))
                    {
                        gradByPoints(0, i) -= penaltyWeight * (g_points(0, i) - diskObstacles(0, j)) / sqrt(dist_2);
                        gradByPoints(1, i) -= penaltyWeight * (g_points(1, i) - diskObstacles(1, j)) / sqrt(dist_2);
                        cost = cost + penaltyWeight * (diskObstacles(2, j) - sqrt(dist_2));
                    }
                }
            }
            matrix2vec(gout, gradByPoints);
            return cost;
        }

        double optimize(CubicCurve &curve,
                               const Eigen::Matrix2Xd &iniInPs,
                               const double &relCostTol)
        {
            //TODO
            points = iniInPs;
            double cost = 100.0f;
            double last_cost = -1.0f;
            int cnt = 0;
            Eigen::VectorXd pointsvec;
            Eigen::VectorXd gradvec;
            while((fabs(cost - last_cost)/MAX(1.0f, fabs(cost))) > lbfgs_params.delta && cnt < 1000)
            {
                
                matrix2vec(pointsvec, points);
                last_cost = cost;
                cost = getGrad(pointsvec, gradvec, curve);
                if(x_buff.size() <= _m)
                {
                    x_buff.push_back(pointsvec);
                    g_buff.push_back(gradvec);
                }
                else
                {
                    Eigen::VectorXd delta_x = pointsvec - x_buff.back();
                    Eigen::VectorXd delta_g = gradvec - g_buff.back();
                    if(delta_g.dot(delta_x) > 1e-6f * gradvec.norm() * delta_x.dot(delta_x))
                    {
                        x_buff.erase(x_buff.begin());
                        g_buff.erase(g_buff.begin());
                        x_buff.push_back(pointsvec);
                        g_buff.push_back(gradvec);
                    }
                }
                Eigen::VectorXd d = gradvec;
                std::vector<double> alpha;
                std::vector<Eigen::VectorXd> s;
                std::vector<Eigen::VectorXd> y;
                for(int i = x_buff.size()-1; i > 0; i--)
                {
                    Eigen::VectorXd tmp_s = x_buff[i] - x_buff[i-1];
                    s.insert(s.begin(), tmp_s);
                    Eigen::VectorXd tmp_y = g_buff[i] - g_buff[i-1];
                    y.insert(y.begin(), tmp_y);
                    double tmp_alpha = tmp_s.dot(d)/tmp_s.dot(tmp_y);
                    alpha.insert(alpha.begin(), tmp_alpha);
                    d = d - tmp_alpha * tmp_y;
                }
                if(y.size() > 0)
                {
                    double gamma = y.back().dot(y.back())/s.back().dot(y.back());
                    d = d/gamma;
                }
                for(unsigned int i = 0; i < s.size(); i++)
                {
                    double beta = y[i].dot(d)/s[i].dot(y[i]);
                    d = d + (alpha[i] - beta) * s[i];
                }
                d = -d;
                Eigen::VectorXd out_pointsvec;
                Eigen::VectorXd out_grad;
                double step = 1.0f;
                int err = line_search_lewisoverton_easy(out_pointsvec, cost, out_grad, step, d, pointsvec
                                              , gradvec, 0.0f, 1000.0f, lbfgs_params, curve);
                if(err < 0)
                {
                    std::cout << "Solver is error: " << err << std::endl;
                    break;
                }
                vec2matrix(points, out_pointsvec);
                cnt++;
            }
            std::cout << "cost is: " << cost << " iterate " << cnt << " times." << std::endl;
            cost = getGrad(pointsvec, gradvec, curve);
            return cost;
        }
    };

}

#endif
