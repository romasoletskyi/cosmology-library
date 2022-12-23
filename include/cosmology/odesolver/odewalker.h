#pragma once

#include <cosmology/odesolver/odesystem.h>

namespace walker {
    template<class Value, class Time>
    class Walker {
    public:
        Walker(matrix::System<Value, Time> system, const std::vector<Time> &grid) : system_(system),
                                                                                    grid_(grid) {
        }

    protected:
        matrix::System<Value, Time> system_;
        const std::vector<Time> &grid_;
    };

    /*
     * Backward Euler discretizes system as:
     * (x^{n+1} - x^n) / dt = A^{n+1}x^{n+1} + b^{n+1} <=> (1 - dt A^{n+1}) x^{n+1} = x^n + dt b^{n+1}
     */
    template<class Value, class Time>
    class BackwardEulerWalker : public Walker<Value, Time> {
    public:
        using Walker<Value, Time>::Walker;

        void next(typename OdeSolver<Value, Time>::Solution &solution, size_t step) {
            using Mat = typename matrix::System<Value, Time>::Mat;

            Time time = this->grid_[step + 1];
            Time dt = time - this->grid_[step];
            Mat mat = Mat::Identity(solution.variableNumber, solution.variableNumber) -
                      dt * this->system_.getMatrix(time, solution.getState(step));
            Mat source = Eigen::Map<const Mat>(solution.getState(step), solution.variableNumber, 1) +
                         dt * this->system_.getSource(time, solution.getState(step));
            Mat nextMatState = mat.colPivHouseholderQr().solve(source);

            for (size_t j = 0; j < solution.variableNumber; ++j) {
                solution(step + 1, j) = nextMatState(j, 0);
            }
        }
    };

    // Standard 4-step Runge-Kutta method
    template<class Value, class Time>
    class RungeKuttaWalker : public Walker<Value, Time> {
    public:
        using Walker<Value, Time>::Walker;

        void next(typename OdeSolver<Value, Time>::Solution &solution, size_t step) {
            Time time = this->grid_[step];
            Time dt = this->grid_[step + 1] - time;

            Mat matState = Eigen::Map<const Mat>(solution.getState(step), solution.variableNumber, 1);
            Mat k1 = getDerivative(matState, time);
            Mat k2 = getDerivative(matState + k1 * dt / 2, time + dt / 2);
            Mat k3 = getDerivative(matState + k2 * dt / 2, time + dt / 2);
            Mat k4 = getDerivative(matState + k3 * dt, time + dt);

            Mat nextMatState = matState + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            for (size_t j = 0; j < solution.variableNumber; ++j) {
                solution(step + 1, j) = nextMatState(j, 0);
            }
        }

    private:
        using Mat = typename matrix::System<Value, Time>::Mat;

        auto getDerivative(const Mat &state, Time time) {
            return this->system_.getMatrix(time, state.data()) * state + this->system_.getSource(time, state.data());
        }
    };
} // namespace walker