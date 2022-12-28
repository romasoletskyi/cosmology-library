#pragma once

#include <cosmology/odesolver/odesystem.h>

template<class Value>
struct Solution;

namespace walker {
    template<class Value, class Time>
    class IWalker {
    protected:
        using Matrix = typename formula::System<Value, Time>::Matrix;

    public:
        explicit IWalker(formula::GeneralSystem<Value, Time> system) : system_(std::move(system)) {
        }

        virtual ~IWalker() = default;

        virtual Matrix nextDynamic(Solution<Value> &solution, const std::vector<Time> &grid, int step) = 0;

        // solution at step' <= step is known, calculates solution at step + 1
        void next(Solution<Value> &solution, const std::vector<Time> &grid, int step) {
            Time time = grid[step + 1];

            auto nextDynamicState = nextDynamic(solution, grid, step);
            auto nextNonDynamicState = system_.nondynamic.getMatrix(time, solution.getState(step)) * nextDynamicState +
                                       system_.nondynamic.getSource(time, solution.getState(step));

            system_.port.storeTotalState(nextDynamicState, nextNonDynamicState, solution.getState(step + 1));
        }

    protected:
        formula::GeneralSystem<Value, Time> system_;
    };

    /*
     * Backward Euler discretizes system as:
     * (x^{n+1} - x^n) / dt = A^{n+1}x^{n+1} + b^{n+1} <=> (1 - dt A^{n+1}) x^{n+1} = x^n + dt b^{n+1}
     */
    template<class Value, class Time>
    class BackwardEulerWalker : public IWalker<Value, Time> {
    private:
        using Matrix = typename IWalker<Value, Time>::Matrix;

    public:
        using IWalker<Value, Time>::IWalker;

        ~BackwardEulerWalker() = default;

        Matrix nextDynamic(Solution<Value> &solution, const std::vector<Time> &grid, int step) {
            Time time = grid[step + 1];
            Time dt = time - grid[step];
            int dynamicVariablesNumber = this->system_.port.getDynamicVariablesNumber();

            Matrix matrix = Matrix::Identity(dynamicVariablesNumber, dynamicVariablesNumber) -
                            dt * this->system_.dynamic.getMatrix(time, solution.getState(step));
            Matrix source = this->system_.port.loadDynamicState(solution.getState(step)) +
                            dt * this->system_.dynamic.getSource(time, solution.getState(step));

            return matrix.colPivHouseholderQr().solve(source);
        }
    };

    // Standard 4-step Runge-Kutta method
    template<class Value, class Time>
    class RungeKuttaWalker : public IWalker<Value, Time> {
    private:
        using Matrix = typename IWalker<Value, Time>::Matrix;

    public:
        using IWalker<Value, Time>::IWalker;

        ~RungeKuttaWalker() = default;

        Matrix nextDynamic(Solution<Value> &solution, const std::vector<Time> &grid, int step) {
            Time time = grid[step];
            Time dt = grid[step + 1] - time;

            Matrix dynamicState = this->system_.port.loadDynamicState(solution.getState(step));

            Matrix k1 = getDerivative(dynamicState, solution.getState(step), time);
            Matrix k2 = getDerivative(dynamicState + k1 * dt / 2, solution.getState(step), time + dt / 2);
            Matrix k3 = getDerivative(dynamicState + k2 * dt / 2, solution.getState(step), time + dt / 2);
            Matrix k4 = getDerivative(dynamicState + k3 * dt, solution.getState(step), time + dt);

            return dynamicState + dt * (k1 + 2 * k2 + 2 * k3 + k4) / 6;
        }

    private:
        auto getDerivative(const Matrix &dynamicState, const Value *totalState, Time time) {
            return this->system_.dynamic.getMatrix(time, totalState) * dynamicState +
                   this->system_.dynamic.getSource(time, totalState);
        }
    };
} // namespace walker