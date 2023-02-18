#pragma once

#include <functional>
#include <cmath>
#include <random>

template<class Value, class Parameters>
class Sampler {
public:
    Sampler(std::vector<Value> values, std::vector<Value> errors,
            std::function<std::vector<Value>(const Parameters &)> prediction,
            std::function<double(const Parameters &)> prior, std::function<Parameters(const Parameters &, std::mt19937& gen)> next,
            Parameters start) :
            values_(std::move(values)), errors_(std::move(errors)), prediction_(prediction), prior_(prior), next_(next),
            parameters_(std::move(start)), logLikelihood_(calculateLogLikelihood(parameters_)) {
    }

    std::pair<Parameters, double> next() {
        auto next = next_(parameters_, gen_);
        double nextLogLikelihood = calculateLogLikelihood(next);
        double accept = std::exp(nextLogLikelihood - logLikelihood_);

        std::uniform_real_distribution<> uniform(0, 1);
        if (uniform(gen_) < accept) {
            parameters_ = next;
            logLikelihood_ = nextLogLikelihood;
            return {next, nextLogLikelihood};
        }

        return {parameters_, nextLogLikelihood};
    }

private:
    double calculateLogLikelihood(const Parameters& parameters) {
        double logLikelihood = std::log(prior_(parameters));
        auto prediction = prediction_(parameters);

        for (size_t i = 0; i < prediction.size(); ++i) {
            logLikelihood -= std::pow(prediction[i] - values_[i], 2) / 2 / std::pow(errors_[i], 2);
        }
        return logLikelihood;
    }

private:
    std::vector<Value> values_;
    std::vector<Value> errors_;
    std::function<std::vector<Value>(const Parameters &)> prediction_;
    std::function<double(const Parameters &)> prior_;
    std::function<Parameters(const Parameters &, std::mt19937& gen)> next_;
    Parameters parameters_;
    double logLikelihood_;
    std::mt19937 gen_;
};
