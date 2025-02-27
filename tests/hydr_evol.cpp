/* -*- indent-tabs-mode: t -*- */

#include <inq/inq.hpp>
#include <perturbations/magnetic.hpp>
#include <perturbations/magnetic_pulse.hpp>

using namespace inq;
using namespace inq::magnitude;

inq::utils::match match(3.0e-4);

int main (int argc, char ** argv){
    auto env = input::environment{};
    vector3<double> b = {0.0, 0.0, 1.0};
    b = b / b.length();
    perturbations::magnetic B{b};

    auto ions = systems::ions(systems::cell::cubic(10.0_b));
    ions.insert("H", {0.0_b, 0.0_b, 0.0_b});

    {
        auto electrons = systems::electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).extra_states(10).spin_non_collinear());
        ground_state::initial_guess(ions, electrons);

        auto result = ground_state::calculate(ions, electrons, options::theory{}.lda(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1.e-8_Ha).max_steps(1000).mixing(0.1), B);
        auto mag = observables::total_magnetization(electrons.spin_density());
        match.check("magnetization direction",        mag[2],         b[2]);
        
        std::vector<double> mz;
        std::vector<double> mx;
        std::vector<double> my;
        auto output = [&](auto data){
            mz.push_back(data.magnetization()[2]);
            mx.push_back(data.magnetization()[0]);
            my.push_back(data.magnetization()[1]);
        };

        vector3<double> b2 = {3.0, 0.0, 0.0};
        perturbations::magnetic B2{b2};
        real_time::propagate<>(ions, electrons, output, options::theory{}.lda(), options::real_time{}.num_steps(200).dt(0.03_atomictime), B2);

        match.check("magnetization module: ", (mx[199]*mx[199] + my[199]*my[199] + mz[199]*mz[199]), mag[2]*mag[2]);
    }

    if (match.fail()) return match.fail();
    
    {
        auto electrons = systems::electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).extra_states(10).spin_non_collinear());
        ground_state::initial_guess(ions, electrons);

        auto result = ground_state::calculate(ions, electrons, options::theory{}.lda(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1.e-8_Ha).max_steps(1000).mixing(0.1), B);
        auto mag = observables::total_magnetization(electrons.spin_density());
        match.check("magnetization direction",      mag[2],        b[2]);
        
        std::vector<double> mz;
        std::vector<double> mx;
        std::vector<double> my;
        auto output = [&](auto data){
            mz.push_back(data.magnetization()[2]);
            mx.push_back(data.magnetization()[0]);
            my.push_back(data.magnetization()[1]);
        };
        
        vector3<double> Db = {6.0, 0.0, -6.0};
        vector3<double> b2 = {1.0, 0.0, -1.0};
        perturbations::magnetic_pulse B2{b, Db, b2, 10.0_ev, 1000.0_attosecond, 1.0_invfs2, 100.0_attosecond, 0.0};
        real_time::propagate<>(ions, electrons, output, options::theory{}.lda(), options::real_time{}.num_steps(200).dt(0.03_atomictime), B2);
        
        match.check("magnetization module: ", (mx[199]*mx[199] + my[199]*my[199] + mz[199]*mz[199]), mag[2]*mag[2]);
    }
    
    if (match.fail()) return match.fail();
}