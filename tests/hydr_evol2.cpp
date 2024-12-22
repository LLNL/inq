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
    auto ev = b;

    auto ions = systems::ions(systems::cell::cubic(10.0_b));
    ions.insert("H", {0.0_b, 0.0_b, 0.0_b});

    auto electrons = systems::electrons(env.par(), ions, options::electrons{}.cutoff(30.0_Ha).extra_states(10).spin_non_collinear());
    ground_state::initial_guess(ions, electrons);

    auto result = ground_state::calculate(ions, electrons, options::theory{}.lda(), inq::options::ground_state{}.steepest_descent().energy_tolerance(1.e-8_Ha).max_steps(1000).mixing(0.1), B);
    auto mag = observables::total_magnetization(electrons.spin_density());
    std::cout << mag << std::endl;
    
    std::vector<double> mz;
    std::vector<double> mx;
    std::vector<double> my;
    std::vector<double> energy;
    std::vector<double> ion;
    std::vector<double> ion_kin;
    std::vector<double> bx;
    std::vector<double> by;
    std::vector<double> bz;
    std::vector<double> time;
    auto output = [&](auto data){
        energy.push_back(data.energy().total());
        ion.push_back(data.energy().ion());
        ion_kin.push_back(data.energy().ion_kinetic());
        mz.push_back(data.magnetization()[2]);
        mx.push_back(data.magnetization()[0]);
        my.push_back(data.magnetization()[1]);
        bx.push_back(data.uniform_magnetic_field()[0]);
        by.push_back(data.uniform_magnetic_field()[1]);
        bz.push_back(data.uniform_magnetic_field()[2]);
        time.push_back(data.time());
    };
    
    vector3<double> Db = {6.0, 0.0, -6.0};
    vector3<double> b2 = {1.0, 0.0, -1.0};
    perturbations::magnetic_pulse B2{b, Db, b2, 10.0_ev, 1000.0_attosecond, 1.0_invfs2, 100.0_attosecond, 0.0};
    real_time::propagate<>(ions, electrons, output, options::theory{}.lda(), options::real_time{}.num_steps(6000).dt(0.03_atomictime), B2);

    if (electrons.kpin_states_comm().rank() == 0) {
        std::ofstream outfile("real-time-mag.txt"); // Create an ofstream object and open the file
        if (outfile.is_open()) {
            for (auto i=0; i<time.size(); i++){
                outfile << time[i] << "            " << mx[i] << "            " << my[i] << "            " << mz[i] << std::endl;
            }
            outfile.close(); // Close the file
            std::cout << "File write successful." << std::endl;
        } else {
            std::cerr << "Error opening file." << std::endl;
        }
        std::ofstream outfile2("real-time-bfield.txt");
        if (outfile2.is_open()) {
            for (auto i=0; i<time.size(); i++){
                outfile2 << time[i] << "            " << bx[i] << "            " << by[i] << "            " << bz[i] << std::endl;
            }
        }
    }
}