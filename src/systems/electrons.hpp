/* -*- indent-tabs-mode: t -*- */

#ifndef SYSTEMS__ELECTRONS
#define SYSTEMS__ELECTRONS

#include <basis/real_space.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <states/ks_states.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <hamiltonian/energy.hpp>
#include <basis/field_set.hpp>
#include <operations/randomize.hpp>
#include <operations/integral.hpp>
#include <operations/io.hpp>
#include <operations/orthogonalize.hpp>
#include <math/complex.hpp>
#include <input/basis.hpp>
#include <input/config.hpp>
#include <input/interaction.hpp>
#include <input/rt.hpp>
#include <input/scf.hpp>
#include <ions/interaction.hpp>
#include <systems/ions.hpp>
#include <density/calculate.hpp>
#include <density/normalize.hpp>

#include <boost/uuid/uuid.hpp>
#include <boost/uuid/uuid_generators.hpp> // uuids::random_generator

#include <cfloat>

#include <boost/archive/iterators/base64_from_binary.hpp>
#include <boost/archive/iterators/transform_width.hpp>

#include<spdlog/spdlog.h>
#include<spdlog/sinks/stdout_color_sinks.h>
#include<spdlog/fmt/ostr.h> // print user defined types

#include <utils/profiling.hpp>

namespace inq {
namespace systems {

	class electrons {
		
	public:
	
		enum class error { NO_ELECTRONS };

		electrons(boost::mpi3::cartesian_communicator<2> cart_comm, const inq::systems::ions & ions, const input::basis arg_basis_input, const input::config & conf = {}):
			full_comm_(cart_comm),
			states_comm_(full_comm_.axis(0)),
			atoms_comm_(states_comm_),
			basis_comm_(full_comm_.axis(1)),
			states_basis_(ions.cell(), arg_basis_input, basis_comm_),
			density_basis_(states_basis_), /* disable the fine density mesh for now density_basis_(states_basis_.refine(arg_basis_input.density_factor(), basis_comm_)), */
			atomic_pot_(ions.geo().num_atoms(), ions.geo().atoms(), states_basis_.gcutoff(), atoms_comm_),
			states_(states::ks_states::spin_config::UNPOLARIZED, atomic_pot_.num_electrons() + conf.excess_charge, conf.extra_states, conf.temperature.in_atomic_units()),
			phi_(states_basis_, states_.num_states(), full_comm_),
			density_(density_basis_),
			occupations_(phi_.local_set_size())
		{

			CALI_CXX_MARK_FUNCTION;
			
			assert(density_basis_.comm().size() == states_basis_.comm().size());

			if(full_comm_.root()){
				logger_ = spdlog::stdout_color_mt("electrons:"+ generate_tiny_uuid());
				logger_->set_level(spdlog::level::trace);
			}

			if(logger()){
				logger()->info("constructed with basis {}", states_basis_);
				logger()->info("constructed with states {}", states_);
			}
			
			if(atomic_pot_.num_electrons() + conf.excess_charge == 0) throw error::NO_ELECTRONS;

			if(logger()){
				logger()->info("constructed with geometry {}", ions.geo_);
				logger()->info("constructed with cell {}", ions.cell_);
			}

			auto myid = gpu::id();
			auto gpuids = full_comm_.all_gather_as<boost::multi::array<decltype(myid), 1>>(myid);

			if(logger()){
				logger()->info("electrons divided among {} processes ({} states x {} domains)", full_comm_.size(), full_comm_.shape()[0], full_comm_.shape()[1]);
#ifdef ENABLE_CUDA
				for(int iproc = 0; iproc < full_comm_.size(); iproc++){
					logger()->info("  process {} has gpu id {}", iproc, gpuids[iproc]);
				}
#else
				logger()->info("  inq is running on the cpu");
#endif
			}
			
		}

		electrons(boost::mpi3::communicator & comm, const inq::systems::ions & ions, const input::basis arg_basis_input, const input::config & conf = {}):
			electrons(boost::mpi3::cartesian_communicator<2>{comm, {1, boost::mpi3::fill}}, ions, arg_basis_input, conf){
		}

		template <typename ArrayType>
		void update_occupations(ArrayType const eigenval) {
			states_.update_occupations(phi_.set_comm(), phi_.set_part(), eigenval, occupations_);
		}

		void save(std::string const & dirname) const {
			operations::io::save(dirname + "/states", phi_);
			if(phi_.basis().comm().root()) operations::io::save(dirname + "/ocupations", phi_.set_comm(), phi_.set_part(), occupations_);
		}
		
		auto load(std::string const & dirname) {
			return operations::io::load(dirname + "/states", phi_)
				and operations::io::load(dirname + "/ocupations", phi_.set_comm(), phi_.set_part(), occupations_);
		}
		
	private:
		static std::string generate_tiny_uuid(){
			auto uuid = boost::uuids::random_generator{}();
			uint32_t tiny = hash_value(uuid) % std::numeric_limits<uint32_t>::max();
			using namespace boost::archive::iterators;
			using it = base64_from_binary<transform_width<unsigned char*, 6, 8>>;
			return std::string(it((unsigned char*)&tiny), it((unsigned char*)&tiny+sizeof(tiny)));//.append((3-sizeof(tiny)%3)%3,'=');
		}

	public: //temporary hack to be able to apply a kick from main and avoid a bug in nvcc

		mutable boost::mpi3::cartesian_communicator<2> full_comm_;
		mutable boost::mpi3::cartesian_communicator<1> states_comm_;
		mutable boost::mpi3::cartesian_communicator<1> atoms_comm_;
		mutable boost::mpi3::cartesian_communicator<1> basis_comm_;
		basis::real_space states_basis_;
		basis::real_space density_basis_;
		hamiltonian::atomic_potential atomic_pot_;
		states::ks_states states_;
		basis::field_set<basis::real_space, complex> phi_;
		basis::field<basis::real_space, double> density_;
		math::array<double, 1> occupations_;

		std::shared_ptr<spdlog::logger> const& logger() const{return logger_;}
	private:
		std::shared_ptr<spdlog::logger> logger_;
	};

}
}

#ifdef INQ_SYSTEMS_ELECTRONS_UNIT_TEST
#undef INQ_SYSTEMS_ELECTRONS_UNIT_TEST

#include <catch2/catch.hpp>

TEST_CASE("class system::electrons", "[system::electrons]") {
	
	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	
	auto comm = boost::mpi3::environment::get_world_instance();
	
	boost::mpi3::cartesian_communicator<2> cart_comm(comm, {});

	std::vector<input::atom> geo;
	geo.push_back( "Si" | math::vector3<double>(0.0,  0.0,  0.0));

	systems::ions ions(input::cell::cubic(5.0_b), geo);

	systems::electrons electrons(cart_comm, ions, input::basis::cutoff_energy(25.0_Ha));
	
	for(int ist = 0; ist < electrons.phi_.set_part().local_size(); ist++){
		auto istg = electrons.phi_.set_part().local_to_global(ist);

		electrons.occupations_[ist] = cos(istg.value());
		
		for(int ip = 0; ip < electrons.phi_.basis().local_size(); ip++){
			auto ipg = electrons.phi_.basis().part().local_to_global(ip);
			electrons.phi_.matrix()[ip][ist] = 20.0*(ipg.value() + 1)*sqrt(istg.value());
		}
	}

	electrons.save("electron_restart");
	
	systems::electrons electrons_read(cart_comm, ions, input::basis::cutoff_energy(25.0_Ha));

	electrons_read.load("electron_restart");

	for(int ist = 0; ist < electrons.phi_.set_part().local_size(); ist++){
		CHECK(electrons.occupations_[ist] == electrons_read.occupations_[ist]);
		for(int ip = 0; ip < electrons.phi_.basis().local_size(); ip++){
			CHECK(electrons.phi_.matrix()[ip][ist] == electrons_read.phi_.matrix()[ip][ist]);
		}
	}

	CHECK(not electrons.load("directory_that_doesnt_exist"));

}

#endif

#endif

