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
#include <input/config.hpp>
#include <input/interaction.hpp>
#include <input/kpoints.hpp>
#include <input/rt.hpp>
#include <input/scf.hpp>
#include <ions/interaction.hpp>
#include <systems/ions.hpp>
#include <states/orbital_set.hpp>

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

	auto & lot() const {
		return lot_;
	}

	auto & lot() {
		return lot_;
	}

	auto & occupations() const {
		return occupations_;
	}

	auto & occupations() {
		return occupations_;
	}
		
	electrons(input::distribution const & dist, const inq::systems::ions & ions, systems::box const & box, const input::config & conf = {}, input::kpoints const & kpts = input::kpoints::gamma()):
		full_comm_(dist.cart_comm()),
		lot_comm_(full_comm_.axis(0)),
		lot_states_comm_(full_comm_.hyperplane(2)),
		states_comm_(full_comm_.axis(1)),
		atoms_comm_(states_comm_),
		states_basis_comm_(full_comm_.hyperplane(0)),
		basis_comm_(full_comm_.axis(2)),
		states_basis_(box, basis_comm_),
		density_basis_(states_basis_), /* disable the fine density mesh for now density_basis_(states_basis_.refine(arg_basis_input.density_factor(), basis_comm_)), */
		atomic_pot_(ions.geo().num_atoms(), ions.geo().atoms(), states_basis_.gcutoff(), atoms_comm_),
		states_(states::ks_states::spin_config::UNPOLARIZED, atomic_pot_.num_electrons() + conf.excess_charge, conf.extra_states, conf.temperature.in_atomic_units(), kpts.num()),
		density_(density_basis_),
		lot_part_(kpts.num(), lot_comm_)
	{

		CALI_CXX_MARK_FUNCTION;

		assert(density_basis_.comm().size() == states_basis_.comm().size());

		std::vector<int> grid_address(3*kpts.num());
		std::vector<int> map(kpts.num());
		std::vector<int> types(ions.geo().num_atoms());
		std::vector<double> positions(3*ions.geo().num_atoms());
		
		for(int iatom = 0; iatom < ions.geo().num_atoms(); iatom++){
			types[iatom] = ions.geo().atoms()[iatom].atomic_number();
			auto pos = ions.cell().cart_to_crystal(ions.cell().position_in_cell(ions.geo().coordinates()[iatom]));
			positions[3*iatom + 0] = pos[0];
			positions[3*iatom + 1] = pos[1];
			positions[3*iatom + 2] = pos[2];
		}

		auto is_shifted = kpts.is_shifted();
		
		spg_get_ir_reciprocal_mesh(reinterpret_cast<int (*)[3]>(grid_address.data()), map.data(), (int const *) &kpts.dims(), (int const *) &is_shifted, 0,
															 reinterpret_cast<double (*)[3]>(const_cast<double *>(ions.cell().amat())),
															 reinterpret_cast<double (*)[3]>(positions.data()), types.data(), ions.geo().num_atoms(), 1e-4);

		lot_weights_.reextent({lot_part_.local_size()});

		max_local_size_ = 0;
		for(int ikpt = 0; ikpt < kpts.num(); ikpt++){
			if(not lot_part_.contains(ikpt)) continue;

			math::vector3<double> kpoint = {(grid_address[3*ikpt + 0] + 0.5*is_shifted[0])/kpts.dims()[0], (grid_address[3*ikpt + 1] + 0.5*is_shifted[1])/kpts.dims()[1], (grid_address[3*ikpt + 2] + 0.5*is_shifted[2])/kpts.dims()[2]};
			kpoint = 2.0*M_PI*ions.cell().cart_to_crystal(kpoint);
			lot_.emplace_back(states_basis_, states_.num_states(), kpoint, states_basis_comm_);

			auto likpt = lot_part_.global_to_local(utils::global_index(ikpt));
			lot_weights_[likpt] = 1.0/kpts.num();
			max_local_size_ = std::max(max_local_size_, lot_[likpt].fields().local_set_size());
		}

		assert(long(lot_.size()) == lot_part_.local_size());
		assert(max_local_size_ > 0);
		
		eigenvalues_.reextent({lot_.size(), max_local_size_});
		occupations_.reextent({lot_.size(), max_local_size_});

		if(atomic_pot_.num_electrons() + conf.excess_charge == 0) throw error::NO_ELECTRONS;
		
		print(ions);

	}

	void print(const inq::systems::ions & ions){
		if(full_comm_.root()){
			logger_ = spdlog::stdout_color_mt("electrons:"+ generate_tiny_uuid());
			logger_->set_level(spdlog::level::trace);
		}

		if(logger()){
			logger()->info("constructed with basis {}", states_basis_);
			logger()->info("constructed with states {}", states_);
		}
			

		if(logger()){
			logger()->info("constructed with geometry {}", ions.geo_);
			logger()->info("constructed with cell {}", ions.cell_);
			logger()->info("system symmetries: " + ions.symmetry_string());	
		}

		auto myid = gpu::id();
		auto gpuids = full_comm_.all_gather_as<boost::multi::array<decltype(myid), 1>>(myid);

		basis::fourier_space fourier_basis(states_basis_);
			
		if(logger()){
			logger()->info("parallelization:");
			logger()->info("  electrons divided among {} processes ({} states x {} domains)", full_comm_.size(), full_comm_.shape()[0], full_comm_.shape()[1]);
#ifdef ENABLE_CUDA
			for(int iproc = 0; iproc < full_comm_.size(); iproc++){
				logger()->info("  process {} has gpu id {}", iproc, gpuids[iproc]);
			}
#else
			logger()->info("  inq is running on the cpu\n");
#endif
			logger()->info("k-point parallelization:");
			logger()->info("  {} k-points divided among {} partitions", lot_part_.size(), lot_part_.comm_size());
			logger()->info("  partition 0 has {} k-points and the last partition has {} k-points\n", lot_part_.local_size(0), lot_part_.local_size(lot_part_.comm_size() - 1));
			
			logger()->info("state parallelization:");
			logger()->info("  {} states divided among {} partitions", lot()[0].fields().set_part().size(), lot()[0].fields().set_part().comm_size());
			logger()->info("  partition 0 has {} states and the last partition has {} states\n", lot()[0].fields().set_part().local_size(0), lot()[0].fields().set_part().local_size(lot()[0].fields().set_part().comm_size() - 1));
				
			logger()->info("real-space parallelization:");
			logger()->info("  {} slices ({} points) divided among {} partitions", states_basis_.cubic_dist(0).size(), states_basis_.part().size(), states_basis_.cubic_dist(0).comm_size());
			logger()->info("  partition 0 has {} slices and the last partition has {} slices ({} and {} points)\n", states_basis_.cubic_dist(0).local_size(0), states_basis_.cubic_dist(0).local_size(states_basis_.part().comm_size() - 1),
										 states_basis_.part().local_size(0), states_basis_.part().local_size(states_basis_.part().comm_size() - 1));

			logger()->info("fourier-space parallelization:");
			logger()->info("  {} slices ({} points) divided among {} partitions", fourier_basis.cubic_dist(2).size(), fourier_basis.part().size(), fourier_basis.cubic_dist(2).comm_size());
			logger()->info("  partition 0 has {} slices and the last partition has {} slices ({} and {} points)\n", fourier_basis.cubic_dist(2).local_size(0), fourier_basis.cubic_dist(2).local_size(fourier_basis.part().comm_size() - 1),
										 fourier_basis.part().local_size(0), fourier_basis.part().local_size(fourier_basis.part().comm_size() - 1));
				
		}		
	}

	template <typename ArrayType>
	void update_occupations(ArrayType const eigenval) {
		states_.update_occupations(lot_states_comm_, eigenval, occupations_);
	}

	void save(std::string const & dirname) const {
		int iphi = 0;
		for(auto & phi : lot()){
			operations::io::save(dirname + "/lot" + operations::io::numstr(iphi + lot_part_.start()) + "/states", phi);
			iphi++;
		}
		if(basis_comm_.root()) operations::io::save(dirname + "/ocupations", lot_states_comm_, occupations().size()*lot()[0].fields().set_part(), occupations());
	}
		
	auto load(std::string const & dirname) {
		auto success = operations::io::load(dirname + "/ocupations", lot_states_comm_, occupations().size()*lot()[0].fields().set_part(), occupations());

		int iphi = 0;
		for(auto & phi : lot()){
			success = success and operations::io::load(dirname + "/lot" + operations::io::numstr(iphi + lot_part_.start()) + "/states", phi.fields());
			iphi++;
		}
		return success;
	}

	auto & eigenvalues() const {
		return eigenvalues_;
	}

	auto & eigenvalues() {
		return eigenvalues_;
	}

	long lot_size() {
		return lot().size();
	}

	auto & lot_weights() const {
		return lot_weights_;
	}

	auto max_local_size() const {
		return max_local_size_;
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

	mutable boost::mpi3::cartesian_communicator<3> full_comm_;
	mutable boost::mpi3::cartesian_communicator<1> lot_comm_;
	mutable boost::mpi3::cartesian_communicator<2> lot_states_comm_;
	mutable boost::mpi3::cartesian_communicator<1> states_comm_;
	mutable boost::mpi3::cartesian_communicator<1> atoms_comm_;
	mutable boost::mpi3::cartesian_communicator<2> states_basis_comm_;	
	mutable boost::mpi3::cartesian_communicator<1> basis_comm_;
	basis::real_space states_basis_;
	basis::real_space density_basis_;
	hamiltonian::atomic_potential atomic_pot_;
	states::ks_states states_;
private:
	std::vector<states::orbital_set<basis::real_space, complex>> lot_;
	math::array<double, 2> eigenvalues_;
	math::array<double, 2> occupations_;
	math::array<double, 1> lot_weights_;
	long max_local_size_;

 	
public:
	basis::field<basis::real_space, double> density_;
	std::shared_ptr<spdlog::logger> const& logger() const{return logger_;}
private:
	std::shared_ptr<spdlog::logger> logger_;

	inq::utils::partition lot_part_;

};

}
}

#ifdef INQ_SYSTEMS_ELECTRONS_UNIT_TEST
#undef INQ_SYSTEMS_ELECTRONS_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE("class system::electrons", "[system::electrons]") {
	
	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	
	auto comm = boost::mpi3::environment::get_world_instance();

	std::vector<input::atom> geo;
	geo.push_back( "Cu" | math::vector3<double>(0.0,  0.0,  0.0));
	geo.push_back( "Cu" | math::vector3<double>(1.0,  0.0,  0.0));

	systems::box box = systems::box::cubic(5.0_b).cutoff_energy(25.0_Ha);
	
	systems::ions ions(box, geo);

	auto dist = input::distribution(comm);
		
	systems::electrons electrons(dist, ions, box);

	CHECK(electrons.states_.num_electrons() == 38.0_a);
	CHECK(electrons.states_.num_states() == 19);

	for(auto & phi : electrons.lot()) {
		
		for(int ist = 0; ist < phi.fields().set_part().local_size(); ist++){
			auto istg = phi.fields().set_part().local_to_global(ist);
			
			electrons.occupations()[0][ist] = cos(istg.value());
			
			for(int ip = 0; ip < phi.fields().basis().local_size(); ip++){
				auto ipg = phi.fields().basis().part().local_to_global(ip);
				phi.fields().matrix()[ip][ist] = 20.0*(ipg.value() + 1)*sqrt(istg.value());
			}
		}
	}
		
	electrons.save("electron_restart");
	
	systems::electrons electrons_read(dist, ions, box);
	
	electrons_read.load("electron_restart");

	CHECK(electrons.lot_size() == electrons_read.lot_size());
	
	int iphi = 0;
	for(auto & phi : electrons.lot()) {

		for(int ist = 0; ist < phi.fields().set_part().local_size(); ist++){
			CHECK(electrons.occupations()[0][ist] == electrons_read.occupations()[0][ist]);
			for(int ip = 0; ip < phi.fields().basis().local_size(); ip++){
				CHECK(phi.fields().matrix()[ip][ist] == electrons_read.lot()[iphi].fields().matrix()[ip][ist]);
			}
		}
		iphi++;
	}
	
	CHECK(not electrons.load("directory_that_doesnt_exist"));

}

#endif

#endif

