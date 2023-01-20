/* -*- indent-tabs-mode: t -*- */

#ifndef SYSTEMS__ELECTRONS
#define SYSTEMS__ELECTRONS

#include <stdexcept>

#include <basis/field_set.hpp>
#include <basis/real_space.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <states/ks_states.hpp>
#include <hamiltonian/ks_hamiltonian.hpp>
#include <hamiltonian/energy.hpp>
#include <ions/brillouin.hpp>
#include <observables/density.hpp>
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
	
	static auto lot_subcomm(parallel::cartesian_communicator<3> & comm){
		return comm.axis(input::parallelization::dimension_kpoints());
	}
	static auto states_subcomm(parallel::cartesian_communicator<3> & comm){
		return comm.axis(input::parallelization::dimension_states());
	}
	static auto basis_subcomm(parallel::cartesian_communicator<3> & comm){
		return comm.axis(input::parallelization::dimension_domains());
	}
	static auto states_basis_subcomm(parallel::cartesian_communicator<3> & comm){
		return comm.plane(input::parallelization::dimension_domains(), input::parallelization::dimension_states());
	}
	static auto lot_states_subcomm(parallel::cartesian_communicator<3> & comm){
		return comm.plane(input::parallelization::dimension_kpoints(), input::parallelization::dimension_states());
	}
	
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



	electrons(input::parallelization const & dist, const inq::systems::ions & ions, systems::box const & box, input::kpoints const & kpts, const input::config & conf = {}):
		electrons(dist, ions, box, conf, kpts)
	{
	}
	
	electrons(input::parallelization const & dist, const inq::systems::ions & ions, systems::box const & box, const input::config & conf = {}, input::kpoints const & kpts = input::kpoints::gamma()):
		brillouin_zone_(ions, kpts),
		full_comm_(dist.cart_comm(conf.num_spin_components_val(), brillouin_zone_.size())),
		lot_comm_(lot_subcomm(full_comm_)),
		lot_states_comm_(lot_states_subcomm(full_comm_)),
		states_comm_(states_subcomm(full_comm_)),
		states_basis_comm_(states_basis_subcomm(full_comm_)),
		states_basis_(box, basis_subcomm(full_comm_)),
		density_basis_(states_basis_), /* disable the fine density mesh for now density_basis_(states_basis_.refine(arg_basis_input.density_factor(), basis_comm_)), */
		atomic_pot_(ions.geo().num_atoms(), ions.geo().atoms(), states_basis_.gcutoff()),
		states_(conf.spin_val(), atomic_pot_.num_electrons() + conf.excess_charge_val(), conf.extra_states_val(), conf.temperature_val(), kpts.num()),
		spin_density_(density_basis_, states_.num_density_components()),
		lot_part_(kpts.num()*states_.num_spin_indices(), lot_comm_)
	{
		CALI_CXX_MARK_FUNCTION;

		assert(lot_part_.local_size() > 0);
		assert(density_basis_.comm().size() == states_basis_.comm().size());

		auto nproc_spin = 1;
		if(states_.num_spin_indices() == 2 and lot_comm_.size()%2 == 0) nproc_spin = 2;

		parallel::cartesian_communicator<2> spin_kpoints_comm(lot_comm_, {nproc_spin, boost::mpi3::fill});

		parallel::partition spin_part(states_.num_spin_indices(), spin_kpoints_comm.axis(0));
		parallel::partition kpts_part(kpts.num(), spin_kpoints_comm.axis(1));

		assert(lot_part_.local_size() == kpts_part.local_size()*spin_part.local_size()); //this is always true because the spin size is either 1 or 2

		lot_weights_.reextent({lot_part_.local_size()});

		max_local_size_ = 0;
		auto ilot = 0;
		for(int ispin = 0; ispin < spin_part.local_size(); ispin++){
			for(int ikpt = 0; ikpt < kpts_part.local_size(); ikpt++){
				lot_weights_[ilot] = brillouin_zone_.kpoint_weight(kpts_part.local_to_global(ikpt).value());
				auto kpoint = brillouin_zone_.kpoint(kpts_part.local_to_global(ikpt).value());
				lot_.emplace_back(states_basis_, states_.num_states(), kpoint, spin_part.local_to_global(ispin).value(), states_basis_comm_);
				max_local_size_ = std::max(max_local_size_, lot_[ikpt].local_set_size());
				ilot++;
			}
		}
		
		assert(long(lot_.size()) == lot_part_.local_size());
		assert(max_local_size_ > 0);
		
		eigenvalues_.reextent({static_cast<boost::multi::size_t>(lot_.size()), max_local_size_});
		occupations_.reextent({static_cast<boost::multi::size_t>(lot_.size()), max_local_size_});

		if(atomic_pot_.num_electrons() + conf.excess_charge_val() == 0) throw std::runtime_error("inq error: the system does not have any electrons");
		
		print(ions);

	}

	electrons(electrons && old_el, input::parallelization const & new_dist):
		brillouin_zone_(std::move(old_el.brillouin_zone_)),
		full_comm_(new_dist.cart_comm(old_el.states_.num_spin_indices(), brillouin_zone_.size())),
		lot_comm_(lot_subcomm(full_comm_)),
		lot_states_comm_(lot_states_subcomm(full_comm_)),
		states_comm_(states_subcomm(full_comm_)),
		states_basis_comm_(states_basis_subcomm(full_comm_)),
		states_basis_(std::move(old_el.states_basis_), basis_subcomm(full_comm_)),
		density_basis_(std::move(old_el.density_basis_), basis_subcomm(full_comm_)),
		atomic_pot_(std::move(old_el.atomic_pot_)),
		states_(std::move(old_el.states_)),
		lot_weights_(std::move(old_el.lot_weights_)),
		max_local_size_(std::move(old_el.max_local_size_)),
		spin_density_(std::move(old_el.spin_density_), density_basis_.comm()),
		logger_(std::move(old_el.logger_)),
		lot_part_(std::move(old_el.lot_part_))
	{

		assert(lot_comm_ == old_el.lot_comm_); //resizing of k points not supported for the moment

		max_local_size_ = 0;
		for(auto & oldphi : old_el.lot_){
			lot_.emplace_back(std::move(oldphi), states_basis_comm_);
			max_local_size_ = std::max(max_local_size_, lot_.back().local_set_size());
		}

		assert(lot_.size() == old_el.lot_.size());

		eigenvalues_.reextent({static_cast<boost::multi::size_t>(lot_.size()), max_local_size_});
		occupations_.reextent({static_cast<boost::multi::size_t>(lot_.size()), max_local_size_});
		
		for(unsigned ilot = 0; ilot < lot_.size(); ilot++){

			parallel::partition part(lot_[ilot].set_size(), states_subcomm(old_el.full_comm_));
			
			parallel::array_iterator eigit(part, states_subcomm(old_el.full_comm_), +old_el.eigenvalues_[ilot]);
			parallel::array_iterator occit(part, states_subcomm(old_el.full_comm_), +old_el.occupations_[ilot]);

			for(; eigit != eigit.end(); ++eigit){
				
				for(int ist = 0; ist < eigenvalues_[ilot].size(); ist++){
					auto istg = lot_[ilot].set_part().local_to_global(ist);
					if(part.contains(istg.value())){
						eigenvalues_[ilot][ist] = (*eigit)[part.global_to_local(istg)];
						occupations_[ilot][ist] = (*occit)[part.global_to_local(istg)];
					}
				}
				
				++occit;
			}
			
		}		
		
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
			if(ions.geo().num_atoms() > 0) logger()->info("system symmetries: " + ions.symmetry_string());	
		}

		auto myid = gpu::id();
		auto gpuids = full_comm_.all_gather_as<boost::multi::array<decltype(myid), 1>>(myid);

		basis::fourier_space fourier_basis(states_basis_);
			
		if(logger()){
			logger()->info("parallelization:");
			logger()->info("  electrons divided among {} processes ({} kpoints x {} domains x {} states)", full_comm_.size(), full_comm_.shape()[0], full_comm_.shape()[1], full_comm_.shape()[2]);
#ifdef ENABLE_CUDA
			for(int iproc = 0; iproc < full_comm_.size(); iproc++){
				logger()->info("  process {} has gpu id {}", iproc, gpuids[iproc]);
			}
#else
			logger()->info("  inq is running on the cpu\n");
#endif
			logger()->info("k-point parallelization:");
			logger()->info("  {} k-points/spin indices divided among {} partitions", lot_part_.size(), lot_part_.comm_size());
			logger()->info("  partition 0 has {} k-points/spin indices and the last partition has {}\n", lot_part_.local_size(0), lot_part_.local_size(lot_part_.comm_size() - 1));
			
			logger()->info("real-space parallelization:");
			logger()->info("  {} slices ({} points) divided among {} partitions", states_basis_.cubic_dist(0).size(), states_basis_.part().size(), states_basis_.cubic_dist(0).comm_size());
			logger()->info("  partition 0 has {} slices and the last partition has {} slices ({} and {} points)", states_basis_.cubic_dist(0).local_size(0), states_basis_.cubic_dist(0).local_size(states_basis_.part().comm_size() - 1),
										 states_basis_.part().local_size(0), states_basis_.part().local_size(states_basis_.part().comm_size() - 1));

			logger()->info("fourier-space parallelization:");
			logger()->info("  {} slices ({} points) divided among {} partitions", fourier_basis.cubic_dist(2).size(), fourier_basis.part().size(), fourier_basis.cubic_dist(2).comm_size());
			logger()->info("  partition 0 has {} slices and the last partition has {} slices ({} and {} points)\n", fourier_basis.cubic_dist(2).local_size(0), fourier_basis.cubic_dist(2).local_size(fourier_basis.part().comm_size() - 1),
										 fourier_basis.part().local_size(0), fourier_basis.part().local_size(fourier_basis.part().comm_size() - 1));

			logger()->info("state parallelization:");
			logger()->info("  {} states divided among {} partitions", lot()[0].set_part().size(), lot()[0].set_part().comm_size());
			logger()->info("  partition 0 has {} states and the last partition has {} states\n", lot()[0].set_part().local_size(0), lot()[0].set_part().local_size(lot()[0].set_part().comm_size() - 1));
				
		}
	}

	template <typename ArrayType>
	void update_occupations(ArrayType const eigenval) {
		states_.update_occupations(lot_states_comm_, eigenval, occupations_);
	}

	void save(std::string const & dirname) const {
		int iphi = 0;
		for(auto & phi : lot()){
			auto basedir = dirname + "/lot" + operations::io::numstr(iphi + lot_part_.start());
			operations::io::save(basedir + "/states", phi);
			if(states_basis_.comm().root()) operations::io::save(basedir + "/occupations", states_comm_, lot()[iphi].set_part(), +occupations()[iphi]);	
			iphi++;
		}
	}
		
	auto load(std::string const & dirname) {
		auto success = true;

		int iphi = 0;
		for(auto & phi : lot()){
			auto basedir = dirname + "/lot" + operations::io::numstr(iphi + lot_part_.start());
			success = success and operations::io::load(basedir + "/states", phi.fields());

			math::array<double, 1> tmpocc(lot()[iphi].set_part().local_size());
			success = success and operations::io::load(basedir + "/occupations", states_comm_, lot()[iphi].set_part(), tmpocc);
			occupations()[iphi] = tmpocc;
			
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

	long lot_size() const {
		return lot().size();
	}

	auto & lot_weights() const {
		return lot_weights_;
	}

	auto & lot_part() const {
		return lot_part_;
	}
	
	auto max_local_size() const {
		return max_local_size_;
	}

	auto density() const {
		return observables::density::total(spin_density_);
	}

	auto & spin_density() const {
		return spin_density_;
	}

	auto & spin_density() {
		return spin_density_;
	}

	auto & states() const {
		return states_;
	}
	
private:
	static std::string generate_tiny_uuid(){
		auto uuid = boost::uuids::random_generator{}();
		uint32_t tiny = hash_value(uuid) % std::numeric_limits<uint32_t>::max();
		using namespace boost::archive::iterators;
		using it = base64_from_binary<transform_width<unsigned char*, 6, 8>>;
		return std::string(it((unsigned char*)&tiny), it((unsigned char*)&tiny+sizeof(tiny)));//.append((3-sizeof(tiny)%3)%3,'=');
	}

	inq::ions::brillouin brillouin_zone_;	
	
public: //temporary hack to be able to apply a kick from main and avoid a bug in nvcc

	mutable parallel::cartesian_communicator<3> full_comm_;
	mutable parallel::cartesian_communicator<1> lot_comm_;
	mutable parallel::cartesian_communicator<2> lot_states_comm_;
	mutable parallel::cartesian_communicator<1> states_comm_;
	mutable parallel::cartesian_communicator<2> states_basis_comm_;
	basis::real_space states_basis_;
	basis::real_space density_basis_;
	hamiltonian::atomic_potential atomic_pot_;
private:
	states::ks_states states_;
	std::vector<states::orbital_set<basis::real_space, complex>> lot_;
	math::array<double, 2> eigenvalues_;
	math::array<double, 2> occupations_;
	math::array<double, 1> lot_weights_;
	long max_local_size_;
	basis::field_set<basis::real_space, double> spin_density_;
 	
public:
	std::shared_ptr<spdlog::logger> const& logger() const{return logger_;}
private:
	std::shared_ptr<spdlog::logger> logger_;

	inq::parallel::partition lot_part_;

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

	systems::box box = systems::box::orthorhombic(6.0_b, 1.0_b, 6.0_b).cutoff_energy(15.0_Ha);
	
	systems::ions ions(box);

	ions.insert("Cu", {0.0_b,  0.0_b,  0.0_b});
	ions.insert("Cu", {1.0_b,  0.0_b,  0.0_b});
	
	auto par = input::parallelization(comm);
		
	systems::electrons electrons(par, ions, box);

	CHECK(electrons.states().num_electrons() == 38.0_a);
	CHECK(electrons.states().num_states() == 19);

	int iphi = 0;
	for(auto & phi : electrons.lot()) {
		
		for(int ist = 0; ist < phi.set_part().local_size(); ist++){
			auto istg = phi.set_part().local_to_global(ist);
			
			electrons.occupations()[iphi][ist] = cos(istg.value());
			
			for(int ip = 0; ip < phi.basis().local_size(); ip++){
				auto ipg = phi.basis().part().local_to_global(ip);
				phi.matrix()[ip][ist] = 20.0*(ipg.value() + 1)*sqrt(istg.value());
			}
		}
		iphi++;
	}

	electrons.save("electron_restart");
	
	systems::electrons electrons_read(par, ions, box);
	
	electrons_read.load("electron_restart");
	
	CHECK(electrons.lot_size() == electrons_read.lot_size());
	
	iphi = 0;
	for(auto & phi : electrons.lot()) {
		
		for(int ist = 0; ist < phi.set_part().local_size(); ist++){
			CHECK(electrons.occupations()[iphi][ist] == electrons_read.occupations()[0][ist]);
			for(int ip = 0; ip < phi.basis().local_size(); ip++){
				CHECK(phi.matrix()[ip][ist] == electrons_read.lot()[iphi].matrix()[ip][ist]);
			}
		}
		iphi++;
	}
	
	SECTION("Redistribute"){
		
		systems::electrons newel(std::move(electrons_read), input::parallelization(comm).domains(1).states());
		
		newel.save("newel_restart");
		
		systems::electrons newel_read(par, ions, box);
		newel_read.load("newel_restart");
		
		CHECK(electrons.lot_size() == newel_read.lot_size());
		
		iphi = 0;
		for(auto & phi : electrons.lot()) {
			
			for(int ist = 0; ist < phi.set_part().local_size(); ist++){
				CHECK(electrons.occupations()[iphi][ist] == newel_read.occupations()[0][ist]);
				for(int ip = 0; ip < phi.basis().local_size(); ip++){
					CHECK(phi.matrix()[ip][ist] == newel_read.lot()[iphi].matrix()[ip][ist]);
				}
			}
			iphi++;
		}

	}

	CHECK(not electrons.load("directory_that_doesnt_exist"));
	
}

#endif

#endif

