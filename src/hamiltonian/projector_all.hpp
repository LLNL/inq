/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__HAMILTONIAN__PROJECTOR_ALL
#define INQ__HAMILTONIAN__PROJECTOR_ALL

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <gpu/array.hpp>
#include <math/vector3.hpp>
#include <basis/real_space.hpp>
#include <basis/spherical_grid.hpp>
#include <hamiltonian/atomic_potential.hpp>
#include <utils/profiling.hpp>
#include <utils/raw_pointer_cast.hpp>

namespace inq {
namespace hamiltonian {

class projector_all {

public: // for CUDA
	
	template <typename ProjectorsType>
	void constructor(ProjectorsType const & projectors){

		CALI_CXX_MARK_FUNCTION;

		max_sphere_size_ = 0;
		max_nlm_ = 0;
		for(auto it = projectors.cbegin(); it != projectors.cend(); ++it) {
			max_sphere_size_ = std::max(max_sphere_size_, it->sphere_.size());
			max_nlm_ = std::max(max_nlm_, it->nproj_);			
		}

		points_ = decltype(points_)({nprojs_, max_sphere_size_});
		positions_ = decltype(positions_)({nprojs_, max_sphere_size_});		
    coeff_ = decltype(coeff_)({nprojs_, max_nlm_}, 0.0);
    matrices_ = decltype(matrices_)({nprojs_, max_nlm_, max_sphere_size_});
		
    auto iproj = 0;
    for(auto it = projectors.cbegin(); it != projectors.cend(); ++it) {
			gpu::run(max_sphere_size_,
							 [poi = begin(points_), pos = begin(positions_), sph = it->sphere_.ref(), iproj, npoint = it->sphere_.size()] GPU_LAMBDA (auto ipoint){
								 if(ipoint < unsigned (npoint)){
									 poi[iproj][ipoint] = sph.grid_point(ipoint);	
									 pos[iproj][ipoint] = sph.point_pos(ipoint);							 
								 } else {
									 poi[iproj][ipoint] = {-1, -1, -1};
								 }
							 });
			
			gpu::run(max_sphere_size_, max_nlm_,
							 [mat = begin(matrices_), itmat = begin(it->matrix_), iproj, np = it->sphere_.size(), nlm = it->nproj_] GPU_LAMBDA (auto ipoint, auto ilm){
								 if(ipoint < (unsigned) np and ilm < (unsigned) nlm) {
									 mat[iproj][ilm][ipoint] = itmat[ilm][ipoint];
								 } else {
									 mat[iproj][ilm][ipoint] = 0.0;								 
								 }
							 });
								 
      coeff_[iproj]({0, it->nproj_}) = it->kb_coeff_;

			nlm_[iproj] = it->nproj_;
			iatom_[iproj] = it->iatom_;
			locally_empty_[iproj] = it->locally_empty();
			
      iproj++;
    }
    
	}
	
public:

	projector_all():
		nprojs_(0),
    max_sphere_size_(0),
    max_nlm_(0) {
  }
  
	////////////////////////////////////////////////////////////////////////////////////////////
	
	template <typename ProjectorsType>
	projector_all(ProjectorsType const & projectors):
		nprojs_(projectors.size()),
		nlm_(nprojs_),
		iatom_(nprojs_),
		locally_empty_(nprojs_)
	{
		constructor(projectors);
	}

	////////////////////////////////////////////////////////////////////////////////////////////		
	
	template <typename KpointType>
	gpu::array<complex, 3> project(states::orbital_set<basis::real_space, complex> const & phi, KpointType const & kpoint) const {
    
		gpu::array<complex, 3> sphere_phi_all({nprojs_, max_sphere_size_, phi.local_set_size()});
		gpu::array<complex, 3> projections_all({nprojs_, max_nlm_, phi.local_set_size()}, 0.0);

		{ CALI_CXX_MARK_SCOPE("projector::gather");
				
			gpu::run(phi.local_set_size(), max_sphere_size_, nprojs_,
							 [sgr = begin(sphere_phi_all), gr = begin(phi.hypercubic()), poi = begin(points_), pos = begin(positions_), kpoint] GPU_LAMBDA (auto ist, auto ipoint, auto iproj){
								 if(poi[iproj][ipoint][0] >= 0){
									 auto phase = polar(1.0, dot(kpoint, pos[iproj][ipoint]));
									 sgr[iproj][ipoint][ist] = phase*gr[poi[iproj][ipoint][0]][poi[iproj][ipoint][1]][poi[iproj][ipoint][2]][ist];
								 } else {
									 sgr[iproj][ipoint][ist] = complex(0.0, 0.0);
								 }
							 });
		}

#ifndef ENABLE_CUDA
		for(auto iproj = 0; iproj < nprojs_; iproj++){
			CALI_CXX_MARK_SCOPE("projector_gemm_1");

			if(locally_empty_[iproj]) continue;
			
			namespace blas = boost::multi::blas;
			blas::real_doubled(projections_all[iproj]) = blas::gemm(phi.basis().volume_element(), matrices_[iproj], blas::real_doubled(sphere_phi_all[iproj]));
		}
#else
		if(max_sphere_size_ > 0) {
			CALI_CXX_MARK_SCOPE("projector_gemm_1");			

			const double zero = 0.0;
			const double vol = phi.basis().volume_element();

			auto status = cublasDgemmStridedBatched(/*cublasHandle_t handle = */ boost::multi::cuda::cublas::context::get_instance().get(),
																							/*cublasOperation_t transa = */ CUBLAS_OP_N,
																							/*cublasOperation_t transb = */ CUBLAS_OP_N,
																							/*int m = */ 2*phi.local_set_size(),
																							/*int n = */ max_nlm_,
																							/*int k = */ max_sphere_size_,
																							/*const double *alpha = */ &vol,
																							/*const double *A = */ reinterpret_cast<double const *>(raw_pointer_cast(sphere_phi_all.data_elements())),
																							/*int lda = */ 2*phi.local_set_size(),
																							/*long long int strideA = */ 2*max_sphere_size_*phi.local_set_size(),
																							/*const double *B = */ raw_pointer_cast(matrices_.data_elements()),
																							/*int ldb = */ max_sphere_size_,
																							/*long long int strideB =*/ max_nlm_*max_sphere_size_,
																							/*const double *beta = */ &zero,
																							/*double *C = */ reinterpret_cast<double *>(raw_pointer_cast(projections_all.data_elements())),
																							/*int ldc = */ 2*phi.local_set_size(),
																							/*long long int strideC = */ 2*max_nlm_*phi.local_set_size(),
																							/*int batchCount = */ nprojs_);
			gpu::sync();
			
			assert(status == CUBLAS_STATUS_SUCCESS);
			
		}
#endif

    { CALI_CXX_MARK_SCOPE("projector_scal");
				
      gpu::run(phi.local_set_size(), max_nlm_, nprojs_,
               [proj = begin(projections_all), coe = begin(coeff_)]
               GPU_LAMBDA (auto ist, auto ilm, auto iproj){
                 proj[iproj][ilm][ist] = proj[iproj][ilm][ist]*coe[iproj][ilm];
               });
		}

		if(phi.basis().comm().size() > 1) {
			CALI_CXX_MARK_SCOPE("projector_all::project::reduce");
			phi.basis().comm().all_reduce_in_place_n(raw_pointer_cast(projections_all.data_elements()), projections_all.num_elements(), std::plus<>{});
		}

#ifndef ENABLE_CUDA
		for(auto iproj = 0; iproj < nprojs_; iproj++) {
			CALI_CXX_MARK_SCOPE("projector_gemm_2");

			if(locally_empty_[iproj]) continue;
			
			namespace blas = boost::multi::blas;
			blas::real_doubled(sphere_phi_all[iproj]) = blas::gemm(1., blas::T(matrices_[iproj]), blas::real_doubled(projections_all[iproj]));
		}
#else
		if(max_sphere_size_ > 0) {
			CALI_CXX_MARK_SCOPE("projector_gemm_2");

			const double zero = 0.0;
			const double one = 1.0;
			
			auto status = cublasDgemmStridedBatched(/*cublasHandle_t handle = */ boost::multi::cuda::cublas::context::get_instance().get(),
																							/*cublasOperation_t transa = */ CUBLAS_OP_N,
																							/*cublasOperation_t transb = */ CUBLAS_OP_T,
																							/*int m = */ 2*phi.local_set_size(),
																							/*int n = */ max_sphere_size_,
																							/*int k = */ max_nlm_,
																							/*const double *alpha = */ &one,
																							/*const double *A = */ reinterpret_cast<double const *>(raw_pointer_cast(projections_all.data_elements())),
																							/*int lda = */ 2*phi.local_set_size(),
																							/*long long int strideA = */ 2*max_nlm_*phi.local_set_size(),
																							/*const double *B = */ raw_pointer_cast(matrices_.data_elements()),
																							/*int ldb = */ max_sphere_size_,
																							/*long long int strideB =*/ max_nlm_*max_sphere_size_,
																							/*const double *beta = */ &zero,
																							/*double *C = */ reinterpret_cast<double *>(raw_pointer_cast(sphere_phi_all.data_elements())),
																							/*int ldc = */ 2*phi.local_set_size(),
																							/*long long int strideC = */ 2*max_sphere_size_*phi.local_set_size(),
																							/*int batchCount = */ nprojs_);

			gpu::sync();
			
			assert(status == CUBLAS_STATUS_SUCCESS);
			
		}
#endif

		return sphere_phi_all;
			
	}

	////////////////////////////////////////////////////////////////////////////////////////////		

	template <typename SpherePhiType, typename KpointType>
	void apply(SpherePhiType & sphere_vnlphi, states::orbital_set<basis::real_space, complex> & vnlphi, KpointType const & kpoint) const {

		CALI_CXX_MARK_SCOPE("projector_all::apply");

		gpu::run(vnlphi.local_set_size(), max_sphere_size_, nprojs_,
						 [sgr = begin(sphere_vnlphi), gr = begin(vnlphi.hypercubic()), poi = begin(points_), pos = begin(positions_), kpoint, empty = begin(locally_empty_)] GPU_LAMBDA (auto ist, auto ipoint, auto iproj){
							 if(not empty[iproj] and poi[iproj][ipoint][0] >= 0){
								 auto phase = polar(1.0, -dot(kpoint, pos[iproj][ipoint]));
								 gpu::atomic::add(&gr[poi[iproj][ipoint][0]][poi[iproj][ipoint][1]][poi[iproj][ipoint][2]][ist], phase*sgr[iproj][ipoint][ist]);
							 }
						 });
	}

	////////////////////////////////////////////////////////////////////////////////////////////

	template <typename OcType, typename PhiType, typename GPhiType>
	struct force_term {
		OcType oc;
		PhiType phi;
		GPhiType gphi;
		constexpr auto operator()(int ist, int ip) const {
			return -2.0*oc[ist]*real(phi[ip][ist]*conj(gphi[ip][ist]));
		}
	};

	////////////////////////////////////////////////////////////////////////////////////////////

	template <typename Projections, typename Coeff, typename Occupations>
	struct energy_reduction {
		Projections proj;
		Coeff coe;
		Occupations occ;
		long spinor_size;

		GPU_FUNCTION auto operator()(long ist, long ilm, long iproj) const {
			auto ist_spinor = ist%spinor_size;
			auto pp = proj[iproj][ilm][ist];
			return real(conj(pp)*pp)*coe[iproj][ilm]*occ[ist_spinor];
		}
		
	};
	
	////////////////////////////////////////////////////////////////////////////////////////////

	template <typename KpointType, typename Occupations>
	double energy(states::orbital_set<basis::real_space, complex> const & phi, KpointType const & kpoint, Occupations const & occupations) const {
    
		gpu::array<complex, 3> sphere_phi_all({nprojs_, max_sphere_size_, phi.local_set_size()});
		gpu::array<complex, 3> projections_all({nprojs_, max_nlm_, phi.local_set_size()}, 0.0);

		{ CALI_CXX_MARK_SCOPE("projector::gather");
				
			gpu::run(phi.local_set_size(), max_sphere_size_, nprojs_,
							 [sgr = begin(sphere_phi_all), gr = begin(phi.hypercubic()), poi = begin(points_), pos = begin(positions_), kpoint] GPU_LAMBDA (auto ist, auto ipoint, auto iproj){
								 if(poi[iproj][ipoint][0] >= 0){
									 auto phase = polar(1.0, dot(kpoint, pos[iproj][ipoint]));
									 sgr[iproj][ipoint][ist] = phase*gr[poi[iproj][ipoint][0]][poi[iproj][ipoint][1]][poi[iproj][ipoint][2]][ist];
								 } else {
									 sgr[iproj][ipoint][ist] = complex(0.0, 0.0);
								 }
							 });
		}

#ifndef ENABLE_CUDA
		for(auto iproj = 0; iproj < nprojs_; iproj++){
			CALI_CXX_MARK_SCOPE("projector_gemm_1");

			if(locally_empty_[iproj]) continue;
			
			namespace blas = boost::multi::blas;
			blas::real_doubled(projections_all[iproj]) = blas::gemm(phi.basis().volume_element(), matrices_[iproj], blas::real_doubled(sphere_phi_all[iproj]));
		}
#else
		if(max_sphere_size_ > 0) {
			CALI_CXX_MARK_SCOPE("projector_gemm_1");			

			const double zero = 0.0;
			const double vol = phi.basis().volume_element();

			auto status = cublasDgemmStridedBatched(/*cublasHandle_t handle = */ boost::multi::cuda::cublas::context::get_instance().get(),
																							/*cublasOperation_t transa = */ CUBLAS_OP_N,
																							/*cublasOperation_t transb = */ CUBLAS_OP_N,
																							/*int m = */ 2*phi.local_set_size(),
																							/*int n = */ max_nlm_,
																							/*int k = */ max_sphere_size_,
																							/*const double *alpha = */ &vol,
																							/*const double *A = */ reinterpret_cast<double const *>(raw_pointer_cast(sphere_phi_all.data_elements())),
																							/*int lda = */ 2*phi.local_set_size(),
																							/*long long int strideA = */ 2*max_sphere_size_*phi.local_set_size(),
																							/*const double *B = */ raw_pointer_cast(matrices_.data_elements()),
																							/*int ldb = */ max_sphere_size_,
																							/*long long int strideB =*/ max_nlm_*max_sphere_size_,
																							/*const double *beta = */ &zero,
																							/*double *C = */ reinterpret_cast<double *>(raw_pointer_cast(projections_all.data_elements())),
																							/*int ldc = */ 2*phi.local_set_size(),
																							/*long long int strideC = */ 2*max_nlm_*phi.local_set_size(),
																							/*int batchCount = */ nprojs_);
			gpu::sync();
			
			assert(status == CUBLAS_STATUS_SUCCESS);
			
		}
#endif

		if(phi.basis().comm().size() > 1) {
			CALI_CXX_MARK_SCOPE("projector_all::project::reduce");
			phi.basis().comm().all_reduce_in_place_n(raw_pointer_cast(projections_all.data_elements()), projections_all.num_elements(), std::plus<>{});
		}
		
		auto en = gpu::run(gpu::reduce(phi.local_set_size()), gpu::reduce(max_nlm_), gpu::reduce(nprojs_),
											 energy_reduction<decltype(begin(projections_all)), decltype(begin(coeff_)), decltype(begin(occupations))>
											 {begin(projections_all), begin(coeff_), begin(occupations), phi.local_spinor_set_size()});
		
		if(phi.set_comm().size() > 1) {
			CALI_CXX_MARK_SCOPE("projector_all::project::reduce");
			phi.set_comm().all_reduce_in_place_n(&en, 1, std::plus<>{});
		}

		return en;
			
	}

	////////////////////////////////////////////////////////////////////////////////////////////

	
	template <typename PhiType, typename GPhiType, typename MetricType, typename OccsType>
	void force(PhiType & phi, GPhiType const & gphi, MetricType const & metric,
						 OccsType const & occs, vector3<double, covariant> const & vector_potential, gpu::array<vector3<double>, 1> & forces_non_local) const {

		CALI_CXX_MARK_FUNCTION;

		namespace blas = boost::multi::blas;

		gpu::array<typename PhiType::element_type, 3> sphere_phi_all({nprojs_, max_sphere_size_, phi.local_set_size()});
		gpu::array<typename GPhiType::element_type, 3> sphere_gphi_all({nprojs_, max_sphere_size_, phi.local_set_size()});
		gpu::array<complex, 3> projections_all({nprojs_, max_nlm_, phi.local_set_size()}, 0.0);
 
		{ CALI_CXX_MARK_SCOPE("projector_all::force::gather");
				
			gpu::run(phi.local_set_size(), max_sphere_size_, nprojs_,
							 [sgr = begin(sphere_phi_all), gsgr = begin(sphere_gphi_all), gr = begin(phi.hypercubic()), ggr = begin(gphi.hypercubic()), poi = begin(points_), pos = begin(positions_), kpoint = phi.kpoint() + vector_potential]
							 GPU_LAMBDA (auto ist, auto ipoint, auto iproj){
								 if(poi[iproj][ipoint][0] >= 0){
									 auto phase = polar(1.0, dot(kpoint, pos[iproj][ipoint]));
									 sgr[iproj][ipoint][ist] = phase*gr[poi[iproj][ipoint][0]][poi[iproj][ipoint][1]][poi[iproj][ipoint][2]][ist];
									 gsgr[iproj][ipoint][ist] = phase*ggr[poi[iproj][ipoint][0]][poi[iproj][ipoint][1]][poi[iproj][ipoint][2]][ist];
								 } else {
									 sgr[iproj][ipoint][ist] = 0.0;
									 gsgr[iproj][ipoint][ist] = {complex(0.0), complex(0.0), complex(0.0)};
								 }
							 });
		}

		for(auto iproj = 0; iproj < nprojs_; iproj++){
			if(locally_empty_[iproj]) continue;
			
			blas::real_doubled(projections_all[iproj]) = blas::gemm(phi.basis().volume_element(), matrices_[iproj], blas::real_doubled(sphere_phi_all[iproj]));
		}
			
		{ CALI_CXX_MARK_SCOPE("projector_force_scal"); 
			
			gpu::run(phi.local_set_size(), max_nlm_, nprojs_,
							 [proj = begin(projections_all), coeff = begin(coeff_)] GPU_LAMBDA (auto ist, auto ipj, auto iproj){
								 proj[iproj][ipj][ist] *= coeff[iproj][ipj];
							 });
		}

		if(phi.basis().comm().size() > 1) {
			phi.basis().comm().all_reduce_in_place_n(raw_pointer_cast(projections_all.data_elements()), projections_all.num_elements(), std::plus<>{});
		}
		
		for(auto iproj = 0; iproj < nprojs_; iproj++){		
			blas::real_doubled(sphere_phi_all[iproj]) = blas::gemm(1.0, blas::T(matrices_[iproj]), blas::real_doubled(projections_all[iproj]));
		}

		gpu::array<vector3<double, covariant>, 1> force(nprojs_, {0.0, 0.0, 0.0});
			
		for(auto iproj = 0; iproj < nprojs_; iproj++) {
			if(locally_empty_[iproj]) continue;
			
				CALI_CXX_MARK_SCOPE("projector_force_sum");
				force[iproj] = gpu::run(gpu::reduce(phi.local_set_size()), gpu::reduce(max_sphere_size_),
																force_term<decltype(begin(occs)), decltype(begin(sphere_phi_all[iproj])), decltype(begin(sphere_gphi_all[iproj]))>{begin(occs), begin(sphere_phi_all[iproj]), begin(sphere_gphi_all[iproj])});
		}

		for(auto iproj = 0; iproj < nprojs_; iproj++) {
			if(locally_empty_[iproj]) continue;
			
			forces_non_local[iatom_[iproj]] += phi.basis().volume_element()*metric.to_cartesian(force[iproj]);
		}
		
	}

	////////////////////////////////////////////////////////////////////////////////////////////
	// Calculates |cphi> += [Vnl, r] | phi>
	////////////////////////////////////////////////////////////////////////////////////////////	
	template <typename KpointType>
	void position_commutator(states::orbital_set<basis::real_space, complex> const & phi, states::orbital_set<basis::real_space, vector3<complex, covariant>> & cphi, KpointType const & kpoint) const {
		
		gpu::array<complex, 3> sphere_phi_all({nprojs_, max_sphere_size_, phi.local_set_size()});
		gpu::array<vector3<complex, contravariant>, 3> sphere_rphi_all({nprojs_, max_sphere_size_, phi.local_set_size()});		

		gpu::array<complex, 3> projections_all({nprojs_, max_nlm_, phi.local_set_size()}, 0.0);
		gpu::array<vector3<complex, contravariant>, 3> rprojections_all({nprojs_, max_nlm_, phi.local_set_size()}, zero<vector3<complex, contravariant>>());

		{ CALI_CXX_MARK_SCOPE("position_commutator::gather");
				
			gpu::run(phi.local_set_size(), max_sphere_size_, nprojs_,
							 [sphi = begin(sphere_phi_all), srphi = begin(sphere_rphi_all), gr = begin(phi.hypercubic()), poi = begin(points_), pos = begin(positions_), kpoint] GPU_LAMBDA (auto ist, auto ipoint, auto iproj){
								 if(poi[iproj][ipoint][0] >= 0){
									 auto rr = static_cast<vector3<double, contravariant>>(pos[iproj][ipoint]);
									 auto phase = polar(1.0, dot(kpoint, rr));
									 sphi[iproj][ipoint][ist] = phase*gr[poi[iproj][ipoint][0]][poi[iproj][ipoint][1]][poi[iproj][ipoint][2]][ist];
									 srphi[iproj][ipoint][ist] = rr*sphi[iproj][ipoint][ist];
								 } else {
									 sphi[iproj][ipoint][ist]     = complex(0.0, 0.0);
									 srphi[iproj][ipoint][ist][0] = complex(0.0, 0.0);
									 srphi[iproj][ipoint][ist][1] = complex(0.0, 0.0);
									 srphi[iproj][ipoint][ist][2] = complex(0.0, 0.0);
								 }
							 });
		}

	 	for(auto iproj = 0; iproj < nprojs_; iproj++){
			CALI_CXX_MARK_SCOPE("position_commutator_gemm_1");

			if(locally_empty_[iproj]) continue;
			
			namespace blas = boost::multi::blas;
			blas::real_doubled(projections_all[iproj]) = blas::gemm(phi.basis().volume_element(), matrices_[iproj], blas::real_doubled(sphere_phi_all[iproj]));

			auto rpa = rprojections_all[iproj].template reinterpret_array_cast<complex>(3).rotated().flatted().unrotated();
			auto sra = sphere_rphi_all[iproj].template reinterpret_array_cast<complex>(3).rotated().flatted().unrotated();

			blas::real_doubled(rpa) = blas::gemm(phi.basis().volume_element(), matrices_[iproj], blas::real_doubled(sra));
		}

    { CALI_CXX_MARK_SCOPE("position_commutator_scal");
				
      gpu::run(phi.local_set_size(), max_nlm_, nprojs_,
               [proj = begin(projections_all), rproj = begin(rprojections_all), coe = begin(coeff_)]
               GPU_LAMBDA (auto ist, auto ilm, auto iproj){
                 proj[iproj][ilm][ist] *= coe[iproj][ilm];
                 proj[iproj][ilm][ist] *= coe[iproj][ilm];								 
               });
		}

		if(phi.basis().comm().size() > 1) {
			phi.basis().comm().all_reduce_in_place_n(raw_pointer_cast(projections_all.data_elements()), projections_all.num_elements(), std::plus<>{});
			phi.basis().comm().all_reduce_in_place_n(raw_pointer_cast(rprojections_all.data_elements()), rprojections_all.num_elements(), std::plus<>{});
		}
		
		for(auto iproj = 0; iproj < nprojs_; iproj++) {
			CALI_CXX_MARK_SCOPE("position_commutator_gemm_2");

			if(locally_empty_[iproj]) continue;
			
			namespace blas = boost::multi::blas;
			blas::real_doubled(sphere_phi_all[iproj]) = blas::gemm(1., blas::T(matrices_[iproj]), blas::real_doubled(projections_all[iproj]));

			auto rpa = rprojections_all[iproj].template reinterpret_array_cast<complex>(3).rotated().flatted().unrotated();
			auto sra = sphere_rphi_all[iproj].template reinterpret_array_cast<complex>(3).rotated().flatted().unrotated();
			blas::real_doubled(sra) = blas::gemm(1., blas::T(matrices_[iproj]), blas::real_doubled(rpa));			
		}

		for(auto iproj = 0; iproj < nprojs_; iproj++){

			if(locally_empty_[iproj]) continue;
			
			gpu::run(phi.local_set_size(), max_sphere_size_,
							 [sgr = begin(sphere_phi_all), srphi = begin(sphere_rphi_all), gr = begin(cphi.hypercubic()), poi = begin(points_), iproj, pos = begin(positions_), kpoint, metric = phi.basis().cell().metric()]
							 GPU_LAMBDA (auto ist, auto ipoint){
								 if(poi[iproj][ipoint][0] >= 0){
									 auto rr = static_cast<vector3<double, contravariant>>(pos[iproj][ipoint]);
									 auto phase = polar(1.0, -dot(kpoint, rr));
									 auto commutator = phase*metric.to_covariant(srphi[iproj][ipoint][ist] - rr*sgr[iproj][ipoint][ist]);
									 gpu::atomic::add(&gr[poi[iproj][ipoint][0]][poi[iproj][ipoint][1]][poi[iproj][ipoint][2]][ist], commutator);
								 }
							 });
		}
	}

	////////////////////////////////////////////////////////////////////////////////////////////		


private:
			
	int nprojs_;
	long max_sphere_size_;
	int max_nlm_;
	gpu::array<vector3<int>, 2> points_;
	gpu::array<vector3<float, contravariant>, 2> positions_;
	gpu::array<double, 2> coeff_;
	gpu::array<double, 3> matrices_;
	gpu::array<int, 1> nlm_;
	gpu::array<int, 1> iatom_;
	gpu::array<bool, 1> locally_empty_;
  
};
  
}
}
#endif

#ifdef INQ_HAMILTONIAN_PROJECTOR_ALL_UNIT_TEST
#undef INQ_HAMILTONIAN_PROJECTOR_ALL_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace inq::magnitude;
	using namespace Catch::literals;
	
}
#endif
