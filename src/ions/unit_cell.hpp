/* -*- indent-tabs-mode: t -*- */

#ifndef INQ__IONS__UNIT_CELL
#define INQ__IONS__UNIT_CELL

// Copyright (C) 2019-2023 Lawrence Livermore National Security, LLC., Xavier Andrade, Alfredo A. Correa
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include <stdexcept>

#include <math/vector3.hpp>
#include <valarray>
#include <array>

namespace inq {
namespace ions {

  class unit_cell{
		using vector_type = vector3<double>;
  private:
    vector_type lattice_[3];
    vector_type reciprocal_[3];
    double volume_;
		int periodicity_;
		
  public:
		
    unit_cell(vector3<double> const& a0, vector3<double> const& a1, vector3<double> const& a2, int arg_periodicity = 3){
			
			periodicity_ = arg_periodicity;
		
			lattice_[0] = a0;
			lattice_[1] = a1;
			lattice_[2] = a2;
			volume_ = dot(a0, cross(a1, a2));

			if(volume_ < 1e-10) throw std::runtime_error("inq error: the lattice volume " + std::to_string(volume_) + " is not positive");
		
			reciprocal_[0] = 2.0*M_PI/volume_*cross(a1, a2);
			reciprocal_[1] = 2.0*M_PI/volume_*cross(a2, a0);
			reciprocal_[2] = 2.0*M_PI/volume_*cross(a0, a1);
		}
		
		template<class lat_type>
		unit_cell(const lat_type & lat, int periodicity = 3):
			unit_cell(vector_type{lat[0][0], lat[0][1], lat[0][2]}, vector_type{lat[1][0], lat[1][1], lat[1][2]}, vector_type{lat[2][0], lat[2][1], lat[2][2]}, periodicity){
		}

		vector_type const& operator[](int ii) const {
			assert(ii >= 0 and ii < 3);         
			return lattice_[ii];
		}
    
    vector_type const& lattice(int ii) const {
			assert(ii >= 0 and ii < 3);
			return lattice_[ii];
		}
		
		vector_type const& reciprocal(int ii) const {
			assert(ii >= 0 and ii < 3);
			return reciprocal_[ii];
		}
		
		auto enlarge(int factor) const {
			return unit_cell(factor*lattice_[0], factor*lattice_[1], factor*lattice_[2], periodicity_);
		}

		auto enlarge(vector3<int> factor) const {
			return unit_cell(factor[0]*lattice_[0], factor[1]*lattice_[1], factor[2]*lattice_[2], periodicity_);
		}
			
    double volume() const { return volume_; }

    template <class output_stream>
    void info(output_stream & out) const {
      out << "UNIT CELL:" << std::endl;
      out << "  Lattice vectors [b] = " << lattice_[0] << std::endl;
      out << "                        " << lattice_[1] << std::endl;
      out << "                        " << lattice_[2] << std::endl;
      out << "  Volume [b^3]        = " << volume_ << std::endl;
      out << std::endl;
    }

		template<class OStream>
		friend OStream& operator<<(OStream& os, unit_cell const& self){
			self.info(os);
			return os;
	  }
		
		////////////////////////////////////////////////////////////////////////////////

    bool operator==(const unit_cell& c) const {
			return ( lattice_[0]==c.lattice_[0] && lattice_[1]==c.lattice_[1] && lattice_[2]==c.lattice_[2] );
		}
			
    bool operator!=(const unit_cell& c) const {
			return ! ( *this == c );
		}

		auto diagonal_length() const {
			return sqrt(norm(lattice_[0]) + norm(lattice_[1]) + norm(lattice_[2]));
		}

		auto periodicity() const {
			return periodicity_;
		}
		
		class cell_metric {

			vector3<double, cartesian> lat_[3];
			vector3<double, cartesian> rlat_[3];

		public:

			template <class AType, class BType>
			cell_metric(AType const & aa, BType const & bb){
				for(int ii = 0; ii < 3; ii++){
					for(int jj = 0; jj < 3; jj++){
						lat_[ii][jj] = aa[ii][jj];
						rlat_[ii][jj] = bb[ii][jj];
					}
				}

			}

			template <class Type1, class Space1, class Type2, class Space2>
			GPU_FUNCTION auto dot(vector3<Type1, Space1> const & vv1, vector3<Type2, Space2> const & vv2) const {
				return to_cartesian(vv1).dot(to_cartesian(vv2));
			}
			
			template <class Type, class Space>
			GPU_FUNCTION auto norm(vector3<Type, Space> const & vv) const {
				return real(dot(vv, vv));
			}

			template <class Type, class Space>
			GPU_FUNCTION auto length(vector3<Type, Space> const & vv) const {
				return sqrt(norm(vv));
			}

			template <class Type>
			GPU_FUNCTION auto to_covariant(vector3<Type, covariant> const & vv) const {
				return vv;
			}
			
			template <class Type>
			GPU_FUNCTION auto to_covariant(vector3<Type, cartesian> const & vv) const {
				return vector3<Type, covariant>{lat_[0].dot(vv), lat_[1].dot(vv), lat_[2].dot(vv)};
			}

			template <class Type>
			GPU_FUNCTION auto to_covariant(vector3<Type, contravariant> const & vv) const {
				return to_covariant(to_cartesian(vv));
			}

			template <class Type>
			GPU_FUNCTION auto to_contravariant(vector3<Type, contravariant> const & vv) const {
				return vv;
			}
			
			template <class Type>
			GPU_FUNCTION auto to_contravariant(vector3<Type, covariant> const & vv) const {
				return to_contravariant(to_cartesian(vv));
			}

			template <class Type>
			GPU_FUNCTION auto to_contravariant(vector3<Type, cartesian> const & vv) const {
				return vector3<Type, contravariant>{rlat_[0].dot(vv), rlat_[1].dot(vv), rlat_[2].dot(vv)}/(2.0*M_PI);
			}

			template <class Type>
			GPU_FUNCTION auto to_cartesian(vector3<Type, cartesian> const & vv) const {
				return vv;
			}
			
			template <class Type>
			GPU_FUNCTION auto to_cartesian(vector3<Type, contravariant> const & vv) const {
				return lat_[0]*vv[0] + lat_[1]*vv[1] + lat_[2]*vv[2];
			}
			
			template <class Type>
			GPU_FUNCTION auto to_cartesian(vector3<Type, covariant> const & vv) const {
				return (rlat_[0]*vv[0] + rlat_[1]*vv[1] + rlat_[2]*vv[2])/(2.0*M_PI);
			}
			
		};

		auto metric() const {
			return cell_metric{lattice_, reciprocal_};
		}

    bool contains(vector3<double, contravariant> point) const {
			return point[0] >= -0.5 && point[0] < 0.5 && point[1] >= -0.5 && point[1] < 0.5 && point[2] >= -0.5 && point[2] <= 0.5;
		}

		bool contains(vector3<double> point) const {
			return contains(metric().to_contravariant(point));
		}
		
		auto position_in_cell(vector3<double, contravariant> crystal_pos) const {
			for(int idir = 0; idir < 3; idir++) {
				crystal_pos[idir] -= floor(crystal_pos[idir]);
				if(crystal_pos[idir] >= 0.5) crystal_pos[idir] -= 1.0;
			}
			assert(contains(crystal_pos));
			return crystal_pos;
		}
		
		auto position_in_cell(vector3<double> const & pos) const {
			auto crystal_pos = metric().to_contravariant(pos);
			return metric().to_cartesian(position_in_cell(crystal_pos));
		}

		auto is_orthogonal() const {
			return fabs(dot(lattice_[0], lattice_[1])) < 1e-8 and fabs(dot(lattice_[1], lattice_[2])) < 1e-8 and fabs(dot(lattice_[2], lattice_[0])) < 1e-8;
		}

		auto is_cartesian() const {
			return is_orthogonal() and lattice_[0][1] < 1e-8 and lattice_[0][2] < 1e-8;
		}
		
  };

}
}
#endif

#ifdef INQ_IONS_UNIT_CELL_UNIT_TEST
#undef INQ_IONS_UNIT_CELL_UNIT_TEST

#include <catch2/catch_all.hpp>

TEST_CASE(INQ_TEST_FILE, INQ_TEST_TAG) {

	using namespace inq;
	using namespace Catch::literals;
	using Catch::Approx;
	
  {
    
    SECTION("Cubic cell"){
    
      ions::unit_cell cell(vector3<double>(10.0, 0.0, 0.0), vector3<double>(0.0, 10.0, 0.0), vector3<double>(0.0, 0.0, 10.0));

      CHECK(cell[0][0] == 10.0_a);
      CHECK(cell[0][1] ==  0.0_a);
      CHECK(cell[0][2] ==  0.0_a);
      CHECK(cell[1][0] ==  0.0_a);
      CHECK(cell[1][1] == 10.0_a);
      CHECK(cell[1][2] ==  0.0_a);
      CHECK(cell[2][0] ==  0.0_a);
      CHECK(cell[2][1] ==  0.0_a);
      CHECK(cell[2][2] == 10.0_a);
      
      CHECK(cell.lattice(0)[0] == 10.0_a);
      CHECK(cell.lattice(0)[1] ==  0.0_a);
      CHECK(cell.lattice(0)[2] ==  0.0_a);
      CHECK(cell.lattice(1)[0] ==  0.0_a);
      CHECK(cell.lattice(1)[1] == 10.0_a);
      CHECK(cell.lattice(1)[2] ==  0.0_a);
      CHECK(cell.lattice(2)[0] ==  0.0_a);
      CHECK(cell.lattice(2)[1] ==  0.0_a);
      CHECK(cell.lattice(2)[2] == 10.0_a);
    
      CHECK(cell.reciprocal(0)[0] == 0.6283185307_a);
      CHECK(cell.reciprocal(0)[1] ==  0.0_a);
      CHECK(cell.reciprocal(0)[2] ==  0.0_a);
      CHECK(cell.reciprocal(1)[0] ==  0.0_a);
      CHECK(cell.reciprocal(1)[1] == 0.6283185307_a);
      CHECK(cell.reciprocal(1)[2] ==  0.0_a);
      CHECK(cell.reciprocal(2)[0] ==  0.0_a);
      CHECK(cell.reciprocal(2)[1] ==  0.0_a);
      CHECK(cell.reciprocal(2)[2] == 0.6283185307_a);
    
      CHECK(cell.volume() == 1000.0_a);

      CHECK(!cell.contains(vector3<double>(5.0, 5.0, 5.0)));
      CHECK(!cell.contains(vector3<double>(-5.0, 5.0, 5.0)));
      CHECK(!cell.contains(vector3<double>(5.0, -5.0, 5.0)));
      CHECK(!cell.contains(vector3<double>(5.0, 5.0, -5.0)));

      CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(0.2, -0.5, 0.867))[0] == 2.0_a);
      CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(0.2, -0.5, 0.867))[1] == -5.0_a);
      CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(0.2, -0.5, 0.867))[2] == 8.67_a);

      CHECK(cell.metric().to_contravariant(vector3<double>(6.66, -3.77, 27.2))[0] == 0.666_a);
      CHECK(cell.metric().to_contravariant(vector3<double>(6.66, -3.77, 27.2))[1] == -0.377_a);
      CHECK(cell.metric().to_contravariant(vector3<double>(6.66, -3.77, 27.2))[2] == 2.72_a);

			auto in_cell = cell.position_in_cell(vector3<double>(6.66, 25.0, -18.33));

			CHECK(in_cell[0] == -3.34_a);
			CHECK(in_cell[1] == -5.00_a);
			CHECK(in_cell[2] == 1.67_a);

			{
				auto vv = cell.metric().to_contravariant(vector3<double, cartesian>{10.0, 10.0, 10.0});
				CHECK(vv[0] == 1.0_a);      
				CHECK(vv[1] == 1.0_a);
				CHECK(vv[2] == 1.0_a);
			}

			{
				auto vv = cell.metric().to_covariant(vector3<double, cartesian>{M_PI/10.0, M_PI/10.0, M_PI/10.0});
				CHECK(vv[0] == Approx(M_PI));       
				CHECK(vv[1] == Approx(M_PI));
				CHECK(vv[2] == Approx(M_PI));
			}

			{
				auto vv = vector3<double, contravariant>{1.0, 0.0, 0.0};
				CHECK(cell.metric().length(vv) == 10.0);            
				CHECK(dot(cell.metric().to_covariant(vv), vv) == 100.0_a);
				CHECK(cell.metric().dot(cell.metric().to_covariant(vv), vv) == 100.0_a);

			}

			CHECK(cell.is_orthogonal());
			CHECK(cell.is_cartesian());         
    }

    SECTION("Parallelepipedic cell"){

			ions::unit_cell cell(vector3<double>(28.62, 0.0, 0.0), vector3<double>(0.0, 90.14, 0.0), vector3<double>(0.0, 0.0, 12.31));
			
      CHECK(cell.lattice(0)[0] == 28.62_a);
      CHECK(cell.lattice(0)[1] ==  0.0_a);
      CHECK(cell.lattice(0)[2] ==  0.0_a);
      CHECK(cell.lattice(1)[0] ==  0.0_a);
      CHECK(cell.lattice(1)[1] == 90.14_a);
      CHECK(cell.lattice(1)[2] ==  0.0_a);
      CHECK(cell.lattice(2)[0] ==  0.0_a);
      CHECK(cell.lattice(2)[1] ==  0.0_a);
      CHECK(cell.lattice(2)[2] == 12.31_a);
    
      CHECK(cell.reciprocal(0)[0] == 0.2195382707_a);
      CHECK(cell.reciprocal(0)[1] ==  0.0_a);
      CHECK(cell.reciprocal(0)[2] ==  0.0_a);
      CHECK(cell.reciprocal(1)[0] ==  0.0_a);
      CHECK(cell.reciprocal(1)[1] == 0.0697047405_a);
      CHECK(cell.reciprocal(1)[2] ==  0.0_a);
      CHECK(cell.reciprocal(2)[0] ==  0.0_a);
      CHECK(cell.reciprocal(2)[1] ==  0.0_a);
      CHECK(cell.reciprocal(2)[2] == 0.5104131038_a);
    
      CHECK(cell.volume() == 31757.421708_a);

      CHECK(cell.contains(vector3<double>(5.0, 5.0, 5.0)));
      CHECK(cell.contains(vector3<double>(-5.0, 5.0, 5.0)));
      CHECK(cell.contains(vector3<double>(5.0, -5.0, 5.0)));
      CHECK(cell.contains(vector3<double>(5.0, 5.0, -5.0)));

      CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(0.2, -0.5, 0.867))[0] == 5.724_a);
      CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(0.2, -0.5, 0.867))[1] == -45.07_a);
      CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(0.2, -0.5, 0.867))[2] == 10.67277_a);

      CHECK(cell.metric().to_contravariant(vector3<double>(6.66, -203.77, 927.2))[0] == 0.2327044025_a);
      CHECK(cell.metric().to_contravariant(vector3<double>(6.66, -203.77, 927.2))[1] == -2.2605946306_a);
      CHECK(cell.metric().to_contravariant(vector3<double>(6.66, -203.77, 927.2))[2] == 75.3208773355_a);

			auto in_cell = cell.position_in_cell(vector3<double>(6.66, 225.0, -18.33));

			CHECK(in_cell[0] == 6.66_a);
			CHECK(in_cell[1] == 44.72_a);
			CHECK(in_cell[2] == -6.02_a);
			
			{
				auto vv = cell.metric().to_contravariant(vector3<double, cartesian>{28.62, 90.14, 12.31});
				CHECK(vv[0] == 1.0_a);      
				CHECK(vv[1] == 1.0_a);
				CHECK(vv[2] == 1.0_a);

				auto vv2 = cell.metric().to_cartesian(vv);
				CHECK(vv2[0] == 28.62_a);       
				CHECK(vv2[1] == 90.14_a);
				CHECK(vv2[2] == 12.31_a);

				CHECK(norm(vv2) == Approx(dot(vv, cell.metric().to_covariant(vv))));
			}

			{
				auto vv1 = cell.metric().to_covariant(vector3<double, cartesian>{M_PI/28.62, 0.0, 0.0});
				CHECK(vv1[0] == Approx(M_PI));      
				CHECK(vv1[1] == Approx(0.0));
				CHECK(vv1[2] == Approx(0.0));

				auto vv2 = cell.metric().to_covariant(vector3<double, cartesian>{0.0, M_PI/90.14, 0.0});
				CHECK(vv2[0] == Approx(0.0));       
				CHECK(vv2[1] == Approx(M_PI));
				CHECK(vv2[2] == Approx(0.0));

				auto vv3 = cell.metric().to_covariant(vector3<double, cartesian>{0.0, 0.0, M_PI/12.31});
				CHECK(vv3[0] == Approx(0.0));       
				CHECK(vv3[1] == Approx(0.0));
				CHECK(vv3[2] == Approx(M_PI));

				CHECK(cell.metric().dot(vv1, vector3<double, cartesian>{28.62, 0.0, 0.0}) == Approx(M_PI));
				
			}
						
			CHECK(cell.metric().length(vector3<double, contravariant>{1.0, 0.0, 0.0}) == 28.62);

			CHECK(cell.is_orthogonal());
			CHECK(cell.is_cartesian());         
		}

    SECTION("Rotated cell"){

			double ll = 10;
			double lv[3][3] = {{ll/sqrt(2.0), ll/sqrt(2.0), 0.0}, {-ll/sqrt(2.0), ll/sqrt(2.0), 0.0}, {0.0, 0.0, ll}};
      ions::unit_cell cell(lv);
      
      CHECK(cell.lattice(0)[0] == 7.0710678119_a);
      CHECK(cell.lattice(0)[1] == 7.0710678119_a);
      CHECK(cell.lattice(0)[2] == 0.0_a);
      CHECK(cell.lattice(1)[0] == -7.0710678119_a);
      CHECK(cell.lattice(1)[1] == 7.0710678119_a);
      CHECK(cell.lattice(1)[2] == 0.0_a);
      CHECK(cell.lattice(2)[0] == 0.0_a);
      CHECK(cell.lattice(2)[1] == 0.0_a);
      CHECK(cell.lattice(2)[2] == 10.0_a);
			
      CHECK(cell.reciprocal(0)[0] == 0.4442882938_a);
      CHECK(cell.reciprocal(0)[1] == 0.4442882938_a);
      CHECK(cell.reciprocal(0)[2] == 0.0_a);
      CHECK(cell.reciprocal(1)[0] == -0.4442882938_a);
      CHECK(cell.reciprocal(1)[1] == 0.4442882938_a);
      CHECK(cell.reciprocal(1)[2] == 0.0_a);
      CHECK(cell.reciprocal(2)[0] == 0.0_a);
      CHECK(cell.reciprocal(2)[1] == 0.0_a);
      CHECK(cell.reciprocal(2)[2] == 0.6283185307_a);

      CHECK(cell.volume() == 1000.0_a);

			CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(1.0, 0.0, 0.0))[0] == 7.0710678119_a);
			CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(1.0, 0.0, 0.0))[1] == 7.0710678119_a);
			CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(1.0, 0.0, 0.0))[2] == 0.0_a);

			CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(0.0, 1.0, 0.0))[0] == -7.0710678119_a);
			CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(0.0, 1.0, 0.0))[1] == 7.0710678119_a);
			CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(0.0, 1.0, 0.0))[2] == 0.0_a);

			CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(0.0, 0.0, 1.0))[0] == 0.0_a);
			CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(0.0, 0.0, 1.0))[1] == 0.0_a);
			CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(0.0, 0.0, 1.0))[2] == 10.0_a);

			CHECK(dot(cell.reciprocal(0), cell.lattice(0)) == Approx(2.0*M_PI));
			CHECK(dot(cell.reciprocal(1), cell.lattice(0)) < 1e-14);
			CHECK(dot(cell.reciprocal(2), cell.lattice(0)) < 1e-14);
			CHECK(dot(cell.reciprocal(0), cell.lattice(1)) < 1e-14);
			CHECK(dot(cell.reciprocal(1), cell.lattice(1)) == 2.0*M_PI);
			CHECK(dot(cell.reciprocal(2), cell.lattice(1)) < 1e-14);
			CHECK(dot(cell.reciprocal(0), cell.lattice(2)) < 1e-14);
			CHECK(dot(cell.reciprocal(1), cell.lattice(2)) < 1e-14);
			CHECK(dot(cell.reciprocal(2), cell.lattice(2)) == 2.0*M_PI);        

		CHECK(cell.metric().to_contravariant(vector3<double>(7.0710678119, 7.0710678119, 0.0))[0] == 1.0_a);
		CHECK(cell.metric().to_contravariant(vector3<double>(7.0710678119, 7.0710678119, 0.0))[1] == 0.0_a.margin(1e-12) );
		CHECK(cell.metric().to_contravariant(vector3<double>(7.0710678119, 7.0710678119, 0.0))[2] == 0.0_a.margin(1e-12) );

			CHECK(cell.is_orthogonal());
			CHECK(not cell.is_cartesian());
    }
		
    SECTION("Non-orthogonal cell"){

			double lv[3][3] = {{6.942, 8.799, 4.759}, {9.627, 7.092, 4.819}, {4.091, 0.721, 1.043}};
      ions::unit_cell cell(lv);
      
      CHECK(cell.lattice(0)[0] == 6.942_a);
      CHECK(cell.lattice(0)[1] == 8.799_a);
      CHECK(cell.lattice(0)[2] == 4.759_a);
      CHECK(cell.lattice(1)[0] == 9.627_a);
      CHECK(cell.lattice(1)[1] == 7.092_a);
      CHECK(cell.lattice(1)[2] == 4.819_a);
      CHECK(cell.lattice(2)[0] == 4.091_a);
      CHECK(cell.lattice(2)[1] == 0.721_a);
      CHECK(cell.lattice(2)[2] == 1.043_a);
    
      CHECK(cell.reciprocal(0)[0] == 3.3736397602_a);
      CHECK(cell.reciprocal(0)[1] == 8.3200742872_a);
      CHECK(cell.reciprocal(0)[2] == -18.9840209206_a);
      CHECK(cell.reciprocal(1)[0] == -4.942140131_a);
      CHECK(cell.reciprocal(1)[1] == -10.517582818_a);
      CHECK(cell.reciprocal(1)[2] == 26.6552948109_a);
      CHECK(cell.reciprocal(2)[0] == 7.4410562534_a);
      CHECK(cell.reciprocal(2)[1] == 10.6318294029_a);
      CHECK(cell.reciprocal(2)[2] == -30.5117208294_a);
    
      CHECK(cell.volume() == 7.305321831_a);

			CHECK(dot(cell.reciprocal(0), cell.lattice(0)) == Approx(2.0*M_PI));
			CHECK(dot(cell.reciprocal(1), cell.lattice(0)) < 1e-12);
			CHECK(dot(cell.reciprocal(2), cell.lattice(0)) < 1e-12);
			CHECK(dot(cell.reciprocal(0), cell.lattice(1)) < 1e-12);
			CHECK(dot(cell.reciprocal(1), cell.lattice(1)) == Approx(2.0*M_PI));
			CHECK(dot(cell.reciprocal(2), cell.lattice(1)) < 1e-12);
			CHECK(dot(cell.reciprocal(0), cell.lattice(2)) < 1e-12);
			CHECK(dot(cell.reciprocal(1), cell.lattice(2)) < 1e-12);
			CHECK(dot(cell.reciprocal(2), cell.lattice(2)) == Approx(2.0*M_PI));
			
      CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(0.2, -0.5, 0.867))[0] == 0.121797_a);
      CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(0.2, -0.5, 0.867))[1] == -1.161093_a);
      CHECK(cell.metric().to_cartesian(vector3<double, contravariant>(0.2, -0.5, 0.867))[2] == -0.553419_a);

      CHECK(cell.metric().to_contravariant(vector3<double>(0.66, -23.77, 2.72))[0] == -39.3396165136_a);
      CHECK(cell.metric().to_contravariant(vector3<double>(0.66, -23.77, 2.72))[1] == 50.8091863243_a);
      CHECK(cell.metric().to_contravariant(vector3<double>(0.66, -23.77, 2.72))[2] == -52.6483546581_a);

      CHECK(!cell.contains(vector3<double, contravariant>(0.5, 0.5, 0.5)));
			CHECK(!cell.contains(vector3<double, contravariant>(1.5, 0.5, 0.5)));
      CHECK(!cell.contains(vector3<double, contravariant>(0.5, -0.1, 0.0)));
      CHECK(!cell.contains(vector3<double, contravariant>(0.5, 0.5, -1.0)));

			{
				auto vv = cell.metric().to_contravariant(vector3<double, cartesian>{9.627, 7.092, 4.819});
				CHECK(fabs(vv[0]) < 1e-12);     
				CHECK(vv[1] == 1.0_a);
				CHECK(fabs(vv[2]) < 1e-12);

				auto vv2 = cell.metric().to_cartesian(vv);
				CHECK(vv2[0] == 9.627_a);       
				CHECK(vv2[1] == 7.092_a);
				CHECK(vv2[2] == 4.819_a);

				CHECK(norm(vv2) == Approx(dot(vv, cell.metric().to_covariant(vv))));
			}

			CHECK(not cell.is_orthogonal());
			CHECK(not cell.is_cartesian());
    }
  }
}
#endif
