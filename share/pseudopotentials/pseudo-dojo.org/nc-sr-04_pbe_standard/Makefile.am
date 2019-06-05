## Process this file with automake to produce Makefile.in

## Copyright (C) 2015 X. Andrade
##
## This program is free software; you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation; either version 2, or (at your option)
## any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## You should have received a copy of the GNU General Public License
## along with this program; if not, write to the Free Software
## Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
## 02110-1301, USA.
##

pdpbestddir = $(pkgdatadir)/pseudopotentials/pseudo-dojo.org/nc-sr-04_pbe_standard/

dist_pdpbestd_DATA =                \
	Ag.upf                      \
	Al.upf                      \
	Al.upf.spacing              \
	Ar.upf                      \
	Ar.upf.spacing              \
	As.upf                      \
	Au.upf                      \
	Ba.upf                      \
	Be.upf                      \
	Be.upf.spacing              \
	Bi.upf                      \
	Br.upf                      \
	B.upf                       \
	B.upf.spacing               \
	Ca.upf                      \
	Cd.upf                      \
	Cl.upf                      \
	Cl.upf.spacing              \
	Co.upf                      \
	Cr.upf                      \
	Cs.upf                      \
	C.upf                       \
	C.upf.spacing               \
	Cu.upf                      \
	D.upf                       \
	D.upf.spacing               \
	Fe.upf                      \
	F.upf                       \
	F.upf.spacing               \
	Ga.upf                      \
	Ga.upf.spacing              \
	Ge.upf                      \
	He.upf                      \
	He.upf.spacing              \
	Hf.upf                      \
	Hg.upf                      \
	H.upf                       \
	H.upf.spacing               \
	In.upf                      \
	Ir.upf                      \
	I.upf                       \
	Kr.upf                      \
	K.upf                       \
	La.upf                      \
	Li.upf                      \
	Li.upf.spacing              \
	Lu.upf                      \
	Mg.upf                      \
	Mn.upf                      \
	Mo.upf                      \
	Na.upf                      \
	Na.upf.spacing              \
	Nb.upf                      \
	Ne.upf                      \
	Ne.upf.spacing              \
	Ni.upf                      \
	N.upf                       \
	N.upf.spacing               \
	Os.upf                      \
	O.upf                       \
	O.upf.spacing               \
	Pb.upf                      \
	Pd.upf                      \
	Po.upf                      \
	Pt.upf                      \
	P.upf                       \
	Rb.upf                      \
	Re.upf                      \
	Rh.upf                      \
	Rn.upf                      \
	Ru.upf                      \
	Sb.upf                      \
	Sc.upf                      \
	Se.upf                      \
	Si.upf                      \
	Si.upf.spacing              \
	Sn.upf                      \
	Sr.upf                      \
	S.upf                       \
	S.upf.spacing               \
	T.upf                       \
	T.upf.spacing               \
	Ta.upf                      \
	Tc.upf                      \
	Te.upf                      \
	Ti.upf                      \
	Tl.upf                      \
	V.upf                       \
	W.upf                       \
	Xe.upf                      \
	Y.upf                       \
	Zn.upf                      \
	Zr.upf

# In order to get `make check' working with a `VPATH' build (e. g. to use
# the `distcheck' target) the pseudopotential files have to be accessible
# in the build tree.
check:
	for f in $(dist_pdpbestd_DATA); do \
	  if [ ! -e $$f -o $$f -ot $(srcdir)/$$f ]; then \
	    cp $(srcdir)/$$f $$f; \
	  fi \
	done

clean-local:
	@if [ -n "$(top_builddir)" -a "$(top_builddir)" != "$(top_srcdir)" ]; then \
	  rm -f $(dist_pdpbestd_DATA); \
	fi

.PHONY: check
