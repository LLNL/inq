bin_SCRIPTS = inc++
CLEANFILES = $(bin_SCRIPTS)

edit = sed \
  -e 's|@libdir[@]|$(pkgdatadir)|g' \
  -e 's|@pkgdatadir[@]|$(pkgdatadir)|g'

inc++: inc++.bash
	$(edit) inc++.bash > inc++
