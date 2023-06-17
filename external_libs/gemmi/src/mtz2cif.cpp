// Copyright 2020-2022 Global Phasing Ltd.

#include <gemmi/mtz2cif.hpp>

#include <climits>       // for INT_MIN
#include <algorithm>     // for all_of
#include <set>
#include <unordered_map>

#include <gemmi/asudata.hpp>   // for calculate_hkl_value_correlation
#include <gemmi/eig3.hpp>      // for eigen_decomposition
#include <gemmi/sprintf.hpp>   // for snprintf_z, to_str
#include <gemmi/atox.hpp>      // for read_word
#include <gemmi/version.hpp>   // for GEMMI_VERSION

namespace gemmi {

namespace {

int find_column_index(const std::string& column, const Mtz& mtz) {
  int idx = -1;
  for (const std::string& label : split_str(column, '|')) {
    for (size_t i = 0; i != mtz.columns.size(); ++i) {
      if (mtz.columns[i].label == label) {
        if (idx == -1)
          idx = (int) i;
        else
          mtz.warn("Column label duplicated: " + label);
      }
    }
    if (idx != -1)
      break;
  }
  return idx;
}

int check_format(const std::string& fmt) {
  // expected format: [#_+-]?\d*(\.\d+)?[fFgGeEc]
  int min_width = 0;
  if (fmt.find('%') != std::string::npos)
    fail("Specify format without %. Got: " + fmt);
  const char* p = fmt.c_str();
  if (*p == '_' || *p == '+' || *p == '-' || *p == '#')
   ++p;
  if (is_digit(*p)) {
    min_width = *p++ - '0';
    if (is_digit(*p)) // two digits of width number max
      min_width = min_width * 10 + (*p++ - '0');
  }
  if (*p == '.' && is_digit(*(p+1))) {
    p += 2;
    if (is_digit(*p)) // two digits of precision numbers max
      ++p;
  }
  if (!std::isalpha(*p) || *(p+1) != '\0')
    fail("wrong format : " + fmt + "\nCorrect examples: g, .4f, 12.5e");
  char c = alpha_up(*p);
  if (c != 'F' && c != 'G' && c != 'E')
    fail("expected floating-point format, got: " + fmt);
  if (min_width > 32)
    fail("the width exceeds 32: " + fmt);
  return min_width;
}

#define WRITE(...) os.write(buf, snprintf_z(buf, 255, __VA_ARGS__))

void write_cell_and_symmetry(const std::string& entry_id,
                             const UnitCell& cell, double* rmsds, const SpaceGroup* sg,
                             char* buf, std::ostream& os) {
  os << "_cell.entry_id " << entry_id << '\n';
  WRITE("_cell.length_a    %8.4f\n", cell.a);
  if (rmsds && rmsds[0] != 0.)
    WRITE("_cell.length_a_esd %7.3f\n", rmsds[0]);
  WRITE("_cell.length_b    %8.4f\n", cell.b);
  if (rmsds && rmsds[1] != 0.)
    WRITE("_cell.length_b_esd %7.3f\n", rmsds[1]);
  WRITE("_cell.length_c    %8.4f\n", cell.c);
  if (rmsds && rmsds[2] != 0.)
    WRITE("_cell.length_c_esd %7.3f\n", rmsds[2]);
  WRITE("_cell.angle_alpha %8.4f\n", cell.alpha);
  if (rmsds && rmsds[3] != 0.)
    WRITE("_cell.angle_alpha_esd %7.3f\n", rmsds[3]);
  WRITE("_cell.angle_beta  %8.4f\n", cell.beta);
  if (rmsds && rmsds[4] != 0.)
    WRITE("_cell.angle_beta_esd %8.3f\n", rmsds[4]);
  WRITE("_cell.angle_gamma %8.4f\n", cell.gamma);
  if (rmsds && rmsds[5] != 0.)
    WRITE("_cell.angle_gamma_esd %7.3f\n", rmsds[5]);
  if (sg) {
    os << "\n_symmetry.entry_id " << entry_id << "\n"
          "_symmetry.space_group_name_H-M '" << sg->pdb_name() << "'\n"
          "_symmetry.Int_Tables_number " << sg->number << '\n';
    // could write _symmetry_equiv.pos_as_xyz, but would it be useful?
  }
}

enum Var { Dot=-1, Qmark=-2, Counter=-3, DatasetId=-4, Image=-5 };

// describes which MTZ column is to be translated to what mmCIF column
struct Trans {
  int col_idx;
  bool is_status = false;
  std::string tag;  // excluding category
  std::string format = "%g";
  int min_width = 0;
};

// state for parse_spec_line
struct SpecParserState {
  size_t verified_spec_size = 0;
  bool discard_next_line = false;
};

// adds results to recipe
void parse_spec_line(const MtzToCif& m2c,
                     std::vector<Trans>& recipe, const char* line,
                     const Mtz& mtz, SpecParserState& state) {
  Trans tr;
  const char* p = line;
  if (*p == '&') {
    if (state.discard_next_line)
      return;
  } else {
    state.verified_spec_size = recipe.size();
    state.discard_next_line = false;
  }
  bool optional = (*p == '?' || *p == '&');
  if (optional)
    ++p;
  std::string column = read_word(p, &p);

  char ctype = '\0';
  if (column[0] != '$') {
    std::string type = read_word(p, &p);
    if (type.size() != 1)
      fail("Spec error: MTZ type '" + type + "' is not one character,"
           "\nin line: " + line);
    ctype = type[0];
    if (m2c.less_anomalous > 0) {
      bool is_i_ano = (ctype == 'K' || ctype == 'M');
      if (m2c.less_anomalous == 1
          ? is_i_ano && mtz.count_type('J') != 0 && mtz.count_type('G') > 1 &&
                        m2c.staraniso_version.empty()
          : is_i_ano || ctype == 'G' || ctype == 'D' || ctype == 'L') {
        recipe.resize(state.verified_spec_size);
        state.discard_next_line = true;
        return;
      }
    }
  }

  tr.tag = read_word(p, &p);
  if (tr.tag[0] == '_' || tr.tag.find('.') != std::string::npos)
    fail("Spec error: expected tag part after _refln., got: " +
         tr.tag + "\nin line: " + line);

  if (column[0] == '$') {
    if (column.size() == 2 && column[1] == '.')
      tr.col_idx = Var::Dot;
    else if (column.size() == 2 && column[1] == '?')
      tr.col_idx = Var::Qmark;
    else if (mtz.is_merged())
      fail("Variables other than $. and $? can be used only for unmerged files");
    else if (column.compare(1, 7, "counter") == 0)
      tr.col_idx = Var::Counter;
    else if (column.compare(1, 7, "dataset") == 0)
      tr.col_idx = Var::DatasetId;
    else if (column.compare(1, 5, "image") == 0)
      tr.col_idx = Var::Image;
    else
      fail("Unknown variable in the spec: " + column);
  } else {
    if (column.find('{') != std::string::npos &&
        !recipe.empty() && recipe.back().col_idx >= 0)
      replace_all(column, "{prev}", mtz.columns[recipe.back().col_idx].label);
    tr.col_idx = find_column_index(column, mtz);
    if (tr.col_idx == -1) {
      if (!optional)
        fail("Column not found: " + column);
      recipe.resize(state.verified_spec_size);
      state.discard_next_line = true;
      return;
    }
    const Mtz::Column& col = mtz.columns[tr.col_idx];
    if (ctype != '*' && col.type != ctype)
      fail("Column ", col.label, " has type ", col.type, " not ", ctype);
  }

  std::string fmt = read_word(p, &p);
  if (!fmt.empty()) {
    if (fmt.size() == 1 && fmt[0] == 'S') {
      tr.is_status = true;
    } else {
      tr.min_width = check_format(fmt);
      tr.format = "%" + fmt;
      if (tr.format[1] == '_')
        tr.format[1] = ' ';
    }
  }

  recipe.push_back(tr);
}

void prepare_recipe(const MtzToCif& m2c, const Mtz& mtz,
                    std::vector<Trans>& recipe) {
  recipe.clear();
  SpecParserState state;
  if (!m2c.spec_lines.empty()) {
    for (const std::string& line : m2c.spec_lines)
      parse_spec_line(m2c, recipe, line.c_str(), mtz, state);
  } else {
    const char** lines = MtzToCif::default_spec(/*for_merged=*/mtz.batches.empty());
    for (; *lines != nullptr; ++lines)
      parse_spec_line(m2c, recipe, *lines, mtz, state);
  }
  if (recipe.empty())
    fail("empty translation recipe");
  for (size_t i = 0; i != recipe.size(); ++i)
    for (size_t j = i + 1; j != recipe.size(); ++j)
      if (recipe[i].tag == recipe[j].tag)
        fail("duplicated output tag: " + recipe[i].tag);
  // H, K, L must be the first columns in MTZ and are required in _refln
  for (int i = 2; i != -1; --i)
    if (!in_vector_f([&](const Trans& t) { return t.col_idx == i; }, recipe)) {
      Trans tr;
      tr.col_idx = i;
      tr.tag = "index_";
      tr.tag += "hkl"[i]; // h, k or l
      recipe.insert(recipe.begin(), tr);
    }
}

// data corresponding to one sweep (dataset) in unmerged MTZ file
struct SweepData {
  int id = 0;
  int batch_count = 0;
  int offset = 0;
  int crystal_id = 1;
  const Mtz::Batch* first_batch = nullptr;
  const Mtz::Dataset* dataset = nullptr;
};

struct SweepInfo {
  std::vector<SweepData> sweeps;
  std::unordered_map<int, int> sweep_indices;

  bool gather_sweep_data(const Mtz& mtz) {
    bool ok = true;
    int prev_number = INT_MIN;
    int prev_dataset = INT_MIN;
    SweepData sweep;
    for (const Mtz::Batch& batch : mtz.batches) {
      int dataset_id = batch.dataset_id();
      if (dataset_id != prev_dataset || batch.number != prev_number + 1) {
        prev_dataset = dataset_id;
        if (sweep.id != 0)
          sweeps.push_back(sweep);
        sweep.id++;
        sweep.batch_count = 0;
        sweep.first_batch = &batch;
        sweep.offset = batch.number - (batch.number % 100);
        try {
          sweep.dataset = &mtz.dataset(dataset_id);
        } catch (std::runtime_error&) {
          sweep.dataset = nullptr;
          ok = false;
          mtz.warn("Reference to absent dataset: " + std::to_string(dataset_id));
        }
      }
      sweep_indices.emplace(batch.number, sweep.id - 1);
      sweep.batch_count++;
      prev_number = batch.number;
    }
    if (sweep.id != 0)
      sweeps.push_back(sweep);
    return ok;
  }

  const SweepData* get_sweep_data(int batch_number) const {
    return &sweeps[sweep_indices.at(batch_number)];
  }
};

void write_special_marker_if_requested(const MtzToCif& m2c,
                                       std::ostream& os, bool merged) {
  if (!m2c.write_special_marker_for_pdb)
    return;
  os << "### IF YOU MODIFY THIS FILE, REMOVE THIS SIGNATURE: ###\n";
  std::string desc;
  if (!m2c.gemmi_run_from.empty())
    desc = " 'run from " + m2c.gemmi_run_from + "'";
  if (!merged || m2c.staraniso_version.empty()) {
    os << "_software.pdbx_ordinal 1\n"
          "_software.classification 'data extraction'\n"
          "_software.name gemmi\n"
          "_software.version " GEMMI_VERSION "\n";
    if (!desc.empty())
      os << "_software.description" << desc << '\n';
  } else {
    os << "loop_\n"
          "_software.pdbx_ordinal\n"
          "_software.classification\n"
          "_software.name\n"
          "_software.version\n";
    if (!desc.empty())
      os << "_software.description\n";
    os << "1 'data extraction' gemmi " GEMMI_VERSION << desc << '\n';
    // STARANISO here tells that intensities were scaled anisotropically.
    os << "2 'data scaling' STARANISO '" << m2c.staraniso_version
       << (!desc.empty() ? "' .\n" : "'\n");
  }
  os << "_pdbx_audit_conform.dict_name mmcif_pdbx.dic\n"
        "_pdbx_audit_conform.dict_version 5.339\n"
        "_pdbx_audit_conform.dict_location "
        "https://mmcif.wwpdb.org/dictionaries/ascii/mmcif_pdbx_v50.dic\n"
        "### END OF SIGNATURE ###\n\n";
}

// Reorder eigenvalues and change signs of eigenvectors in the STARANISO way.
// It minimises the rotation angle from the basis vectors to the eigenvectors,
// which is equivalent to maximising the trace of the eigenvector matrix.
void reorder_staraniso_eigensystem(Mat33& vectors, double (&values)[3]) {
  const int8_t permut[6][3] = {{0,1,2}, {1,2,0}, {2,0,1}, {1,0,2}, {2,1,0}, {0,2,1}};
  const int8_t signs[8][3] = {{1,1,1}, {1,-1,-1}, {-1,1,-1}, {-1,-1,1},
                              {-1,-1,-1}, {-1,1,1}, {1,-1,1}, {1,1,-1}};
  double max_trace = -INFINITY;
  int permut_pos = 0, sign_pos = 0;
  bool det_neg = std::signbit(vectors.determinant());
  for (int i = 0; i < 6; ++i) {
    int jbase = det_neg == (i > 2) ? 0 : 4;
    const int8_t (&p)[3] = permut[i];
    for (int j = jbase; j < jbase+4; ++j) {
      double trace = 0.;
      for (int k = 0; k < 3; ++k)
        trace += signs[j][k] * vectors[k][p[k]];
      if (trace > max_trace) {
        max_trace = trace;
        permut_pos = i;
        sign_pos = j;
      }
    }
  }
  const int8_t (&p)[3] = permut[permut_pos];
  double tmp[3];
  for (int i = 0; i < 3; ++i) {
    for (int j = 0; j < 3; ++j)
      tmp[j] = signs[sign_pos][j] * vectors.a[i][p[j]];
    std::memcpy(vectors.a[i], tmp, sizeof(tmp));
  }
  for (int j = 0; j < 3; ++j)
    tmp[j] = values[p[j]];
  std::memcpy(values, tmp, sizeof(tmp));
}

// Get the first (non-zero) DWAVEL corresponding to an intensity, amplitude
// or sigma column from the template.
double get_wavelength(const Mtz& mtz, const std::vector<Trans>& spec) {
  for (const Trans& tr : spec) {
    if (tr.col_idx < 0)
      continue;
    const Mtz::Column& col = mtz.columns.at(tr.col_idx);
    if (col.type == 'F' || col.type == 'J' ||
        col.type == 'K' || col.type == 'G' ||
        col.type == 'Q' || col.type == 'M' || col.type == 'L') { // sigma
      double wavelength = mtz.dataset(col.dataset_id).wavelength;
      if (wavelength != 0.)
        return wavelength;
    }
  }
  return 0.;
}

const Trans* get_status_translation(const std::vector<Trans>& recipe) {
  for (const Trans& t: recipe)
    if (t.is_status)
      return &t;
  return nullptr;
}

void write_main_loop(const MtzToCif& m2c, const SweepInfo& sweep_info,
                     const Mtz& mtz, const std::vector<Trans>& recipe,
                     char (&buf)[256], std::ostream& os) {
  // prepare indices
  std::vector<int> value_indices;  // used for --skip_empty
  std::vector<int> sigma_indices;  // used for status 'x' and --skip-negative-sigi
  for (const Trans& tr : recipe) {
    if (tr.col_idx < 0)
      continue;
    const Mtz::Column& col = mtz.columns[tr.col_idx];
    if (m2c.skip_empty) {
      if (m2c.skip_empty_cols.empty() ? col.type != 'H' && col.type != 'I'
                                      : is_in_list(col.label, m2c.skip_empty_cols))
        value_indices.push_back(tr.col_idx);
    }
    if (col.type == 'Q' || col.type == 'L' || col.type == 'M')
      sigma_indices.push_back(tr.col_idx);
  }

  bool unmerged = !mtz.is_merged();
  os << "\nloop_\n";
  for (const Trans& tr : recipe) {
    os << (unmerged ? "_diffrn_refln." : "_refln.");
    if (m2c.with_comments && tr.col_idx >= 0) {
      WRITE("%-26s # ", tr.tag.c_str());
      const Mtz::Column& col = mtz.columns.at(tr.col_idx);
      const Mtz::Dataset& ds = mtz.dataset(col.dataset_id);
      // dataset is assigned to column only in merged MTZ
      if (unmerged)
        os << col.label;
      else
        WRITE("%-14s from dataset %s", col.label.c_str(), ds.dataset_name.c_str());
    } else {
      os << tr.tag;
    }
    os << '\n';
  }
  int batch_idx = find_column_index("BATCH", mtz);
  std::unordered_map<int, const Mtz::Batch*> batch_by_number;
  for (const Mtz::Batch& b : mtz.batches)
    batch_by_number.emplace(b.number, &b);

  int free_flag_value = m2c.free_flag_value;
  if (free_flag_value < 0) {
    // CCP4 uses flags 0,...N-1 (usually N=20), with default free set 0
    // PHENIX uses 0/1 flags with free set 1
    if (const Trans* tr_status = get_status_translation(recipe)) {
      int count = 0;
      for (float val : mtz.columns[tr_status->col_idx])
        if (val == 0.f)
          count++;
      free_flag_value = count < mtz.nreflections / 2 ? 0 : 1;
    }
  }

  auto write_int = [](char* p, int num) {
    //return snprintf_z(p, 32, "%d", num);
    std::string s = std::to_string(num);
    std::memcpy(p, s.data(), s.size());
    return s.size();
  };

  char* ptr = buf;
  for (int i = 0, idx = 0; i != mtz.nreflections; ++i) {
    const float* row = &mtz.data[i * mtz.columns.size()];
    if (m2c.trim > 0) {
      if (row[0] < -m2c.trim || row[0] > m2c.trim ||
          row[1] < -m2c.trim || row[1] > m2c.trim ||
          row[2] < -m2c.trim || row[2] > m2c.trim) {
        continue;
      }
    }
    if (!value_indices.empty())
      if (std::all_of(value_indices.begin(), value_indices.end(),
                      [&](int n) { return std::isnan(row[n]); }))
        continue;
    int batch_number = 0;
    const SweepData* sweep = nullptr;
    if (unmerged) {
      if (m2c.skip_negative_sigi &&
          std::any_of(sigma_indices.begin(), sigma_indices.end(),
                      [&](int n) { return row[n] < 0; }))
        continue;
      if (batch_idx == -1)
        fail("BATCH column not found");
      batch_number = (int) row[batch_idx];
      auto it = batch_by_number.find(batch_number);
      if (it == batch_by_number.end())
        fail("unexpected values in column BATCH");
      const Mtz::Batch& batch = *it->second;
      sweep = sweep_info.get_sweep_data(batch.number);
    }
    bool first = true;
    for (const Trans& tr : recipe) {
      if (first)
        first = false;
      else
        *ptr++ = ' ';
      static_assert(sizeof(buf) == 256, "sizeof buf");
      if (ptr - buf > 256 - 36) {
        os.write(buf, ptr - buf);
        ptr = buf;
      }
      if (tr.col_idx < 0) {
        switch (tr.col_idx) {
          case Var::Dot: *ptr++ = '.'; break;
          case Var::Qmark: *ptr++ = '?'; break;
          case Var::Counter: ptr += write_int(ptr, ++idx); break;
          case Var::DatasetId: ptr += write_int(ptr, sweep->id); break;
          case Var::Image: ptr += write_int(ptr, batch_number - sweep->offset); break;
        }
      } else {
        float v = row[tr.col_idx];
        if (tr.is_status) {
          char status = 'x';
          if (sigma_indices.empty() ||
              !std::all_of(sigma_indices.begin(), sigma_indices.end(),
                           [&](int n) { return std::isnan(row[n]); }))
            status = int(v) == free_flag_value ? 'f' : 'o';
          *ptr++ = status;
        } else if (std::isnan(v)) {
          // we checked that min_width <= 32
          for (int j = 1; j < tr.min_width; ++j)
            *ptr++ = ' ';
          *ptr++ = '?';
        } else {
#if defined(__GNUC__)
# pragma GCC diagnostic push
# pragma GCC diagnostic ignored "-Wformat-nonliteral"
#endif
          ptr += snprintf_z(ptr, 32, tr.format.c_str(), v);
#if defined(__GNUC__)
# pragma GCC diagnostic pop
#endif
        }
      }
    }
    *ptr++ = '\n';
  }
  os.write(buf, ptr - buf);
}

}  // anonymous namespace

void write_staraniso_b_in_mmcif(const SMat33<double>& b,
                                const std::string& entry_id,
                                char* buf,
                                std::ostream& os) {
  double eigenvalues[3];
  Mat33 eigenvectors = eigen_decomposition(b, eigenvalues);
  reorder_staraniso_eigensystem(eigenvectors, eigenvalues);
  // three mandatory items in _reflns
  os << "\n_reflns.entry_id " << entry_id
     << "\n_reflns.pdbx_ordinal 1"
        "\n_reflns.pdbx_diffrn_id 1";
  const char* prefix = "\n_reflns.pdbx_aniso_B_tensor_eigen";
  // STARANISO team decided that "the best course is to subtract the smallest
  // eigenvalue of the best-fit tensor from all eigenvalues. The smallest
  // eigenvalue is then exactly zero, and the tensor represents the anisotropy
  // correction that is actually applied to the data, there being no correction
  // in the strongest diffracting direction (this is the same as is done on the
  // UCLA server).". Actually here we have the minimum value already subtracted,
  // but due to the limited precision of the values stored as text in MTZ it's
  // not exactly 0, so we subtract the minimum value again.
  double min_eigv = std::min(eigenvalues[0], std::min(eigenvalues[1], eigenvalues[2]));
  for (int i = 0; i < 3; ++i) {
    WRITE("%svalue_%d %.5g", prefix, i+1, eigenvalues[i] - min_eigv);
    for (int j = 0; j < 3; ++j)
      WRITE("%svector_%d_ortho[%d] %.5g", prefix, i+1, j+1, eigenvectors[j][i]);
  }
  os << '\n';
}

void MtzToCif::write_cif(const Mtz& mtz, const Mtz* mtz2,
                         SMat33<double>* staraniso_b, std::ostream& os) {
  SweepInfo sweep_info;
  if (mtz2 && mtz.is_merged() == mtz2->is_merged())
    fail("If two MTZ files are given, one must be merged and one unmerged,\n"
         "got two ", mtz.is_merged() ? "merged" : "unmerged");
  const Mtz* merged = mtz.is_merged() ? &mtz : mtz2;
  const Mtz* unmerged = mtz.is_merged() ? mtz2 : &mtz;

  char buf[256];
  if (with_comments) {
    os << "# Converted by gemmi-mtz2cif " GEMMI_VERSION "\n";
    for (const Mtz* m : {merged, unmerged})
      if (m) {
        os << "# from " << (m->is_merged() ? "" : "un") << "merged MTZ: "
           << m->source_path << '\n';
        os << "#   title: " << m->title << '\n';
        if (with_history)
          for (size_t i = 0; i != m->history.size(); ++i)
            os << "#   history #" << i << ": " << m->history[i] << '\n';
      }
  }
  os << "data_" << (block_name ? block_name : "mtz");

  os << "\n\n_entry.id " << entry_id << "\n\n";

  // If we have user-provided spec file, we don't take responsibility
  // for the result. The spec in such case (write_special_marker_for_pdb==true)
  // is used for merged data only, so we can add the marker in unmerged block.
  if (spec_lines.empty() || !merged)
    write_special_marker_if_requested(*this, os, merged);

  std::vector<Trans> recipe;
  if (merged)
    prepare_recipe(*this, *merged, recipe);

  if (unmerged) {
    bool ok = sweep_info.gather_sweep_data(*unmerged);
    std::set<const Mtz::Dataset*> used_datasets;
    std::vector<std::string> crystal_names;
    // Prepare used_datasets, crystal_names and set SweepData::crystal_id.
    for (SweepData& sweep : sweep_info.sweeps)
      if (sweep.dataset) {
        used_datasets.insert(sweep.dataset);
        if (ok) {
          auto it = std::find(crystal_names.begin(), crystal_names.end(),
                              sweep.dataset->crystal_name);
          sweep.crystal_id = int(it - crystal_names.begin()) + 1;
          if (it == crystal_names.end())
            crystal_names.push_back(sweep.dataset->crystal_name);
        }
      }

    if (crystal_names.empty()) {
      os << "_exptl_crystal.id 1\n";
    } else {
      os << "loop_\n_exptl_crystal.id\n";
      for (size_t i = 0; i < crystal_names.size(); ++i)
        os << i+1 << " # " << crystal_names[i] << '\n';
    }

    bool scaled = (unmerged->column_with_label("SCALEUSED") != nullptr);

    os << "\nloop_\n"
          "_diffrn.id\n_diffrn.crystal_id\n_diffrn.details\n";
    for (const SweepData& sweep : sweep_info.sweeps) {
      os << sweep.id << ' ' << sweep.crystal_id << " '";
      if (scaled)
        os << "scaled ";
      os << "unmerged data'\n";
    }

    os << "\nloop_\n"
          "_diffrn_measurement.diffrn_id\n"
          "_diffrn_measurement.details\n";
    for (const SweepData& sweep : sweep_info.sweeps)
      os << sweep.id << " '" << sweep.batch_count << " frames'\n";

    os << "\nloop_\n"
          "_diffrn_radiation.diffrn_id\n"
          "_diffrn_radiation.wavelength_id\n";
    for (const SweepData& sweep : sweep_info.sweeps) {
      os << sweep.id << ' ';
      if (sweep.dataset)
        os << sweep.dataset->id;
      else
        os << '?';
      os << '\n';
    }

    os << "\nloop_\n"
          "_diffrn_radiation_wavelength.id\n"
          "_diffrn_radiation_wavelength.wavelength\n";
    for (const Mtz::Dataset* ds : used_datasets)
      os << ds->id << ' ' << ds->wavelength << '\n';
    os << '\n';

    if (enable_UB) {
      os << "loop_\n"
            "_diffrn_orient_matrix.diffrn_id\n"
            "_diffrn_orient_matrix.UB[1][1]\n"
            "_diffrn_orient_matrix.UB[1][2]\n"
            "_diffrn_orient_matrix.UB[1][3]\n"
            "_diffrn_orient_matrix.UB[2][1]\n"
            "_diffrn_orient_matrix.UB[2][2]\n"
            "_diffrn_orient_matrix.UB[2][3]\n"
            "_diffrn_orient_matrix.UB[3][1]\n"
            "_diffrn_orient_matrix.UB[3][2]\n"
            "_diffrn_orient_matrix.UB[3][3]\n";
      for (const SweepData& sweep : sweep_info.sweeps) {
        Mat33 u = sweep.first_batch->matrix_U();
        Mat33 b = sweep.first_batch->get_cell().calculate_matrix_B();
        Mat33 ub = u.multiply(b);
        WRITE("%d  %#g %#g %#g  %#g %#g %#g  %#g %#g %#g\n", sweep.id,
              ub.a[0][0], ub.a[0][1], ub.a[0][2],
              ub.a[1][0], ub.a[1][1], ub.a[1][2],
              ub.a[2][0], ub.a[2][1], ub.a[2][2]);
      }
      os << '\n';
    }
  } else {  // not unmerged, i.e. only merged data
    double w = std::isnan(wavelength) ? get_wavelength(mtz, recipe) : wavelength;
    if (w > 0.)
      os << "_diffrn_radiation_wavelength.id 1\n"
         << "_diffrn_radiation_wavelength.wavelength " << to_str(w) << "\n\n";
  }

  if (merged) {
    write_cell_and_symmetry(entry_id, mtz.get_cell(), nullptr, mtz.spacegroup, buf, os);
  } else {
    double rmsds[6];
    UnitCell cell = unmerged->get_average_cell_from_batch_headers(rmsds);
    write_cell_and_symmetry(entry_id, cell, rmsds, mtz.spacegroup, buf, os);
  }

  if (write_staraniso_tensor && staraniso_b)
    write_staraniso_b_in_mmcif(*staraniso_b, entry_id, buf, os);

  if (merged)
    write_main_loop(*this, sweep_info, *merged, recipe, buf, os);
  if (unmerged) {
    // if --depo flag is used, the spec file is for the merged data only
    if (write_special_marker_for_pdb)
      spec_lines.clear();
    prepare_recipe(*this, *unmerged, recipe);
    write_main_loop(*this, sweep_info, *unmerged, recipe, buf, os);
  }
}

void MtzToCif::write_cif_from_xds(const XdsAscii& xds, std::ostream& os) {
  char buf[256];
  if (with_comments) {
    os << "# Converted by gemmi-mtz2cif " GEMMI_VERSION "\n";
    os << "# from scaled unmerged XDS_ASCII: " << xds.source_path << '\n';
  }
  os << "data_" << (block_name ? block_name : "xds");
  os << "\n\n_entry.id " << entry_id << "\n\n";

  write_special_marker_if_requested(*this, os, false);

  os << "_exptl_crystal.id 1\n";

  os << "\nloop_\n_diffrn.id\n_diffrn.crystal_id\n_diffrn.details\n";
  for (const XdsAscii::Iset& iset : xds.isets)
    // We could write iset.input_info as details, but then local paths
    // could end up in the deposition.
    os << iset.id << " 1 ?\n";
  os << '\n';

  os << "loop_\n"
        "_diffrn_measurement.diffrn_id\n"
        "_diffrn_measurement.details\n";
  for (const XdsAscii::Iset& iset : xds.isets) {
    os << iset.id;
    if (iset.frame_count >= 0)
      os << " '" << iset.frame_count << " frames'\n";
    else
      os << " ?\n";
  }
  os << '\n';

  double w_all = wavelength;
  if (std::isnan(w_all) &&
      std::all_of(xds.isets.begin(), xds.isets.end(),
                  [&](const XdsAscii::Iset& i) { return i.wavelength == xds.wavelength; }))
    w_all = xds.wavelength;

  os << "loop_\n"
        "_diffrn_radiation.diffrn_id\n"
        "_diffrn_radiation.wavelength_id\n";
  for (const XdsAscii::Iset& iset : xds.isets)
    os << iset.id << ' ' << (std::isnan(w_all) ? iset.id : 1) << '\n';
  os << '\n';

  os << "loop_\n"
        "_diffrn_radiation_wavelength.id\n"
        "_diffrn_radiation_wavelength.wavelength\n";
  auto number_or_dot = [](double d) { return std::isnan(d) ? "." : to_str(d); };
  if (std::isnan(w_all)) {
    for (const XdsAscii::Iset& iset : xds.isets)
      os << iset.id << ' ' << number_or_dot(iset.wavelength) << '\n';
  } else {
    os << '1' << ' ' << number_or_dot(w_all) << '\n';
  }
  os << '\n';

  if (enable_UB && xds.has_cell_axes()) {
    os << "loop_\n"
          "_diffrn_orient_matrix.diffrn_id\n"
          "_diffrn_orient_matrix.UB[1][1]\n"
          "_diffrn_orient_matrix.UB[1][2]\n"
          "_diffrn_orient_matrix.UB[1][3]\n"
          "_diffrn_orient_matrix.UB[2][1]\n"
          "_diffrn_orient_matrix.UB[2][2]\n"
          "_diffrn_orient_matrix.UB[2][3]\n"
          "_diffrn_orient_matrix.UB[3][1]\n"
          "_diffrn_orient_matrix.UB[3][2]\n"
          "_diffrn_orient_matrix.UB[3][3]\n";
    Mat33 q = xds.calculate_conversion_from_cambridge().inverse();
    Mat33 ub = q.multiply(xds.cell_axes.inverse());
    // AFAICS if we have geometry information we don't have ISET records,
    // so we'll have only one row here.
    for (const XdsAscii::Iset& iset : xds.isets) {
      WRITE("%d  %#g %#g %#g  %#g %#g %#g  %#g %#g %#g\n", iset.id,
            ub.a[0][0], ub.a[0][1], ub.a[0][2],
            ub.a[1][0], ub.a[1][1], ub.a[1][2],
            ub.a[2][0], ub.a[2][1], ub.a[2][2]);
    }
    os << '\n';
  }

  const SpaceGroup* sg = find_spacegroup_by_number(xds.spacegroup_number);
  double rmsds[6] = {0., 0., 0., 0., 0., 0.};
  if (xds.isets.size() > 1) {
    double mean[6] = {xds.unit_cell.a, xds.unit_cell.b, xds.unit_cell.c,
                      xds.unit_cell.alpha, xds.unit_cell.beta, xds.unit_cell.gamma};
    int n = 0;
    for (const XdsAscii::Iset& iset : xds.isets) {
      for (int i = 0; i < 6; ++i)
        rmsds[i] += sq(mean[i] - iset.cell_constants[i]) * iset.frame_count;
      n += iset.frame_count;
    }
    for (int i = 0; i < 6; ++i)
      rmsds[i] = std::sqrt(rmsds[i] / n);
  }
  write_cell_and_symmetry(entry_id, xds.unit_cell, rmsds, sg, buf, os);

  os << "\nloop_"
        "\n_diffrn_refln.diffrn_id"
        "\n_diffrn_refln.id"
        "\n_diffrn_refln.index_h"
        "\n_diffrn_refln.index_k"
        "\n_diffrn_refln.index_l"
        "\n_diffrn_refln.intensity_net"
        "\n_diffrn_refln.intensity_sigma";
  if (xds.oscillation_range != 0.)
    os << "\n_diffrn_refln.pdbx_scan_angle";
  os << "\n_diffrn_refln.pdbx_image_id\n";
  int idx = 0;
  for (const XdsAscii::Refl& refl : xds.data) {
    if (refl.sigma < 0 && skip_negative_sigi)  // misfit
      continue;
    char* ptr = buf;
    ptr += snprintf_z(ptr, 128, "%d %d %d %d %d %g %.5g ",
                      refl.iset, ++idx, refl.hkl[0], refl.hkl[1], refl.hkl[2],
                      refl.iobs, refl.sigma);
    if (xds.oscillation_range != 0.) {
      double angle = xds.rot_angle(refl);
      ptr += snprintf_z(ptr, 16, "%.5g ", angle);
    }
    ptr += snprintf_z(ptr, 16, "%d\n", refl.frame());
    os.write(buf, ptr - buf);
  }
}

#undef WRITE

void remove_appendix_from_column_names(Mtz& mtz, std::ostream& out) {
  std::string appendix;
  for (char type : {'J', 'F'}) {  // intensity, amplitude
    auto cols = mtz.columns_with_type(type);
    // Remove the appendix only in clear cases. If mtz has multiple
    // intensity columns they could have IMEAN_this and IMEAN.
    if (cols.size() == 1) {
      auto pos = cols[0]->label.find('_');
      if (pos == std::string::npos)
        return;
      appendix = cols[0]->label.substr(pos);
      break;
    }
  }
  if (appendix.empty())
    return;
  out << "Ignoring '" << appendix << "' appended to column names.\n";
  for (Mtz::Column& col : mtz.columns) {
    size_t from_end  = appendix.size();
    // the appendix can before (+)/(-), i.e. I_appendix(+)
    if (!col.label.empty() && col.label.back() == ')')
      from_end += 3;
    if (from_end < col.label.size() &&
        col.label.compare(col.label.size() - from_end, appendix.size(), appendix) == 0)
      col.label.erase(col.label.size() - from_end, appendix.size());
  }
}

bool validate_merged_mtz_deposition_columns(const Mtz& mtz, std::ostream& out) {
  bool ok = true;
  if (!mtz.rfree_column()) {
    out << "ERROR. Merged file is missing free-set flag.\n";
    ok = false;
  }
  if (!mtz.imean_column() && !mtz.iplus_column()) {
    out << "ERROR. Merged file is missing intensities.\n";
    ok = false;
  }
  if (!mtz.column_with_one_of_labels({"F", "FP", "FOBS", "F-obs",
                                      "F(+)", "FOBS(+)", "F-obs(+)"})) {
    out << "Merged file is missing amplitudes\n"
           "(which is fine if intensities were used for refinement)\n";
  }
  if (!ok) {
    out << "Columns in the merged file:";
    for (const Mtz::Column& col : mtz.columns)
      out << ' ' << col.label;
    out << '\n';
  }
  return ok;
}

bool validate_merged_intensities(Intensities& mi, Intensities& ui,
                                 bool relaxed_check, std::ostream& out) {
  // XDS files have 4 significant digits. Using accuracy 5x the precision.
  const double max_diff = 0.005;
  out << "Checking if both files match...\n";
  bool ok = true;
  if (ui.spacegroup == mi.spacegroup) {
    out << "The same space group: " << mi.spacegroup_str() << '\n';
  } else {
    GroupOps gops1 = ui.spacegroup->operations();
    GroupOps gops2 = mi.spacegroup->operations();
    if (!gops1.has_same_centring(gops2) || !gops1.has_same_rotations(gops2))
      ok = false;
    out << (ok ? "WARNING" : "ERROR")
        << ". Different space groups in merged and unmerged files:\n"
        << mi.spacegroup_str() << " and " << ui.spacegroup_str() << '\n';
    if (!ok)
      out << "(in the future, this app may recognize compatible space groups\n"
             "and reindex unmerged data if needed; for now, it's on you)\n";
  }

  auto eq = [](double x, double y, double rmsd) { return std::fabs(x - y) < rmsd + 0.02; };
  if(eq(mi.unit_cell.a,     ui.unit_cell.a,     ui.unit_cell_rmsd[0]) &&
     eq(mi.unit_cell.b,     ui.unit_cell.b,     ui.unit_cell_rmsd[1]) &&
     eq(mi.unit_cell.c,     ui.unit_cell.c,     ui.unit_cell_rmsd[2]) &&
     eq(mi.unit_cell.alpha, ui.unit_cell.alpha, ui.unit_cell_rmsd[3]) &&
     eq(mi.unit_cell.beta,  ui.unit_cell.beta,  ui.unit_cell_rmsd[4]) &&
     eq(mi.unit_cell.gamma, ui.unit_cell.gamma, ui.unit_cell_rmsd[5])) {
    out << "The same unit cell parameters.\n";
  } else {
    const UnitCell& mc = mi.unit_cell;
    const UnitCell& uc = ui.unit_cell;
    out << "Unit cell parameters differ:";
    out << "\n    merged: " << mc.a << ' ' << mc.b << ' ' << mc.c << "  "
                            << mc.alpha << ' ' << mc.beta << ' ' << mc.gamma;
    out << "\n  unmerged: " << uc.a << ' ' << uc.b << ' ' << uc.c << "  "
                            << uc.alpha << ' ' << uc.beta << ' ' << uc.gamma;
    out << '\n';
    ok = false;
  }

  size_t ui_size1 = ui.data.size();
  ui.merge_in_place(mi.type);  // it also sorts
  size_t ui_size2 = ui.data.size();
  ui.remove_systematic_absences();
  out << "Unmerged reflections: " << ui_size1 << " (" << ui_size2 << " merged "
      << mi.type_str() << ", " << ui.data.size() << " w/o sysabs)\n";
  mi.switch_to_asu_indices(/*merged=*/true);
  mi.sort();
  size_t mi_size1 = mi.data.size();
  mi.remove_systematic_absences();
  out << "Merged reflections: " << mi_size1 << ' ' << mi.type_str()
      << " (" << mi.data.size() << " w/o sysabs)\n";

  if (mi.staraniso_b.ok()) {
    out << "Taking into account the anisotropy tensor that was used for scaling.\n";
    for (Intensities::Refl& refl : ui.data)
      refl.value *= mi.staraniso_b.scale(refl.hkl, ui.unit_cell);
  }

  // first pass - calculate CC and scale
  gemmi::Correlation corr = calculate_hkl_value_correlation(ui.data, mi.data);
  out << "Corr. coef. of " << corr.n << ' ' << mi.type_str() << " values: "
      << 100 * corr.coefficient() << "%\n";
  double scale = corr.mean_ratio();
  out << "Ratio of compared intensities (merged : unmerged): " << scale << '\n';

  // second pass - check that all reflections match
  double max_weighted_sq_diff = 0.;
  const Intensities::Refl* max_diff_r1 = nullptr;
  const Intensities::Refl* max_diff_r2 = nullptr;
  int differ_count = 0;
  int missing_count = 0;
  auto r1 = ui.data.begin();
  auto r2 = mi.data.begin();
  auto refln_str = [](const Intensities::Refl& r) {
    return r.intensity_label() + std::string(" ") + miller_str(r.hkl);
  };
  while (r1 != ui.data.end() && r2 != mi.data.end()) {
    if (r1->hkl == r2->hkl && r1->isign == r2->isign) {
      if (!relaxed_check) {
        double value1 = scale * r1->value;
        double sigma1 = scale * r1->sigma; // is this approximately correct
        double sq_max = std::max(sq(value1), sq(r2->value));
        double sq_diff = sq(value1 - r2->value);
        // Intensities may happen to be rounded to two decimal places,
        // so if the absolute difference is <0.01 it's OK.
        if (sq_diff > 1e-4 && sq_diff > sq(max_diff) * sq_max) {
          if (differ_count == 0) {
            out << "First difference: " << r1->hkl_label()
                << ' ' << value1 << " vs " << r2->value << '\n';
          }
          ++differ_count;
          double weighted_sq_diff = sq_diff / (sq(sigma1) + sq(r2->sigma));
          if (weighted_sq_diff > max_weighted_sq_diff) {
            max_weighted_sq_diff = weighted_sq_diff;
            max_diff_r1 = &*r1;
            max_diff_r2 = &*r2;
          }
        }
      }
      ++r1;
      ++r2;
    } else if (*r1 < *r2) {
      ++r1;
    } else {
      if (missing_count == 0)
        out << "First missing reflection in unmerged data: " << refln_str(*r1) << '\n';
      ++missing_count;
      ++r2;
    }
  }

  if (differ_count != 0) {
    out << "Most significant difference: " << refln_str(*max_diff_r1) << ' '
        << scale * max_diff_r1->value << " vs " << max_diff_r2->value << '\n';
    out << differ_count << " of " << corr.n << " intensities differ too much (by >"
        << to_str(max_diff * 100) << "%).\n";
    if (differ_count >= 0.001 * corr.n)
      ok = false;
    else
      out << "(less than 0.1% of all intensities -"
          << " probably outlier rejection during merging)\n";
  }
  if (missing_count != 0) {
    out << missing_count << " out of " << mi.data.size()
        << " reflections in the merged file not found in unmerged data\n";
    ok = false;
  }
  if (!relaxed_check) {
    if (differ_count == 0 && missing_count == 0) {
      out << "Intensities match.";
      if (!ok)
        out << " But other problems were found (see above).";
      out << '\n';
    } else {
      out << (ok ? "OK. Intensities almost match.\n"
                 : "ERROR. Intensities do not match.\n");
    }
  }
  return ok;
}

} // namespace gemmi
