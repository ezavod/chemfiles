// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_FORMAT_NCRST_HPP
#define CHEMFILES_FORMAT_NCRST_HPP

#include <string>
#include <vector>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"
#include "chemfiles/files/NcFile.hpp"

namespace chemfiles {

class Frame;
class UnitCell;
class Vector3D;
template <class T> class span;

/// [Amber Restart][NetCDF] file format reader.
///
/// [NetCDF]: http://ambermd.org/netcdf/nctraj.xhtml
class AmberRestartFormat final : public Format {
  public:
    AmberRestartFormat(std::string path, File::Mode mode, File::Compression compression);

    void read_step(size_t step, Frame& frame) override;
    void read(Frame& frame) override;
    void write(const Frame& frame) override;

    size_t nsteps() override;

  private:
    /// Read the unit cell, the file is assumed to be valid.
    UnitCell read_cell();
    /// Generic function to read an std::vector<Vector3D>, the file is assumed to be valid.
    void read_array(span<Vector3D> array, const std::string& name);

    /// Write an std::vector<Vector3D> to the file, as a variable with the name `name`.
    void write_array(const std::vector<Vector3D>& array, const std::string& name);
    /// Write an UnitCell to the file.
    void write_cell(const UnitCell& cell);

    /// Associated NetCDF file.
    NcFile file_;
    /// Has the step been read?
    bool step_done_;
    /// Was the associated file validated?
    bool validated_;
};

template <> FormatInfo format_information<AmberRestartFormat>();

} // namespace chemfiles

#endif
