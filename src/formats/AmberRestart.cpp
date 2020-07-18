// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <array>
#include <cassert>
#include <string>
#include <vector>

#include "chemfiles/File.hpp"
#include "chemfiles/Format.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/UnitCell.hpp"
#include "chemfiles/files/NcFile.hpp"

#include "chemfiles/config.h"
#include "chemfiles/error_fmt.hpp"
#include "chemfiles/external/optional.hpp"
#include "chemfiles/external/span.hpp"
#include "chemfiles/types.hpp"
#include "chemfiles/warnings.hpp"

#include "chemfiles/formats/AmberRestart.hpp"

using namespace chemfiles;

template <> FormatInfo chemfiles::format_information<AmberRestartFormat>() {
    return FormatInfo("Amber Restart")
        .with_extension(".ncrst")
        .description("Amber convention for binary NetCDF Restart files");
}

//! Check the validity of a NetCDF file
static bool is_valid(const NcFile& file_, size_t natoms) {
    bool writing = (natoms != static_cast<size_t>(-1));

    if (file_.global_attribute("Conventions") != "AMBERRESTART") {
        if (!writing) {
            warning("Amber Restart reader", "we can only read AMBER convention");
        }
        return false;
    }

    if (file_.global_attribute("ConventionVersion") != "1.0") {
        if (!writing) {
            warning("Amber Restart reader", "we can only read version 1.0 of AMBER convention");
        }
        return false;
    }

    if (file_.dimension("spatial") != 3) {
        if (!writing) {
            warning("Amber Restart reader", "wrong size for spatial dimension: should be 3, is {}",
                    file_.dimension("spatial"));
        }
        return false;
    }

    if (writing) {
        if (file_.dimension("atom") != natoms) {
            warning("Amber Restart writer", "wrong size for atoms dimension: should be {}, is {}",
                    natoms, file_.dimension("atom"));
            return false;
        }
    }
    return true;
}

AmberRestartFormat::AmberRestartFormat(std::string path, File::Mode mode,
                                       File::Compression compression)
    : file_(std::move(path), mode), step_done_(false), validated_(false) {
    if (file_.mode() == File::READ) {
        if (!is_valid(file_, static_cast<size_t>(-1))) {
            throw format_error("invalid AMBER Restart file at '{}'", file_.path());
        }
        validated_ = true;
    } else if (file_.mode() == File::APPEND) {
        throw format_error("append mode ('a') is not supported with AMBER Restart format");
    }
    if (compression != File::DEFAULT) {
        throw format_error("compression is not supported with NetCDF format");
    }
}

size_t AmberRestartFormat::nsteps() { return 1; }

void AmberRestartFormat::read_step(const size_t step, Frame& frame) {
    if (step != 0) {
        throw format_error("AMBER Restart format only supports reading one frame");
    }

    frame.set_cell(read_cell());

    frame.resize(file_.dimension("atom"));
    read_array(frame.positions(), "coordinates");
    if (file_.variable_exists("velocities")) {
        frame.add_velocities();
        read_array(*frame.velocities(), "velocities");
    }

    step_done_ = true;
}

void AmberRestartFormat::read(Frame& frame) {
    if (step_done_) {
        throw format_error("AMBER Restart format only supports reading one frame");
    }

    read_step(0, frame);
}

UnitCell AmberRestartFormat::read_cell() {
    if (!file_.variable_exists("cell_lengths") || !file_.variable_exists("cell_angles")) {
        return {}; // No UnitCell information
    }

    if (file_.optional_dimension("cell_spatial", 0) != 3 ||
        file_.optional_dimension("cell_angular", 0) != 3) {
        return {}; // No UnitCell information
    }

    auto length_var = file_.variable<nc::NcDouble>("cell_lengths");
    auto angles_var = file_.variable<nc::NcDouble>("cell_angles");

    std::vector<size_t> start{0};
    std::vector<size_t> count{3};

    auto length = length_var.get(start, count);
    auto angles = angles_var.get(start, count);

    assert(length.size() == 3);
    assert(angles.size() == 3);

    if (length_var.attribute_exists("scale_factor")) {
        float scale_factor = length_var.float_attribute("scale_factor");
        length[0] *= scale_factor;
        length[1] *= scale_factor;
        length[2] *= scale_factor;
    }

    if (angles_var.attribute_exists("scale_factor")) {
        float scale_factor = angles_var.float_attribute("scale_factor");
        angles[0] *= scale_factor;
        angles[1] *= scale_factor;
        angles[2] *= scale_factor;
    }

    return {length[0], length[1], length[2], angles[0], angles[1], angles[2]};
}

void AmberRestartFormat::read_array(span<Vector3D> array, const std::string& name) {
    auto array_var = file_.variable<nc::NcDouble>(name);
    auto natoms = file_.dimension("atom");
    assert(array.size() == natoms);

    std::vector<size_t> start{0, 0};
    std::vector<size_t> count{natoms, 3};
    auto data = array_var.get(start, count);

    if (array_var.attribute_exists("scale_factor")) {
        float scale_factor = array_var.float_attribute("scale_factor");
        for (auto& value : data) {
            value *= scale_factor;
        }
    }

    for (size_t i = 0; i < natoms; i++) {
        array[i][0] = data[3 * i + 0];
        array[i][1] = data[3 * i + 1];
        array[i][2] = data[3 * i + 2];
    }
}

// Initialize a file, assuming that it is empty
static void initialize(NcFile& file, size_t natoms, bool with_velocities) {
    file.set_nc_mode(NcFile::DEFINE);

    file.add_global_attribute("Conventions", "AMBERRESTART");
    file.add_global_attribute("ConventionVersion", "1.0");
    file.add_global_attribute("program", "Chemfiles");
    file.add_global_attribute("programVersion", CHEMFILES_VERSION);

    file.add_dimension("spatial", 3);
    file.add_dimension("atom", natoms);
    file.add_dimension("cell_spatial", 3);
    file.add_dimension("cell_angular", 3);
    file.add_dimension("label", nc::STRING_MAXLEN);

    auto spatial = file.add_variable<nc::NcChar>("spatial", "spatial");
    auto cell_spatial = file.add_variable<nc::NcChar>("cell_spatial", "cell_spatial");
    auto cell_angular = file.add_variable<nc::NcChar>("cell_angular", "cell_angular", "label");

    auto coordinates = file.add_variable<nc::NcDouble>("coordinates", "atom", "spatial");
    coordinates.add_string_attribute("units", "angstrom");

    auto cell_lenght = file.add_variable<nc::NcDouble>("cell_lengths", "cell_spatial");
    cell_lenght.add_string_attribute("units", "angstrom");

    auto cell_angles = file.add_variable<nc::NcDouble>("cell_angles", "cell_angular");
    cell_angles.add_string_attribute("units", "degree");

    if (with_velocities) {
        auto velocities = file.add_variable<nc::NcDouble>("velocities", "atom", "spatial");
        velocities.add_string_attribute("units", "angstrom/picosecond");
    }
    file.set_nc_mode(NcFile::DATA);

    spatial.add("xyz");
    cell_spatial.add("abc");
    cell_angular.add({"alpha", "beta", "gamma"});
}

void AmberRestartFormat::write(const Frame& frame) {
    if (step_done_) {
        throw format_error("AMBER Restart format only supports writing one frame");
    }

    auto natoms = frame.size();
    // If we created the file, let's initialize it.
    if (!validated_) {
        initialize(file_, natoms, bool(frame.velocities()));
        assert(is_valid(file_, natoms));
        validated_ = true;
    }
    write_cell(frame.cell());
    write_array(frame.positions(), "coordinates");
    auto velocities = frame.velocities();
    if (velocities) {
        write_array(*velocities, "velocities");
    }

    step_done_ = true;
}

void AmberRestartFormat::write_array(const std::vector<Vector3D>& array, const std::string& name) {
    auto var = file_.variable<nc::NcDouble>(name);
    auto natoms = array.size();
    std::vector<size_t> start{0, 0};
    std::vector<size_t> count{natoms, 3};

    auto data = std::vector<double>(natoms * 3);
    for (size_t i = 0; i < natoms; i++) {
        data[3 * i + 0] = array[i][0];
        data[3 * i + 1] = array[i][1];
        data[3 * i + 2] = array[i][2];
    }
    var.add(start, count, data);
}

void AmberRestartFormat::write_cell(const UnitCell& cell) {
    auto length = file_.variable<nc::NcDouble>("cell_lengths");
    auto angles = file_.variable<nc::NcDouble>("cell_angles");

    auto length_data = std::vector<double>{cell.a(), cell.b(), cell.c()};

    auto angles_data = std::vector<double>{cell.alpha(), cell.beta(), cell.gamma()};

    std::vector<size_t> start{0};
    std::vector<size_t> count{3};
    length.add(start, count, length_data);
    angles.add(start, count, angles_data);
}
