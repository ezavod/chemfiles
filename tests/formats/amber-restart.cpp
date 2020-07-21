// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include "catch.hpp"
#include "chemfiles.hpp"
#include "helpers.hpp"
using namespace chemfiles;

TEST_CASE("Read files in Amber Restart format") {
    SECTION("Water") {
        auto file = Trajectory("data/netcdf/water.ncrst");
        CHECK(file.nsteps() == 1);
        auto frame = file.read();
        CHECK(frame.size() == 297);

        // Check cell
        auto cell = frame.cell();
        CHECK(cell.shape() == UnitCell::ORTHORHOMBIC);
        CHECK(approx_eq(cell.a(), 15.0, 1e-5));
        CHECK(approx_eq(cell.b(), 15.0, 1e-5));
        CHECK(approx_eq(cell.c(), 15.0, 1e-5));

        // Check positions
        auto positions = frame.positions();
        CHECK(approx_eq(positions[0], Vector3D(0.4172191, 8.303366, 11.73717), 1e-4));
        CHECK(approx_eq(positions[296], Vector3D(6.664049, 11.61418, 12.96149), 1e-4));
    }

    SECTION("Missing unit cell") {
        auto file = Trajectory("data/netcdf/no-cell.ncrst");
        // Check `read_step`
        auto frame = file.read_step(0);
        CHECK(frame.size() == 1989);
        CHECK(frame.cell() == UnitCell());
    }

    SECTION("Scale factor") {
        auto file = Trajectory("data/netcdf/scaled_traj.ncrst");
        auto frame = file.read();
        CHECK(frame.size() == 1938);

        // Check cell
        auto cell = frame.cell();
        CHECK(cell.shape() == UnitCell::ORTHORHOMBIC);
        CHECK(approx_eq(cell.a(), 60.9682 * 1.765, 1e-4));
        CHECK(approx_eq(cell.b(), 60.9682 * 1.765, 1e-4));
        CHECK(cell.c() == 0);

        // Check positions
        auto positions = frame.positions();
        CHECK(approx_eq(positions[0], Vector3D(1.39, 1.39, 0) * 0.455, 1e-4));
        CHECK(approx_eq(positions[296], Vector3D(29.10, 37.41, 0) * 0.455, 1e-4));

        // Check velocities
        auto velocities = *frame.velocities();
        CHECK(
            approx_eq(velocities[1400], Vector3D(-0.042603, -0.146347, 12.803150) * -0.856, 1e-4));
        CHECK(approx_eq(velocities[1600], Vector3D(0.002168, 0.125240, 4.188500) * -0.856, 1e-4));
    }
}

TEST_CASE("Write files in Amber Restart format") {
    auto tmpfile = NamedTempPath(".ncrst");

    auto file = Trajectory(tmpfile, 'w');
    Frame frame;
    frame.resize(4);
    auto positions = frame.positions();
    frame.add_velocities();
    auto velocities = *frame.velocities();
    for (size_t i = 0; i < 4; i++) {
        positions[i] = Vector3D(1, 2, 3);
        velocities[i] = Vector3D(-3, -2, -1);
    }
    file.write(frame);

    CHECK_THROWS_WITH(file.write(frame), "AMBER Restart format only supports writing one frame");

    file.close();

    Trajectory check(tmpfile, 'r');
    frame = check.read();

    positions = frame.positions();
    CHECK(approx_eq(positions[0], Vector3D(1, 2, 3), 1e-4));
    CHECK(approx_eq(positions[1], Vector3D(1, 2, 3), 1e-4));
    CHECK(approx_eq(positions[2], Vector3D(1, 2, 3), 1e-4));
    CHECK(approx_eq(positions[3], Vector3D(1, 2, 3), 1e-4));

    CHECK(frame.velocities());
    velocities = *frame.velocities();
    CHECK(approx_eq(velocities[0], Vector3D(-3, -2, -1), 1e-4));
    CHECK(approx_eq(velocities[1], Vector3D(-3, -2, -1), 1e-4));
    CHECK(approx_eq(velocities[2], Vector3D(-3, -2, -1), 1e-4));
    CHECK(approx_eq(velocities[3], Vector3D(-3, -2, -1), 1e-4));
}
