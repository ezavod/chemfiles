// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cstdint>
#include <cstdlib>

#include <string>
#include <vector>
#include <functional>

#include "chemfiles/config.h"
#include "chemfiles/capi/types.h"
#include "chemfiles/capi/misc.h"
#include "chemfiles/capi/utils.hpp"
#include "chemfiles/capi/shared_allocator.hpp"

#include "chemfiles/misc.hpp"
#include "chemfiles/external/optional.hpp"

#include "chemfiles/FormatMetadata.hpp"

using namespace chemfiles;

static_assert(sizeof(chfl_status) == sizeof(int), "Wrong size for chfl_status enum");

static CHFL_THREAD_LOCAL std::string CAPI_LAST_ERROR;

void chemfiles::set_last_error(const std::string& message) {
    CAPI_LAST_ERROR = message;
}

extern "C" void chfl_free(const void* const object) {
    shared_allocator::free(object); // NOLINT: we can use shared_allocator::free to implement chfl_free
}

extern "C" const char* chfl_version(void) {
    return CHEMFILES_VERSION;
}

extern "C" const char* chfl_last_error(void) {
    return CAPI_LAST_ERROR.c_str();
}

extern "C" chfl_status chfl_clear_errors(void) {
    CHFL_ERROR_CATCH(
        CAPI_LAST_ERROR = "";
    )
}

extern "C" chfl_status chfl_set_warning_callback(chfl_warning_callback callback) {
    CHFL_ERROR_CATCH(
        set_warning_callback([callback](std::string message) {
            callback(message.c_str());
        });
    )
}

extern "C" chfl_status chfl_add_configuration(const char* path) {
    CHFL_ERROR_CATCH(
        add_configuration(path);
    )
}

extern "C" chfl_status chfl_formats_list(chfl_format_metadata** metadata, uint64_t* count) {
    CHECK_POINTER(metadata);
    CHECK_POINTER(count);
    CHFL_ERROR_CATCH(
        auto formats = formats_list();
        *count = static_cast<uint64_t>(formats.size());
        *metadata = static_cast<chfl_format_metadata*>(malloc(*count * sizeof(chfl_format_metadata)));  // NOLINT: use malloc instead of chfl_allocator to allow users to call free
        for (uint64_t i=0; i<*count; i++) {
            // here we rely on the fact that only one instance of each matdata
            // exist, and that they come from static storage, allowing us to
            // take reference to data directly inside `formats[i]` even if is it
            // a local temporary.
            const auto& meta = formats[i].get();
            (*metadata)[i].name = meta.name;
            (*metadata)[i].extension = meta.extension.value_or(nullptr);
            (*metadata)[i].description = meta.description;
            (*metadata)[i].reference = meta.reference;

            (*metadata)[i].read = meta.read;
            (*metadata)[i].write = meta.write;
            (*metadata)[i].memory = meta.memory;

            (*metadata)[i].positions = meta.positions;
            (*metadata)[i].velocities = meta.velocities;
            (*metadata)[i].unit_cell = meta.unit_cell;
            (*metadata)[i].atoms = meta.atoms;
            (*metadata)[i].bonds = meta.bonds;
            (*metadata)[i].residues = meta.residues;
        }
    )
}
