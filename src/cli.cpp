#include "cli.hpp"
#include "configuration.hpp"
#include "build_metadata.hpp"
#include "sdbd/export.hpp"
#include <algorithm>
#include <chrono>
#include <locale>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string_view>

namespace sdbd::cli {
namespace {
// Print commands or focused help for a known operation; no configuration is loaded.
void help(std::string_view command = {}) {
    std::cout<<"SALTED — steepest-descent ballistic deposition\n";
    if (command.empty()) {
        std::cout<<"  deposit run CONFIG [--output DIRECTORY]       Deposit particles\n"
            <<"  deposit generate CONFIG [--output FILE]       Prepare STL collision geometry\n"
            <<"  deposit validate --config CONFIG              Check inputs without running\n\n"
            <<"Start: deposit run examples/cube/config.conf\n"
            <<"View deposition.vtk in ParaView with 3D Glyphs (Sphere).\n"
            <<"Use deposit COMMAND --help for settings and advanced overrides.\n"
            <<"Other actions: --help, --version, --print-example\n";
        return;
    }
    if (command=="run") {
        std::cout<<"Usage: deposit run CONFIG [--output DIRECTORY] [--overrides FILE]\n"
            <<"Settings: particle = sphere|spiky|cube|stick|FILE, count, seed, box, units.\n"
            <<"Use size for built-ins; length, radius and spheres for a stick.\n"
            <<"Defaults: sphere, count 24, seed 84328, box 720 720 6680, units mm.\n"
            <<"Outputs: run.conf, deposition.vtk, plus stl_poses.vtk when applicable.\n"
            <<"Sides are periodic; box height limits accepted surface height, without a lid.\n";
    } else if (command=="generate") {
        std::cout<<"Usage: deposit generate CONFIG [--output FILE] [--overrides FILE]\n"
            <<"Required: stl (binary mesh), declared units, and output particle file.\n"
            <<"Controls: scale, voxel_resolution, maximum_spheres, accuracy.\n"
            <<"Requires a GEMSS-enabled build. Generate once, then run with particle = FILE.\n";
    } else {
        std::cout<<"Usage: deposit validate --config CONFIG [--overrides FILE]\n"
            <<"The positional form, deposit validate CONFIG, is also accepted.\n"
            <<"Checks run inputs or STL options when stl is present; output is optional.\n"
            <<"No deposition, voxelization, reconstruction or output files.\n";
    }
    std::cout<<"Configuration: key = value; # comments; quote paths containing spaces.\n"
        <<"Precedence: defaults < CONFIG < --overrides FILE < --output.\n"
        <<"Config paths are file-relative; --output is cwd-relative. Existing outputs are refused.\n"
        <<"See README.md and docs/configuration.md for advanced controls and limits.\n";
}
// Append concise replay-ignored metadata to the one run record; throw on failed I/O.
void metadata(const std::filesystem::path& record,const Configuration& c,
              const std::optional<ParticleGeometry>& geometry,const std::string& state,double seconds) {
    std::ofstream out(record,std::ios::app);
    out.imbue(std::locale::classic());
    out<<std::setprecision(std::numeric_limits<long double>::max_digits10);
    out<<"# status "<<std::quoted(state)<<"\n# SALTED "<<SDBD_VERSION<<"\n# compiler "<<build::compiler
       <<"\n# GSL "<<SDBD_GSL_VERSION<<"\n# source_sha256 "<<build::source_sha256
       <<"\n# build_configuration "<<build::configuration<<"\n# compile_flags "<<std::quoted(build::flags)
       <<"\n# long_double_digits "<<std::numeric_limits<long double>::digits
       <<"\n# wall_seconds "<<seconds;
    const auto& g=c.simulation.grid;
    out<<"\n# grid "<<g.nx<<' '<<g.ny<<' '<<g.nz<<' '<<g.dx<<' '<<g.zo<<'\n';
    if (geometry && geometry->physical_properties()) {
        const auto& p=*geometry->physical_properties();
        out<<"# target_volume "<<p.volume<<"\n# density "<<p.density<<"\n# target_com "<<geometry->center_of_mass()
           <<"\n# inertia_per_density";
        for (auto v:p.inertia_per_density) out<<' '<<v;
        out<<"\n# backend "<<std::quoted(p.backend_version)<<"\n# generation_parameters "<<std::quoted(p.parameters)
           <<"\n# source_checksum "<<std::quoted(p.checksum)<<"\n# source_units "<<std::quoted(p.source_units)
           <<"\n# source_scale "<<p.source_scale<<"\n# voxel_size "<<p.voxel_size<<'\n';
    }
    if (geometry && geometry->stl_source()) {
        auto source=geometry->stl_source()->path.lexically_relative(std::filesystem::absolute(record).parent_path());
        if (source.empty()) source=geometry->stl_source()->path;
        out<<"# source_stl "<<std::quoted(source.generic_string())<<"\n# quaternion_order WXYZ\n# transform P+R*(scale*v)\n";
    }
    out.close();
    if (!out) throw std::runtime_error("Cannot append run metadata");
}
// Validate a new destination without replacing any previous simulation or particle file.
void new_destination(const std::filesystem::path& path) {
    if (path.filename().empty()) throw std::invalid_argument("Output must name a new file or directory");
    if (std::filesystem::exists(path)) throw std::invalid_argument("Output already exists: "+path.string()+"; choose a new output");
}
// Describe resolved geometry and exact physical domain without running deposition.
void summary(const Configuration& c,const std::optional<ParticleGeometry>& geometry) {
    std::cout<<"Particle: "<<c.particle<<"; count: "<<c.simulation.count<<"; seed: "<<c.simulation.seed
        <<"; RNG: "<<c.simulation.rng<<"\nBox: "<<c.box<<' '<<c.units<<" (x/y periodic; open top)\n";
    if (geometry) std::cout<<"Spheres per particle: "<<geometry->spheres().size()<<"; physical COM: "<<geometry->center_of_mass()<<'\n';
}
}
int command_line(int argc,char** argv) {
    if (argc==2 && std::string_view(argv[1])=="--version") {
        std::cout<<"SALTED "<<SDBD_VERSION<<" (C++20; GEMSS "
#ifdef SDBD_CLI_GEMSS
                 <<"enabled";
#else
                 <<"disabled";
#endif
        std::cout<<")\n";
        return 0;
    }
    if (argc==2 && std::string_view(argv[1])=="--print-example") {
        std::cout<<"particle = sphere\nsize = 60\ncount = 24\nbox = 720 720 6680\nunits = mm\nseed = 84328\noutput = results/spheres\n";
        return 0;
    }
    if (argc==1 || (argc==2 && std::string_view(argv[1])=="--help")) {
        help();
        return argc==1 ? 1 : 0;
    }
    const std::string_view command=argv[1];
    if (command!="run" && command!="generate" && command!="validate") throw std::invalid_argument("Expected run, generate or validate; see --help");
    if (argc==3 && std::string_view(argv[2])=="--help") {
        help(command);
        return 0;
    }

    // Keep scientific settings in the configuration; CLI paths select inputs and output.
    std::optional<std::filesystem::path> config,output,overrides;
    for (int i=2;i<argc;++i) {
        const std::string_view arg=argv[i];
        if (arg=="--output" || arg=="--config" || arg=="--overrides") {
            if (++i>=argc || std::string_view(argv[i]).starts_with("--")) throw std::invalid_argument("Missing value for "+std::string(arg));
            auto& target=arg=="--output" ? output : arg=="--overrides" ? overrides : config;
            if (target) throw std::invalid_argument("Duplicate "+std::string(arg));
            target=argv[i];
        } else if (arg.starts_with("--")) throw std::invalid_argument("Set simulation/generation parameters in CONFIG, not "+std::string(arg));
        else if (config) throw std::invalid_argument("Supply one configuration file");
        else config=argv[i];
    }
    if (!config) throw std::invalid_argument("A configuration file is required; see --print-example");

    // Resolve precedence before validating either operation's inputs.
    auto fields=read_configuration(*config);
    if (overrides) {
        const auto replacements=read_configuration(*overrides);
        // A new geometry invalidates the previous file checksum; an explicitly
        // supplied replacement checksum is merged and checked normally.
        if (replacements.contains("particle")) fields.erase("particle_checksum");
        for (const auto& [key,value]:replacements) fields.insert_or_assign(key,value);
    }
    if (output) fields.insert_or_assign("output",Value{output->string(),std::filesystem::current_path(),"--output"});
    const bool validation=command=="validate",generation=command=="generate" || (validation && fields.contains("stl"));
    auto c=resolve(fields,generation,!validation);
    if (!validation) new_destination(c.output);

    // Prepare reusable collision geometry independently of deposition.
    if (generation) {
#ifdef SDBD_CLI_GEMSS
        const auto& g=*c.generation;
        if (validation) {
            validate_stl(g.stl,g.options);
            std::cout<<"Mesh/options valid. No voxelization, reconstruction or files written.\n";
            return 0;
        }
        const auto particle=generate_stl(g.stl,g.options);
        std::filesystem::create_directories(c.output.parent_path());
        save_particle(c.output,particle);
        const auto& p=*particle.physical_properties();
        const auto begin=p.parameters.find("achieved=");
        const auto achieved=begin==std::string::npos ? "unknown"
            : p.parameters.substr(begin+9,p.parameters.find(';',begin)-(begin+9));
        if (!c.quiet) std::cout<<"Generated "<<particle.spheres().size()<<" spheres in "<<c.output.string()
            <<"\nVoxel spacing: "<<p.voxel_size<<' '<<c.units<<"; target COM: "<<particle.center_of_mass()
            <<"\nRequested Dice: "<<g.options.reconstruction_accuracy<<"; achieved: "<<achieved
            <<"\nMass properties originate from the target grid (compute_physics=2).\n";
        if (achieved!="unknown" && std::stold(achieved)<g.options.reconstruction_accuracy)
            std::cerr<<"Reconstruction did not meet requested accuracy; inspect geometry before deposition.\n";
        return 0;
#else
        throw std::invalid_argument("GEMSS is disabled; generation/mesh validation needs SDBD_ENABLE_GEMSS=ON");
#endif
    }

    // Resolve the particle and physical domain before creating any simulation output.
    std::optional<ParticleGeometry> geometry;
    if (c.particle_file) {
        geometry=load_particle(*c.particle_file);
        const auto checksum=file_checksum(*c.particle_file);
        if (!c.particle_checksum.empty() && c.particle_checksum!=checksum) throw std::invalid_argument("Prepared particle checksum differs from saved run");
        c.particle_checksum=checksum;
    } else if (c.particle=="stick") geometry=generate_stick(c.length,c.radius,c.spheres);
    configure_domain(c,geometry);
    const bool poses=c.stl_poses && geometry && geometry->stl_source();
    if (poses) validate_stl_source(*geometry);
    if (!c.quiet) summary(c,geometry);
    if (validation) {
        std::cout<<"Validation successful. No simulation or files written.\n";
        return 0;
    }

    std::filesystem::create_directories(c.output.parent_path());
    if (!std::filesystem::create_directory(c.output)) throw std::invalid_argument("Output already exists");
    const auto record=c.output/"run.conf";
    write_configuration(record,c);
    const auto start=std::chrono::steady_clock::now();
    try {
        if (!c.quiet) c.simulation.progress=[](int n,int total) {
            if (n==total || n==1 || n%std::max(1,total/10)==0)
                std::cerr<<"Deposited "<<n<<'/'<<total<<'\n';
        };
        const auto legacy=[&](const char* name) {
            return c.legacy_output ? c.output/name : std::filesystem::path{};
        };
        simulate(c.simulation,{legacy("positions.txt"),legacy("height.txt"),legacy("details.txt"),legacy("particles.txt"),
            c.output/"deposition.vtk",poses ? c.output/"stl_poses.vtk" : std::filesystem::path{},record},geometry);
        metadata(record,c,geometry,"completed",std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
    } catch (const std::exception& e) {
        try {
            metadata(record,c,geometry,std::string("failed: ")+e.what(),
                     std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count());
        } catch (...) {}
        throw;
    }
    if (!c.quiet) std::cout<<"Output: "<<c.output.string()<<'\n';
    return 0;
}
}
