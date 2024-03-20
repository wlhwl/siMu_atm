/*
 * (c) Copyright 2018 CORSIKA Project, corsika-project@lists.kit.edu
 *
 * This software is distributed under the terms of the GNU General Public
 * Licence version 3 (GPL Version 3). See file LICENSE for a full version of
 * the license.
 */

/* clang-format off */
// InteractionCounter used boost/histogram, which
// fails if boost/type_traits have been included before. Thus, we have
// to include it first...
#include <corsika/framework/process/InteractionCounter.hpp>
/* clang-format on */
#include <corsika/framework/core/Cascade.hpp>
#include <corsika/framework/core/EnergyMomentumOperations.hpp>
#include <corsika/framework/core/Logging.hpp>
#include <corsika/framework/core/PhysicalUnits.hpp>
#include <corsika/framework/geometry/PhysicalGeometry.hpp>
#include <corsika/framework/geometry/Plane.hpp>
#include <corsika/framework/geometry/Sphere.hpp>
#include <corsika/framework/process/DynamicInteractionProcess.hpp>
#include <corsika/framework/process/ProcessSequence.hpp>
#include <corsika/framework/process/SwitchProcessSequence.hpp>
#include <corsika/framework/random/RNGManager.hpp>
#include <corsika/framework/utility/CorsikaFenv.hpp>

#include <corsika/modules/writers/EnergyLossWriter.hpp>
#include <corsika/modules/writers/LongitudinalWriter.hpp>
#include <corsika/modules/writers/PrimaryWriter.hpp>
#include <corsika/modules/writers/SubWriter.hpp>
#include <corsika/output/OutputManager.hpp>

#include <corsika/media/CORSIKA7Atmospheres.hpp>
#include <corsika/media/Environment.hpp>
#include <corsika/media/GeomagneticModel.hpp>
#include <corsika/media/GladstoneDaleRefractiveIndex.hpp>
#include <corsika/media/HomogeneousMedium.hpp>
#include <corsika/media/IMagneticFieldModel.hpp>
#include <corsika/media/LayeredSphericalAtmosphereBuilder.hpp>
#include <corsika/media/MediumPropertyModel.hpp>
#include <corsika/media/MediumProperties.hpp>
#include <corsika/media/NuclearComposition.hpp>
#include <corsika/media/ShowerAxis.hpp>
#include <corsika/media/UniformMagneticField.hpp>
#include <corsika/media/IMediumModel.hpp>
#include <corsika/media/IMediumPropertyModel.hpp>

#include <corsika/modules/BetheBlochPDG.hpp>
#include <corsika/modules/Epos.hpp>
#include <corsika/modules/ObservationPlane.hpp>
#include <corsika/modules/PROPOSAL.hpp>
#include <corsika/modules/ParticleCut.hpp>
#include <corsika/modules/Pythia8.hpp>
#include <corsika/modules/QGSJetII.hpp>
#include <corsika/modules/Sibyll.hpp>
#include <corsika/modules/Sophia.hpp>
#include <corsika/modules/StackInspector.hpp>
#include <corsika/modules/UrQMD.hpp>

#include <corsika/setup/SetupStack.hpp>
#include <corsika/setup/SetupTrajectory.hpp>

#include <boost/filesystem.hpp>

#include <CLI/App.hpp>
#include <CLI/Config.hpp>
#include <CLI/Formatter.hpp>

#include <cstdlib>
#include <iomanip>
#include <limits>
#include <string>

#include "TRandom3.h"
#include <TTree.h>
#include <TFile.h>
#include "TRandom.h"
#include "TF1.h"
#include "TMath.h"
#include "TString.h"
#include "TList.h"
#include "initial_energy_generator.hpp"

using namespace corsika;
using namespace std;

using IMediumType = IMediumPropertyModel<IMagneticFieldModel<IMediumModel>>;
using EnvType = Environment<IMediumType>;
template <typename TInterface>
using MyExtraEnv = MediumPropertyModel<UniformMagneticField<TInterface>>;
using StackType = setup::Stack<EnvType>;
using TrackingType = setup::Tracking;
using Particle = StackType::particle_type;

//
// This is the main example script which runs EAS with fairly standard settings
// w.r.t. what was implemented in CORSIKA 7. Users may want to change some of the
// specifics (observation altitude, magnetic field, energy cuts, etc.), but this
// example is the most physics-complete one and should be used for full simulations
// of particle cascades in air
//

template <typename TEnvironmentInterface, template <typename> typename TExtraEnv,
        typename TEnvironment, typename... TArgs>
void create_my_5layer_atmosphere(TEnvironment& env, AtmosphereId const atmId,
                                 Point const& center, TArgs... args) {

    // construct the atmosphere builder
    auto builder = make_layered_spherical_atmosphere_builder<
                   TEnvironmentInterface, TExtraEnv>::create(center, constants::EarthRadius::Mean,
                   std::forward<TArgs>(args)...);
    
    builder.setNuclearComposition({{Code::Hydrogen, Code::Oxygen}, {0.11, 0.89}});
    builder.addLinearLayer(1.02_g/square(1_cm), 1_cm, 0_km);                      // composition values from AIRES manual
    builder.setNuclearComposition({{Code::Nitrogen, Code::Oxygen, Code::Argon,},{0.7847, 1. - 0.7847 - 0.0047, 0.0047}});                                           
    auto const params = atmosphereParameterList[static_cast<uint8_t>(atmId)];
    for (int i = 0; i < 4; ++i) {
         builder.addExponentialLayer(params[i].offset, params[i].scaleHeight,
                                     params[i].altitude);
    }
    builder.addLinearLayer(params[4].offset, params[4].scaleHeight, params[4].altitude);
    builder.assemble(env);
}

long registerRandomStreams(long seed) {
  RNGManager<>::getInstance().registerRandomStream("cascade");
  RNGManager<>::getInstance().registerRandomStream("qgsjet");
  RNGManager<>::getInstance().registerRandomStream("sibyll");
  RNGManager<>::getInstance().registerRandomStream("sophia");
  RNGManager<>::getInstance().registerRandomStream("epos");
  RNGManager<>::getInstance().registerRandomStream("pythia");
  RNGManager<>::getInstance().registerRandomStream("urqmd");
  RNGManager<>::getInstance().registerRandomStream("proposal");
  if (seed == 0) {
    std::random_device rd;
    seed = rd();
    CORSIKA_LOG_INFO("random seed (auto) {}", seed);
  } else {
    CORSIKA_LOG_INFO("random seed {}", seed);
  }
  RNGManager<>::getInstance().setSeed(seed);
  return seed;
}

int main(int argc, char** argv) {

  // the main command line description
  CLI::App app{"Simulate standard (downgoing) showers with CORSIKA 8."};

  // some options that we want to fill in
  int A, Z, nevent = 0;

  // the following section adds the options to the parser

  // we start by definining a sub-group for the primary ID
  auto opt_Z = app.add_option("-Z", Z, "Atomic number for primary")
                   ->check(CLI::Range(0, 26))
                   ->group("Primary");
  auto opt_A = app.add_option("-A", A, "Atomic mass number for primary")
                   ->needs(opt_Z)
                   ->check(CLI::Range(1, 58))
                   ->group("Primary");
  app.add_option("-p,--pdg",
                 "PDG code for primary (p=2212, gamma=22, e-=11, nu_e=12, mu-=13, "
                 "nu_mu=14, tau=15, nu_tau=16).")
      ->excludes(opt_A)
      ->excludes(opt_Z)
      ->group("Primary");
  // the remainding options
  app.add_option("-E,--energy", "Primary energy in GeV")
      ->default_val(0.)
      ->group("Primary");
  app.add_option("--eMin", "Minimum energy of the inject spectrum [GeV]")
      ->check(CLI::Range(50.0, 1e8)) 
      ->default_val(1e5)
      ->group("Primary");
  app.add_option("--eMax", "Maximum energy of the inject spectrum [GeV]")
      ->check(CLI::Range(50.0, 1e10))
      ->default_val(1e8)
      ->group("Primary");
  app.add_option("-z,--zenith", "Primary zenith angle (deg)")
      ->default_val(0.)
      ->check(CLI::Range(0., 90.))
      ->group("Primary");
  app.add_option("-a,--azimuth", "Primary azimuth angle (deg)")
      ->default_val(0.)
      ->check(CLI::Range(0., 360.))
      ->group("Primary");
  app.add_option("--emcut",
                 "Min. kin. energy of photons, electrons and "
                 "positrons in tracking (GeV)")
      ->default_val(0.5e-3)
      ->check(CLI::Range(0.000001, 1.e13))
      ->group("Config");
  app.add_option("--hadcut", "Min. kin. energy of hadrons in tracking (GeV)")
      ->default_val(0.3)
      ->check(CLI::Range(0.000001, 1.e13))
      ->group("Config");
  app.add_option("--mucut", "Min. kin. energy of muons in tracking (GeV)")
      ->default_val(0.3)
      ->check(CLI::Range(0.000001, 1.e13))
      ->group("Config");
  app.add_option("--nucut", "Min. kin. energy of neutrinos in tracking (GeV)")
      ->default_val(0.3)
      ->check(CLI::Range(0.000001, 1.e13))
      ->group("Config");  
  bool track_neutrinos = false;
  app.add_flag("--track-neutrinos", track_neutrinos, "switch on tracking of neutrinos")
      ->group("Config");
  app.add_option("--neutrino-interaction-type",
                 "charged (CC) or neutral current (NC) or both")
      ->default_val("both")
      ->check(CLI::IsMember({"neutral", "NC", "charged", "CC", "both"}))
      ->group("Misc.");
  app.add_option("--observation-level",
                 "Height above earth radius of the observation level (in m)")
      ->default_val(0.)
      ->check(CLI::Range(-5.e3, 1.e5))
      ->group("Config");
  app.add_option("--injection-height",
                 "Height above earth radius of the injection point (in m)")
      ->default_val(112.75e3)
      ->check(CLI::Range(-1.e3, 1.e6))
      ->group("Config");
  app.add_option("-N,--nevent", nevent, "The number of events/showers to run.")
      ->default_val(1)
      ->check(CLI::PositiveNumber)
      ->group("Library/Output");
  app.add_option("-f,--filename", "Filename for output library.")
      ->required()
      ->default_val("corsika_library")
      ->check(CLI::NonexistentPath)
      ->group("Library/Output");
  app.add_option("-d,--dir", "Directory for output library.")
      ->default_val(boost::filesystem::current_path().string())
      ->group("Library/Output");
  app.add_option("-s,--seed", "The random number seed.")
      ->default_val(0)
      ->check(CLI::NonNegativeNumber)
      ->group("Misc.");
  bool force_interaction = false;
  app.add_option("-v,--verbosity", "Verbosity level: warn, info, debug, trace.")
      ->default_val("info")
      ->check(CLI::IsMember({"warn", "info", "debug", "trace"}))
      ->group("Misc.");
  app.add_option("-M,--hadronModel", "High-energy hadronic interaction model")
      ->default_val("SIBYLL-2.3d")
      ->check(CLI::IsMember({"SIBYLL-2.3d", "QGSJet-II.04", "EPOS-LHC"}))
      ->group("Misc.");
  app.add_option("-T,--hadronModelTransitionEnergy",
                 "Transition between high-/low-energy hadronic interaction "
                 "model in GeV")
      ->default_val(std::pow(10, 1.9)) // 79.4 GeV
      ->check(CLI::NonNegativeNumber)
      ->group("Misc.");
  // parse the command line options into the variables
  CLI11_PARSE(app, argc, argv);

  if (app.count("--verbosity")) {
    auto const loglevel = app["--verbosity"]->as<std::string>();
    if (loglevel == "warn") {
      logging::set_level(logging::level::warn);
    } else if (loglevel == "info") {
      logging::set_level(logging::level::info);
    } else if (loglevel == "debug") {
      logging::set_level(logging::level::debug);
    } else if (loglevel == "trace") {
#ifndef _C8_DEBUG_
      CORSIKA_LOG_ERROR("trace log level requires a Debug build.");
      return 1;
#endif
      logging::set_level(logging::level::trace);
    }
  }

  // check that we got either PDG or A/Z
  // this can be done with option_groups but the ordering
  // gets all messed up
  if (app.count("--pdg") == 0) {
    if ((app.count("-A") == 0) || (app.count("-Z") == 0)) {
      CORSIKA_LOG_ERROR("If --pdg is not provided, then both -A and -Z are required.");
      return 1;
    }
  }

  // initialize random number sequence(s)
  auto seed = registerRandomStreams(app["--seed"]->as<long>());

  /* === START: SETUP ENVIRONMENT AND ROOT COORDINATE SYSTEM === */

  EnvType env;
  const CoordinateSystemPtr& rootCS = env.getCoordinateSystem();
  const Point center{rootCS, 0_m, 0_m, 0_m};
  create_my_5layer_atmosphere<IMediumType, MyExtraEnv>(
          env, AtmosphereId::LinsleyUSStd, center, Medium::AirDry1Atm, MagneticFieldVector{rootCS, 50_uT, 0_T, 0_T} );

  /* === START: CONSTRUCT PRIMARY PARTICLE === */

  // parse the primary ID as a PDG or A/Z code
  Code beamCode;

  // check if we want to use a PDG code instead
  if (app.count("--pdg") > 0) {
    beamCode = convert_from_PDG(PDGCode(app["--pdg"]->as<int>()));
  } else {
    // check manually for proton and neutrons
    if ((A == 1) && (Z == 1))
      beamCode = Code::Proton;
    else if ((A == 1) && (Z == 0))
      beamCode = Code::Neutron;
    else
      beamCode = get_nucleus_code(A, Z);
  }
  HEPEnergyType mass = get_mass(beamCode);

  /* === START: CONSTRUCT GEOMETRY === */
  auto const observationHeight =
      app["--observation-level"]->as<double>() * 1_m + constants::EarthRadius::Mean;
  auto const atmosphere_height = app["--injection-height"]->as<double>() * 1_m;
  auto const injectionHeight = atmosphere_height + constants::EarthRadius::Mean;
  /* === END: CONSTRUCT GEOMETRY === */

  std::stringstream args;
  for (int i = 0; i < argc; ++i) { args << argv[i] << " "; }
  // create the output manager that we then register outputs with
  auto const outputDir = boost::filesystem::path(app["--dir"]->as<std::string>());
  OutputManager output(app["--filename"]->as<std::string>(), seed, args.str(), outputDir);

  DynamicInteractionProcess<StackType> heModel;

  // have SIBYLL always for PROPOSAL photo-hadronic interactions
  auto sibyll = std::make_shared<corsika::sibyll::Interaction>(env);

  if (auto const modelStr = app["--hadronModel"]->as<std::string>();
      modelStr == "SIBYLL-2.3d") {
    heModel = DynamicInteractionProcess<StackType>{sibyll};
  } else if (modelStr == "QGSJet-II.04") {
    heModel = DynamicInteractionProcess<StackType>{
        std::make_shared<corsika::qgsjetII::Interaction>()};
  } else if (modelStr == "EPOS-LHC") {
    heModel = DynamicInteractionProcess<StackType>{
        std::make_shared<corsika::epos::Interaction>()};
  } else {
    CORSIKA_LOG_CRITICAL("invalid choice \"{}\"; also check argument parser", modelStr);
    return EXIT_FAILURE;
  }

  InteractionCounter heCounted{heModel};

  corsika::pythia8::Decay decayPythia;

  //neutrino interactions with pythia (options are: NC, CC)
  bool NC = false;
  bool CC = false;
  if (auto const nuIntStr = app["--neutrino-interaction-type"]->as<std::string>();
     nuIntStr == "neutral" || nuIntStr == "NC") {
     NC = true;
     CC = false;
  } else if (nuIntStr == "charged" || nuIntStr == "CC") {
     NC = false;
     CC = true;
  } else if (nuIntStr == "both") {
     NC = true;
     CC = true;
  }
  corsika::pythia8::NeutrinoInteraction neutrinoPrimaryPythia(NC, CC);

  // hadronic photon interactions in resonance region
  corsika::sophia::InteractionModel sophia;

  HEPEnergyType const emcut = 1_GeV * app["--emcut"]->as<double>();
  HEPEnergyType const hadcut = 1_GeV * app["--hadcut"]->as<double>();
  HEPEnergyType const mucut = 1_GeV * app["--mucut"]->as<double>();
  HEPEnergyType const nucut = 1_GeV * app["--nucut"]->as<double>();
  ParticleCut cut(emcut, emcut, hadcut, mucut, !track_neutrinos);

  // tell proposal that we are interested in all energy losses above the particle cut
  set_energy_production_threshold(Code::Electron, std::min({emcut, hadcut, mucut}));
  set_energy_production_threshold(Code::Positron, std::min({emcut, hadcut, mucut}));
  set_energy_production_threshold(Code::Photon, std::min({emcut, hadcut, mucut}));
  set_energy_production_threshold(Code::MuMinus, std::min({emcut, hadcut, mucut}));
  set_energy_production_threshold(Code::MuPlus, std::min({emcut, hadcut, mucut}));
  set_energy_production_threshold(Code::TauMinus, std::min({emcut, hadcut, mucut}));
  set_energy_production_threshold(Code::TauPlus, std::min({emcut, hadcut, mucut}));
  
  set_kinetic_energy_propagation_threshold(Code::NuE, nucut);
  set_kinetic_energy_propagation_threshold(Code::NuEBar, nucut);
  set_kinetic_energy_propagation_threshold(Code::NuMu, nucut);
  set_kinetic_energy_propagation_threshold(Code::NuMuBar, nucut);
  set_kinetic_energy_propagation_threshold(Code::NuTau, nucut);
  set_kinetic_energy_propagation_threshold(Code::NuTauBar, nucut);

  
  // energy threshold for high energy hadronic model. Affects LE/HE switch for
  // hadron interactions and the hadronic photon model in proposal
  HEPEnergyType const heHadronModelThreshold =
      1_GeV * app["--hadronModelTransitionEnergy"]->as<double>();

  corsika::proposal::Interaction emCascade(
      env, sophia, sibyll->getHadronInteractionModel(), heHadronModelThreshold);

  // use BetheBlochPDG for hadronic continuous losses, and proposal otherwise
  corsika::proposal::ContinuousProcess emContinuous(env);

  corsika::urqmd::UrQMD leIntModel{};
  InteractionCounter leIntCounted{leIntModel};

  // assemble all processes into an ordered process list
  struct EnergySwitch {
    HEPEnergyType cutE_;
    EnergySwitch(HEPEnergyType cutE)
        : cutE_(cutE) {}
    bool operator()(const Particle& p) const { return (p.getKineticEnergy() < cutE_); }
  };
  auto hadronSequence =
      make_select(EnergySwitch(heHadronModelThreshold), leIntCounted, heCounted);

  // observation plane
  // observation plane sea level
  Point const seaPlaneCenter = Point(rootCS, {0_m, 0_m, constants::EarthRadius::Mean});
  Plane const seaPlane(seaPlaneCenter, DirectionVector(rootCS, {0., 0., 1.}));
  ObservationPlane<TrackingType, ParticleWriterParquet> seaobservationLevel{
                          seaPlane, DirectionVector(rootCS, {1., 0., 0.}),
                          false,   // plane should not "absorb" particles
                          false}; // do not print z-coordinate
  // register ground particle output
  //output.add("particles_sea_level", seaobservationLevel);
  
  //PrimaryWriter<TrackingType, ParticleWriterParquet> seaprimaryWriter(seaobservationLevel);
  //output.add("primary_sea", seaprimaryWriter);
  
  //observation plane 3km under sea level
  Point const detPlaneCenter = Point(rootCS, {0_m, 0_m, observationHeight});
  Plane const detPlane(detPlaneCenter, DirectionVector(rootCS, {0., 0., 1.}));
  ObservationPlane<TrackingType, ParticleWriterParquet> detobservationLevel{
                           detPlane, DirectionVector(rootCS, {1., 0., 0.}),
                           true,   // plane should "absorb" particles
                           false}; // do not print z-coordinate
  // register ground particle output
  output.add("particles_det_level", detobservationLevel);

  //PrimaryWriter<TrackingType, ParticleWriterParquet> detprimaryWriter(detobservationLevel);
  //output.add("primary_det", detprimaryWriter);
  
  // assemble the final process sequence
  auto sequence = make_sequence(neutrinoPrimaryPythia, hadronSequence,
                                decayPythia, emCascade, emContinuous, /*seaobservationLevel,*/ detobservationLevel, cut);

  /* === END: SETUP PROCESS LIST === */

  // create the cascade object using the default stack and tracking
  // implementation
  TrackingType tracking;
  StackType stack;
  Cascade EAS(env, tracking, sequence, output, stack);

  double eMin = app["--eMin"]->as<double>() * 1e9;
  double eMax = app["--eMax"]->as<double>() * 1e9;
  initial_energy_generator e_generator(eMin, eMax, -2, seed);
  TRandom3* rnd = new TRandom3(seed==0? 0 : (seed+1) );

  // trigger the output manager to open the library for writing
  output.startOfLibrary();

  // Record primary particles
  ofstream primout("Primaries.json");
  primout << "[\n";
  // loop over each shower
  for (int i_shower = 1; i_shower < nevent + 1; i_shower++) {

    CORSIKA_LOG_INFO("Shower {} / {} ", i_shower, nevent);
    // setup particle stack, and add primary particle
    stack.clear();

    //particle energy
    HEPEnergyType const E0 = 1_eV * e_generator.get_E_minus1();
    //HEPEnergyType const E0 = 1_GeV * app["--eMin"]->as<double>();
    std::cout<<"initial energy: "<<E0<<" So far so good"<<std::endl;
    
    Point const showerCore{rootCS, 0_m, 0_m, observationHeight};
    // direction of the shower in (theta, phi) space   
    auto const cos_theta = rnd->Uniform(0.9,1);
    auto const sin_theta = sqrt(1 - pow(cos_theta, 2));
    auto const phi = rnd->Uniform(0, 2 * TMath::Pi());
    DirectionVector inject_direction(rootCS, {-sin_theta*cos(phi), -sin_theta*sin(phi), -cos_theta});

    double geant4_halfZ = 285, geant4_radius = 2000;
    double R_circle = sqrt(geant4_halfZ*geant4_halfZ + geant4_radius*geant4_radius);

    DirectionVector tagent{rootCS, {0, 0, 1}};
    if (abs(inject_direction.getZ(rootCS)) < 0.99)
        tagent = DirectionVector{rootCS, {-inject_direction.getY(rootCS), inject_direction.getX(rootCS), 0}};
    else
        tagent = DirectionVector{rootCS, {-inject_direction.getZ(rootCS), 0, inject_direction.getX(rootCS)}};

    tagent = tagent / tagent.getNorm();
    DirectionVector bitagent = inject_direction.cross(tagent);
    double phi_ = rnd->Uniform(0, 2*TMath::Pi());
    LengthType radius = R_circle * sqrt(rnd->Uniform(0, 1))*1_m;
    Point const destPoint = detPlaneCenter + cos(phi_) * radius * tagent + sin(phi_) * radius * bitagent;
    LengthType dest_x = destPoint.getX(rootCS), dest_y = destPoint.getY(rootCS), dest_z = destPoint.getZ(rootCS);

    auto getInjectorLength = [atmosphere_height, dest_z, cos_theta](){
        double l1 = (dest_z)/1_m;
        double l2 = (atmosphere_height+constants::EarthRadius::Mean) / 1_m;
        double c = l1*l1 - l2*l2;
        double b = 2 * l1 * cos_theta;
        return 1./2 * (-b + sqrt(b*b - 4*c)) * 1_m;
    };
    auto injectorLength = getInjectorLength();
    Point const injectorPos = Point(rootCS, {destPoint.getX(rootCS)+injectorLength*sin_theta*cos(phi), destPoint.getY(rootCS)+injectorLength*sin_theta*sin(phi), dest_z+injectorLength*cos_theta});

    /* === END: CONSTRUCT PRIMARY PARTICLE === */

    // add the desired particle to the stack
    auto const primaryProperties = std::make_tuple(
        beamCode, E0 - get_mass(beamCode),
        inject_direction, injectorPos, 0_ns);
    stack.addParticle(primaryProperties);

    primout << "{\n";
    primout << "  \"shower\": " << i_shower-1 << ",\n";
    primout << "  \"pdg\": " << beamCode << ",\n  \"E\": " << E0/1_GeV << ",\n";
    primout << "  \"nx\": " << -sin_theta*cos(phi)<< ",\n";
    primout << "  \"ny\": " << -sin_theta*sin(phi) << ",\n";
    primout << "  \"nz\": " << -cos_theta << "\n";
    if(i_shower==nevent)
        primout << "}\n";
    else
        primout << "},\n";
    
    //seaprimaryWriter.recordPrimary(primaryProperties);
    //detprimaryWriter.recordPrimary(primaryProperties);
    // run the shower
    EAS.run();
    std::cout<<"Run with flying colours"<<std::endl;
  }
  
  primout << "]";
  primout.close();
  // and finalize the output on disk
  output.endOfLibrary();

  return EXIT_SUCCESS;
}
