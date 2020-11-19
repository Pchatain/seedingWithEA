#include "SeedFinderOptions.hpp"

#include "Acts/Seeding/Seedfinder.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/TelescopeSeeding/TelescopeSeedingAlgorithm.hpp"

#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

using namespace boost::program_options;

void ActsExamples::Options::addSeedFinderOptions(
    boost::program_options::options_description& opt) {
  opt.add_options()("sf-minPt", value<float>(), "Seed finder minimum pT.")(
      "sf-cotThetaMax", value<float>(), "cot of maximum theta angle")(
      "sf-deltaRMin", value<float>(),
      "Minimum distance in mm between two SPs in a seed")(
      "sf-deltaRMax", value<float>(),
      "Maximum distance in mm between two SPs in a seed")(
      "sf-impactMax", value<float>(), "max impact parameter in mm")(
      "sf-sigmaScattering", value<float>(),
      "How many sigmas of scattering to include in seeds")(
      "sf-maxSeedsPerSpM", value<size_t>(),
      "How many seeds can share one middle SpacePoint")(
      "sf-collisionRegionMin", value<float>(),
      "limiting location of collision region in z in mm")(
      "sf-collisionRegionMax", value<float>(),
      "limiting location of collision region in z in mm")(
      "sf-zMin", value<float>(),
      "Minimum z of space points included in algorithm")(
      "sf-zMax", value<float>(),
      "Maximum z of space points included in algorithm")(
      "sf-rMax", value<float>(),
      "Max radius of Space Points included in algorithm in mm")(
      "sf-rMin", value<float>(),
      "Min radius of Space Points included in algorithm in mm")(
      "sf-bFieldInZ", value<float>(), "Magnetic field strength in kiloTesla")(
      "sf-maxPt", value<float>(), "maximum Pt for scattering cut")(
      "sf-radLengthPerSeed", value<float>(), "Average radiation length");
}

void ActsExamples::Options::addSeedPerfOptions(
    boost::program_options::options_description& opt) {
  opt.add_options()("sf-fltPrtEtaMin", value<float>(),
                    "Min eta for particle in efficiency calculation")(
      "sf-fltPrtEtaMax", value<float>(),
      "Max eta for particle in efficiency calculation")(
      "sf-fltPrtPtMin", value<float>(),
      "Min transverse momentum for particles in efficiency calculation")(
      "sf-fltPrt3Hits", value<bool>(),
      "Whether to filter particles based on if they have 3 hits in pixel "
      "layers")("sf-fltPrtOuterHits", value<bool>(),
                "Whether to filter particles based on if they have hits in "
                "outer detector")(
      "sf-fltPrtPtMax", value<float>(),
      "Max transverse momentum for particles in efficiency calculation");
}

void ActsExamples::Options::addMLOutput(
    boost::program_options::options_description& opt) {
  opt.add_options()("output-ML", value<bool>(),
                    "True if we want output to be ML friendly");
}

bool ActsExamples::Options::readMLOutputConfig(
    const boost::program_options::variables_map& vm) {
  bool outputIsML = false;
  if (vm.count("output-ML")) {
    outputIsML = vm["output-ML"].as<bool>();
  }
  return outputIsML;
}

// Read the seed performance config
ActsExamples::TelescopeSeedingAlgorithm::Config
ActsExamples::Options::readSeedPerfConfig(
    const boost::program_options::variables_map& vm) {
  ActsExamples::TelescopeSeedingAlgorithm::Config cfg;
  if (vm.count("sf-fltPrtEtaMin")) {
    cfg.etaMin = vm["sf-fltPrtEtaMin"].as<float>();
  }
  if (vm.count("sf-fltPrtEtaMax")) {
    cfg.etaMax = vm["sf-fltPrtEtaMax"].as<float>();
  }
  if (vm.count("sf-fltPrtPtMin")) {
    cfg.ptMin = vm["sf-fltPrtPtMin"].as<float>();
  }
  if (vm.count("sf-fltPrtPtMax")) {
    cfg.ptMax = vm["sf-fltPrtPtMax"].as<float>();
  }
  if (vm.count("sf-fltPrt3Hits")) {
    cfg.fltPrt3Hits = vm["sf-fltPrt3Hits"].as<bool>();
  }
  if (vm.count("sf-fltPrtOuterHits")) {
    cfg.fltPrtOuterHits = vm["sf-fltPrtOuterHits"].as<bool>();
  }
  return cfg;
}

// Read the seed finder config.
Acts::SeedfinderConfig<LdmxSpacePoint>
ActsExamples::Options::readSeedFinderConfig(
    const boost::program_options::variables_map& vm) {
  Acts::SeedfinderConfig<LdmxSpacePoint> cfg;
  if (vm.count("sf-minPt")) {
    cfg.minPt = vm["sf-minPt"].as<float>();
  }
  if (vm.count("sf-cotThetaMax")) {
    cfg.cotThetaMax = vm["sf-cotThetaMax"].as<float>();
  }
  if (vm.count("sf-deltaRMin")) {
    cfg.deltaRMin = vm["sf-deltaRMin"].as<float>();
  }
  if (vm.count("sf-deltaRMax")) {
    cfg.deltaRMax = vm["sf-deltaRMax"].as<float>();
  }
  if (vm.count("sf-impactMax")) {
    cfg.impactMax = vm["sf-impactMax"].as<float>();
  }
  if (vm.count("sf-sigmaScattering")) {
    cfg.sigmaScattering = vm["sf-sigmaScattering"].as<float>();
  }
  if (vm.count("sf-maxSeedsPerSpM")) {
    cfg.maxSeedsPerSpM = vm["sf-maxSeedsPerSpM"].as<size_t>();
  }
  if (vm.count("sf-collisionRegionMin")) {
    cfg.collisionRegionMin = vm["sf-collisionRegionMin"].as<float>();
  }
  if (vm.count("sf-collisionRegionMax")) {
    cfg.collisionRegionMax = vm["sf-collisionRegionMax"].as<float>();
  }
  if (vm.count("sf-zMin")) {
    cfg.zMin = vm["sf-zMin"].as<float>();
  }
  if (vm.count("sf-zMax")) {
    cfg.zMax = vm["sf-zMax"].as<float>();
  }
  if (vm.count("sf-rMax")) {
    cfg.rMax = vm["sf-rMax"].as<float>();
  }
  if (vm.count("sf-rMin")) {
    cfg.rMin = vm["sf-rMin"].as<float>();
  }
  if (vm.count("sf-bFieldInZ")) {
    cfg.bFieldInZ = vm["sf-bFieldInZ"].as<float>();
  }
  // if (vm.count("sf-beamPos")) {
  //   auto beamPos = vm["sf-beamPos"].template as<read_range>();
  //   if (beamPos.size() != 2) {
  //     throw std::invalid_argument(
  //         "Beam position dimension handed over is wrong. It must be of "
  //         "dimension 2. Beam position specifies the x, y coordinates of the "
  //         "beam.");
  //   }
  //   cfg.beamPos = {beamPos.at(0), beamPos.at(1)};
  // }
  if (vm.count("sf-radLengthPerSeed")) {
    cfg.radLengthPerSeed = vm["sf-radLengthPerSeed"].as<float>();
  }
  // if (vm.count("sf-maxPt")) {
  //   cfg.maxPt = vm["sf-maxPt"].as<float>();
  // }
  return cfg;
}
