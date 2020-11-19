#include "ActsExamples/TelescopeSeeding/TelescopeSeedingAlgorithm.hpp"

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/SeedToTrackParamMaker.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <sstream>

ActsExamples::TelescopeSeedingAlgorithm::TelescopeSeedingAlgorithm(
    ActsExamples::TelescopeSeedingAlgorithm::Config cfg,
    Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("TelescopeSeedingAlgorithm", lvl),
      m_cfg(std::move(cfg)) {
  ACTS_INFO("Constructed TelescopeSeeding Algorithm");
  ACTS_INFO("Input SpacePoint collection:" << m_cfg.inputSpacePoints);
}

ActsExamples::ProcessCode ActsExamples::TelescopeSeedingAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // Configure the tools --- is it necessary to do the configuration in execute?

  Acts::SeedfinderConfig<LdmxSpacePoint> config = m_cfg.seedFinderCfg;

  // Tagger r max
  config.rMax = 1000.;
  // config.deltaRMin = 3.;
  // config.deltaRMax = 220.;
  config.collisionRegionMin = -50;
  config.collisionRegionMax = 50;
  config.zMin = -300;
  config.zMax = 300.;
  // config.maxSeedsPerSpM = 5;

  // More or less the max angle is something of the order of 50 / 600 (assuming
  // that the track hits all the layers) Theta for the seeder is like ATLAS eta,
  // so it's 90-lambda. Max lamba is of the order of ~0.1 so cotThetaMax will
  // be 1./tan(pi/2 - 0.1) ~ 1.4.
  config.cotThetaMax = 1.5;

  // cotThetaMax and deltaRMax matter to choose the binning in z. The bin size
  // is given by cotThetaMax*deltaRMax

  // config.sigmaScattering = 2.25;
  config.minPt = 500.;
  config.bFieldInZ = 1.5e-3;  // in kT
  // config.beamPos = {0, 0};    // units?
  // config.impactMax = 20.;

  // setup the spacepoint grid configuration

  // setup spacepoint grid config
  Acts::SpacePointGridConfig gridConf;
  gridConf.bFieldInZ = config.bFieldInZ;
  gridConf.minPt = config.minPt;
  gridConf.rMax = config.rMax;
  gridConf.zMax = config.zMax;
  gridConf.zMin = config.zMin;
  gridConf.deltaRMax = config.deltaRMax;
  gridConf.cotThetaMax = config.cotThetaMax;

  auto bottomBinFinder = std::make_shared<Acts::BinFinder<LdmxSpacePoint>>(
      Acts::BinFinder<LdmxSpacePoint>());

  auto topBinFinder = std::make_shared<Acts::BinFinder<LdmxSpacePoint>>(
      Acts::BinFinder<LdmxSpacePoint>());

  // The seed finder needs a seed filter instance

  // In the seed finder there is the correction for the beam axis (?), which you
  // could ignore if you set the penalty for high impact parameters. So removing
  // that in the seeder config.

  Acts::SeedFilterConfig seedFilter_cfg;
  seedFilter_cfg.impactWeightFactor = 0.;

  // For the moment no experiment dependent cuts are assigned to the filter
  config.seedFilter = std::make_unique<Acts::SeedFilter<LdmxSpacePoint>>(
      Acts::SeedFilter<LdmxSpacePoint>(seedFilter_cfg));

  Acts::Seedfinder<LdmxSpacePoint> seedFinder(config);

  // Retrieve the SpacePoints
  using ldmx_sps = std::vector<LdmxSpacePoint>;
  const auto& allSps = ctx.eventStore.get<ldmx_sps>(m_cfg.inputSpacePoints);

  // Select the ones to be used for seeding
  // I'll use layer 3,5,7
  int numSeeds = 0;
  int qualityCutFoundPrts = 0;
  float qualityCutPrts = 0;
  int nTrueSeeds = 0;
  int nDuplicateSeeds = 0;
  std::vector<const LdmxSpacePoint*> spVec;
  int toPrint = 5;
  auto within = [](double x, double min, double max) {
    return (min <= x) and (x < max);
  };
  bool hitLay3 = false;
  bool hitLay7 = false;
  bool hitLay9 = false;
  float pTFound = 0;
  float pZFound = 0;
  for (size_t isp = 0; isp < allSps.size(); isp++) {
    if (allSps[isp].layer() == 3 || allSps[isp].layer() == 7 ||
        allSps[isp].layer() == 9) {
      if (allSps[isp].trackID() == 1 &&
          within(allSps[isp].pZ(), m_cfg.ptMin, m_cfg.ptMax)) {
        qualityCutPrts = 1;
        pTFound = allSps[isp].pT();
        pZFound = allSps[isp].pZ();
        if (allSps[isp].layer() == 3) {
          hitLay3 = true;
        }
        if (allSps[isp].layer() == 7) {
          hitLay7 = true;
        }
        if (allSps[isp].layer() == 9) {
          hitLay9 = true;
        }
      }
      spVec.push_back(&(allSps[isp]));
    }
  }
  if (!(hitLay3 && hitLay7 && hitLay9) && m_cfg.fltPrt3Hits) {
    qualityCutPrts = 0;
  }

  // crosscheck
  int printCounter = 0;
  for (size_t isp = 0; isp < spVec.size(); isp++) {
    ACTS_INFO("spVec [" << isp << "] layer " << spVec[isp]->layer()
                        << ".  trackID " << spVec[isp]->trackID()
                        << " and id of " << spVec[isp]->id() << " and pT of "
                        << spVec[isp]->pT());
    ACTS_INFO("(x,y,z) = (" << spVec[isp]->x() << "," << spVec[isp]->y() << ","
                            << spVec[isp]->z() << ")");
    if (printCounter >= toPrint) {
      break;
    }
    printCounter++;
  }

  // create grid with bin sizes according to the configured geometry
  std::unique_ptr<Acts::SpacePointGrid<LdmxSpacePoint>> grid =
      Acts::SpacePointGridCreator::createGrid<LdmxSpacePoint>(gridConf);

  ACTS_DEBUG("Formed grid");

  // covariance tool, sets covariances per spacepoint as required  --- for the
  // moment 0 covariance
  auto ct = [=](const LdmxSpacePoint& sp, float, float,
                float) -> Acts::Vector2D {
    return {0., 0.};
  };

  ACTS_DEBUG("Defined covariance tool");

  // create the space point group
  auto spGroup = Acts::BinnedSPGroup<LdmxSpacePoint>(
      spVec.begin(), spVec.end(), ct, bottomBinFinder, topBinFinder,
      std::move(grid), config);

  ACTS_DEBUG("BinnedSPGroup defined");

  // seed vector
  std::vector<std::vector<Acts::Seed<LdmxSpacePoint>>> seedVector;

  // find the seeds
  auto groupIt = spGroup.begin();
  auto endOfGroups = spGroup.end();

  for (; !(groupIt == endOfGroups); ++groupIt) {
    seedVector.push_back(seedFinder.createSeedsForGroup(
        groupIt.bottom(), groupIt.middle(), groupIt.top()));
  }

  // The seed to track parameters fitter
  Acts::SeedToTrackParamMaker seedToTrackMaker;
  for (auto& outVec : seedVector) {
    numSeeds += outVec.size();
    for (size_t iSeed = 0; iSeed < outVec.size(); iSeed++) {
      std::array<double, 9> outData;
      Acts::Transform3D tP;

      tP.setIdentity();
      tP(0, 0) = 0;
      tP(0, 1) = 0;
      tP(0, 2) = 1;

      tP(1, 0) = 0;
      tP(1, 1) = 1;
      tP(1, 2) = 0;

      tP(2, 0) = 1;
      tP(2, 1) = 0;
      tP(2, 2) = 0;

      tP.translation().x() = -213.1;
      tP.translation().y() = 0.;
      tP.translation().z() = 0.;
      //   std::cout << "PF:: CHECK CHECK \n" << tP.translation() << std::endl;
      seedToTrackMaker.FitSeedAtlas(outVec[iSeed], outData, tP, 0.0015);
      //   std::cout << outData[0] << " " << outData[1] << " " << outData[2] <<
      //   " "
      //             << outData[3] << " " << outData[4] << std::endl;

      seedToTrackMaker.KarimakiFit(outVec[iSeed].sp(), outData);
      // h_p->Fill(outData[4]);
    }
  }

  ACTS_INFO(spVec.size() << " hits, " << seedVector.size() << " regions, "
                         << numSeeds << " seeds");

  int nSeeds = numSeeds;
  // analyze the seeds
  for (auto& outVec : seedVector) {
    for (size_t iSeed = 0; iSeed < outVec.size(); iSeed++) {
      bool trueSeed = outVec[iSeed].sp()[0]->trackID() == 1 &&
                      outVec[iSeed].sp()[1]->trackID() == 1 &&
                      outVec[iSeed].sp()[2]->trackID() == 1;
      if (trueSeed) {
        nTrueSeeds++;
        // pTFound = outVec[iSeed].sp()[2]->pT();
        // pZFound = outVec[iSeed].sp()[2]->pZ();
      }
    }
  }
  if (nTrueSeeds > 0) {
    nDuplicateSeeds = nTrueSeeds - 1;
    qualityCutFoundPrts = 1;
  }
  ACTS_INFO("nTrueSeeds " << nTrueSeeds)
  if (qualityCutPrts == 0) {
    ACTS_INFO("No particles that satiisfy these cuts to be found")
  }
  // else if (nSeeds == 0) {
  //   ACTS_INFO("No seeds found")
  // }
  else {
    std::stringstream outStream;
    outStream << "mlTag"
              << ","
              << "seeds" << nSeeds << ","
              << "eff" << 100 * qualityCutFoundPrts << ","
              << "true" << nTrueSeeds << ","
              << "dup" << nDuplicateSeeds << ","
              << "pT" << pTFound << ","
              << "pZ" << pZFound << std::endl;
    std::cout << outStream.str();  // This is now atomic, so we don't get
                                   // garbled output
  }

  return ActsExamples::ProcessCode::SUCCESS;
}
