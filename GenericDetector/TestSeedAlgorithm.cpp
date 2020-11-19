// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "TestSeedAlgorithm.hpp"

#include "ACTFW/EventData/GeometryContainers.hpp"
#include "ACTFW/EventData/IndexContainers.hpp"
#include "ACTFW/EventData/ProtoTrack.hpp"
#include "ACTFW/EventData/SimHit.hpp"
#include "ACTFW/EventData/SimIdentifier.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/EventData/SimVertex.hpp"
#include "ACTFW/Framework/WhiteBoard.hpp"
#include "ACTFW/Io/Csv/CsvPlanarClusterReader.hpp"
#include "ACTFW/Utilities/Range.hpp"
#include "ACTFW/Validation/ProtoTrackClassification.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Seeding/ATLASCuts.hpp"
#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SpacePoint.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsFatras/EventData/Barcode.hpp"

#include <chrono>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <boost/type_erasure/any_cast.hpp>

using HitParticlesMap = FW::IndexMultimap<ActsFatras::Barcode>;
using ProtoTrack = FW::ProtoTrack;

FW::TestSeedAlgorithm::TestSeedAlgorithm(
    const FW::TestSeedAlgorithm::Config& cfg, Acts::Logging::Level level)
    : FW::BareAlgorithm("TestSeedAlgorithm", level), m_cfg(std::move(cfg)) {
  // Is inputClusters already checked by base constructor?
  // I think this is only true for the writer
  if (m_cfg.inputClusters.empty()) {
    throw std::invalid_argument(
        "Missing clusters input collection with the hits");
  }
  if (m_cfg.inputHitParticlesMap.empty()) {
    throw std::invalid_argument("Missing hit-particles map input collection");
  }
  if (m_cfg.outputSeeds.empty()) {
    throw std::invalid_argument("Missing output collection for seeds");
  }
  if (m_cfg.inputParticles.empty()) {
    throw std::invalid_argument("Missing input particles collection");
  }
}

SpacePoint* FW::TestSeedAlgorithm::readSP(
    std::size_t hit_id, const Acts::GeometryID geoId,
    const Acts::PlanarModuleCluster& cluster,
    const HitParticlesMap& hitParticlesMap, const AlgorithmContext& ctx) const {
  const auto& parameters = cluster.parameters();
  Acts::Vector2D localPos(parameters[0], parameters[1]);
  Acts::Vector3D globalFakeMom(1, 1, 1);
  Acts::Vector3D globalPos(0, 0, 0);
  // transform local into global position information
  cluster.referenceObject().localToGlobal(ctx.geoContext, localPos,
                                          globalFakeMom, globalPos);
  float x, y, z, r, varianceR, varianceZ;
  x = globalPos.x();
  y = globalPos.y();
  z = globalPos.z();
  r = std::sqrt(x * x + y * y);
  varianceR = 0;  // why these values?
  varianceZ = 0;

  // get truth particles that are a part of this space point
  std::vector<FW::ParticleHitCount> particleHitCount;
  for (auto hitParticle : makeRange(hitParticlesMap.equal_range(hit_id))) {
    auto particleId = hitParticle.second;
    // search for existing particle in the existing hit counts
    auto isSameParticle = [=](const ParticleHitCount& phc) {
      return (phc.particleId == particleId);
    };
    auto it = std::find_if(particleHitCount.begin(), particleHitCount.end(),
                           isSameParticle);
    // either increase count if we saw the particle before or add it
    if (it != particleHitCount.end()) {
      it->hitCount += 1;
    } else {
      particleHitCount.push_back({particleId, 1u});
    }
  }

  SpacePoint* SP = new SpacePoint{
      hit_id, x, y, z, r, geoId, varianceR, varianceZ, particleHitCount};
  return SP;
}

ProtoTrack FW::TestSeedAlgorithm::seedToProtoTrack(
    const Acts::Seed<SpacePoint>* seed) const {
  ProtoTrack track;
  track.reserve(seed->sp().size());
  for (std::size_t i = 0; i < seed->sp().size(); i++) {
    track.emplace_back(seed->sp()[i]->Id());
  }
  return track;
}

FW::ProcessCode FW::TestSeedAlgorithm::execute(
    const AlgorithmContext& ctx) const {
  // read in the hits
  const auto& clusters =
      ctx.eventStore.get<FW::GeometryIdMultimap<Acts::PlanarModuleCluster>>(
          m_cfg.inputClusters);
  // read in the map of hitId to particleId truth information
  const HitParticlesMap hitParticlesMap =
      ctx.eventStore.get<HitParticlesMap>(m_cfg.inputHitParticlesMap);
  const auto& particleHitsMap = invertIndexMultimap(hitParticlesMap);
  /*
  const auto& particles =
      ctx.eventStore.get<SimParticleContainer>(m_cfg.inputParticles);*/

  std::size_t nHitsTotal = hitParticlesMap.size();

  // create the space points
  std::size_t clustCounter = 0;
  std::size_t nIgnored = 0;
  std::vector<const SpacePoint*> spVec;
  // since clusters are ordered, we simply count the hit_id as we read
  // clusters. Hit_id isn't stored in a cluster. This is how
  // CsvPlanarClusterWriter did it.
  std::size_t hit_id = 0;
  for (const auto& entry : clusters) {
    Acts::GeometryID geoId = entry.first;
    const Acts::PlanarModuleCluster& cluster = entry.second;
    std::size_t volumeId = geoId.volume();
    std::size_t layerId = geoId.layer();

    // filter out hits that aren't part of the volumes and layers track seeding
    // is supposed to work on.
    if (volumeId == 8) {
      SpacePoint* SP = readSP(hit_id, geoId, cluster, hitParticlesMap, ctx);
      spVec.push_back(SP);
      clustCounter++;

    } else {
      nIgnored++;
    }
    hit_id++;
  }

  // set up the cuts applied to seed finder algorithm
  Acts::SeedfinderConfig<SpacePoint> config = m_cfg.seedFinderCfg; /*
  // silicon detector max
  ACTS_INFO("The minPt passed in is " << m_cfg.seedFinderCfg.minPt)
  config.rMax = 200;  // 200
  config.deltaRMin = 1.;
  config.deltaRMax = 80.;  // 160
  config.collisionRegionMin = -250.;
  config.collisionRegionMax = 250.;
  config.zMin = -2000.;       // -2800 affects numbber of bins
  config.zMax = 2000.;        // 2800
  config.maxSeedsPerSpM = 5;  // 5 treated as +1 by seedFilter.ipp
  // 2.7 eta
  config.cotThetaMax = 7.40627;   // 7.40627. so far no effect on performance
  config.sigmaScattering = 2.25;  // better efficiency, more duplicates
  config.radLengthPerSeed = 0.1;  // default 0.05

  config.minPt = m_cfg.seedFinderCfg.minPt;
  config.bFieldInZ = 0.00199724;  // 0.00199724

  config.beamPos = {0, 0};  // {-.5, -.5} worse performance on {-.5, -.5}
  config.impactMax =
      3;  // 10 at 15 very minor increase found in raw eff, none in eff
*/
  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>());
  auto topBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>());
  Acts::SeedFilterConfig sfconf;
  sfconf.maxSeedsPerSpM = config.maxSeedsPerSpM;
  // sfconf default values listed here
  // maxSeedsPerSpM = 10
  // compatSeedLimit = 2;
  // sfconf.deltaRMin = config.deltaRMin;
  // compatSeedWeight = 200.;
  // impactWeightFactor = 1.;
  // deltaInvHelixDiameter = 0.00003;
  Acts::ATLASCuts<SpacePoint> atlasCuts = Acts::ATLASCuts<SpacePoint>();
  config.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
      Acts::SeedFilter<SpacePoint>(sfconf, &atlasCuts));
  Acts::Seedfinder<SpacePoint> a(config);

  // covariance tool, sets covariances per spacepoint as required
  auto ct = [=](const SpacePoint& SP, float, float, float) -> Acts::Vector2D {
    return {SP.varianceR, SP.varianceZ};
  };

  // setup spacepoint grid config
  Acts::SpacePointGridConfig gridConf;
  gridConf.bFieldInZ = config.bFieldInZ;
  gridConf.minPt = config.minPt;
  gridConf.rMax = config.rMax;
  gridConf.zMax = config.zMax;
  gridConf.zMin = config.zMin;
  gridConf.deltaRMax = config.deltaRMax;
  gridConf.cotThetaMax = config.cotThetaMax;

  // create grid with bin sizes according to the configured geometry
  std::unique_ptr<Acts::SpacePointGrid<SpacePoint>> grid =
      Acts::SpacePointGridCreator::createGrid<SpacePoint>(gridConf);
  ACTS_DEBUG("We are about to make the grid")
  auto spGroup = Acts::BinnedSPGroup<SpacePoint>(spVec.begin(), spVec.end(), ct,
                                                 bottomBinFinder, topBinFinder,
                                                 std::move(grid), config);
  ACTS_DEBUG("We made the grid")
  std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVector;
  ACTS_DEBUG("Briana")
  auto start = std::chrono::system_clock::now();
  ACTS_DEBUG("Peter")
  auto groupIt = spGroup.begin();
  ACTS_DEBUG("James")
  auto endOfGroups = spGroup.end();
  ACTS_DEBUG("Charlie")
  if (groupIt == endOfGroups) {
    ACTS_DEBUG("The grid is empty")
  } else {
    ACTS_DEBUG("Now we execute seed finding algorithm")
    // actually executues the seed finding algoirthm here
    for (; !(groupIt == endOfGroups); ++groupIt) {
      seedVector.push_back(a.createSeedsForGroup(
          groupIt.bottom(), groupIt.middle(), groupIt.top()));
    }
  }
  ACTS_DEBUG("MArk")
  auto end = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_seconds = end - start;
  if (m_cfg.outputIsML) {
    std::cout << "mlTag"
              << ","
              << "time"
              << "," << elapsed_seconds.count() << std::endl;
  } else {
    ACTS_INFO("Time to create seeds: " << elapsed_seconds.count() << "s")
    ACTS_DEBUG("Number of regions: " << seedVector.size())
    ACTS_INFO("Number of hits used is: " << clustCounter << " --- "
                                         << 100 * clustCounter / nHitsTotal
                                         << "% usage")  // some of the hits
  }

  ProtoTrackContainer seeds;
  seeds.reserve(seedVector.size());
  // create prototracks for all input particles
  for (auto& regionVec : seedVector) {
    for (std::size_t iseed = 0; iseed < regionVec.size(); iseed++) {
      const Acts::Seed<SpacePoint>* seed = &regionVec[iseed];
      ProtoTrack track = seedToProtoTrack(seed);
      seeds.emplace_back(std::move(track));
    }
  }
  if (seeds.empty()) {
    ACTS_DEBUG("seeds are empty")
  }
  if (seedVector.empty()) {
    ACTS_DEBUG("Seed vector is empty")
  }
  // store proto tracks to be analyzed by TrackFindingPerforrmanceWriter
  ctx.eventStore.add(m_cfg.outputProtoSeeds, std::move(seeds));
  // store seeds to be analyzed by TrackSeedingPerformanceWriter
  ctx.eventStore.add(m_cfg.outputSeeds, std::move(seedVector));
  return FW::ProcessCode::SUCCESS;
}
