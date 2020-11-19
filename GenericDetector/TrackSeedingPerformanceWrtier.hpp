// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "ACTFW/EventData/GeometryContainers.hpp"
#include "ACTFW/EventData/ProtoTrack.hpp"
#include "ACTFW/EventData/SimParticle.hpp"
#include "ACTFW/Framework/WriterT.hpp"
#include "ACTFW/Validation/EffPlotTool.hpp"
#include "ACTFW/Validation/FakeRatePlotTool.hpp"
#include "ACTFW/Validation/ResPlotTool.hpp"
#include "ACTFW/Validation/TrackSummaryPlotTool.hpp"
#include "Acts/Plugins/Digitization/PlanarModuleCluster.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SpacePoint.hpp"

#include <mutex>

class TFile;
class TTree;

namespace FW {

/// Write out the effeciency of the seed finding algorithm.
///
/// Efficiency here is the number of particles found by a seed divided by the
/// total number of particles that pass quality cuts.
///
/// A common file can be provided for the writer to attach his TTree,
/// this is done by setting the Config::rootFile pointer to an existing file
///
/// Safe to use from multiple writer threads - uses a std::mutex lock.
class TrackSeedingPerformanceWriter final
    : public WriterT<std::vector<std::vector<Acts::Seed<SpacePoint>>>> {
 public:
  struct Config {
    /// Input hit to particles map
    std::string inputHitParticlesMap;
    /// Input truth particles collection.
    std::string inputParticles;
    /// Input seeds to be analyzed.
    std::string inputSeeds;
    /// Input seeds as proto tracks
    std::string inputProtoSeeds;
    /// input Clusters from the event#-hits.csv file for calculating efficiency.
    std::string inputClusters;
    /// Output directory.
    std::string outputDir;
    /// Output filename.
    std::string outputFilename = "performance_track_seeding.root";
    // Controls whether to print out ML friendly output
    bool outputIsML = false;
    // Tag for whether a line is ML output
    std::string mlTag = "mlTag";
    // The quality cuts to be applied when evaluating seed finder efficiency
    /// Whether to filter particles based on if they have 3 hits in pixel layers
    bool fltPrt3Hits = true;
    /// Wehther to filter particles based on if they have hits in outer detector
    bool fltPrtOuterHits = true;
    /// Maximum distance from the origin in the transverse plane
    double rhoMax = std::numeric_limits<double>::max();
    /// Maximum absolute distance from the origin along z
    double absZMax = std::numeric_limits<double>::max();
    // Truth particle kinematic cuts
    double phiMin = std::numeric_limits<double>::lowest();
    double phiMax = std::numeric_limits<double>::max();
    double etaMin = std::numeric_limits<double>::lowest();
    double etaMax = std::numeric_limits<double>::max();
    double absEtaMin = std::numeric_limits<double>::lowest();
    double absEtaMax = std::numeric_limits<double>::max();
    double ptMin = 0.0;
    double ptMax = std::numeric_limits<double>::max();
    /// Keep neutral particles
    bool keepNeutral = false;
    /// Requirement on number of recorded hits
    //@TODO: implement detector-specific requirements
    size_t nHitsMin = 0;
    size_t nHitsMax = std::numeric_limits<size_t>::max();

    /// Plot tool configurations.
    EffPlotTool::Config effPlotToolConfig;
    FakeRatePlotTool::Config fakeRatePlotConfig;
    TrackSummaryPlotTool::Config trackSummaryPlotToolConfig;
  };

  /// @brief Finds all the particles that are in common to all space points in
  /// the seed.
  /// @param seed The seed to be analyzed
  std::set<ActsFatras::Barcode> identifySharedParticles(
      const Acts::Seed<SpacePoint>* seed) const;

  /// @brief Returns true if we expect the seed finder to be able to find this
  /// particle
  /// @param prt The particle to see whether it's findable
  /// @param particleHitsMap inverted map from hitParticlesMap, maps fromo
  /// particle barcodes to hits. TODO: get type information
  /// @param clusters Used to get information on the hits that a particle has
  /// made. TODO: Verify this is a valid filter
  bool prtFindable(
      const ActsFatras::Particle& prt, const auto& particleHitsMap,
      const FW::GeometryIdMultimap<Acts::PlanarModuleCluster>& clusters) const;

  /// @brief Analyzes onse seed. Finds whether or not the seed contains a truth
  /// particle. Technically, space points can have multiple particles that are a
  /// part of them, so analyzeSeed finds how many particles are in common.
  /// @param seed The seed to be processed.
  /// @param hitParticlesMap map from hits to particles
  /// @param truthCount map from particles found to how many seeds found them
  /// @param fakeCount map from particles to how many fake seeds they were a
  /// part of
  std::set<ActsFatras::Barcode> analyzeSeed(
      const Acts::Seed<SpacePoint>* seed,
      const FW::IndexMultimap<ActsFatras::Barcode>& hitParticlesMap,
      std::unordered_map<ActsFatras::Barcode, std::size_t>& truthCount,
      std::unordered_map<ActsFatras::Barcode, std::size_t>& fakeCount) const;

  /// @brief Converts a seed to a proto track of 3 hits
  FW::ProtoTrack seedToProtoTrack(const Acts::Seed<SpacePoint>* seed) const;

  /// @brief Writes the fake rate and efficiency plots
  /// TODO: Add the rests of the fake rate plots.
  void writePlots(
      const std::vector<std::vector<Acts::Seed<SpacePoint>>>& seedVector,
      const IndexMultimap<ActsFatras::Barcode>& hitParticlesMap,
      const SimParticleContainer& particles,
      const FW::GeometryIdMultimap<Acts::PlanarModuleCluster>& clusters);

  //   void analyzeSeeds(
  //       const std::vector<std::vector<Acts::Seed<SpacePoint>>>& seedVector,
  //       const IndexMultimap<ActsFatras::Barcode>& hitParticlesMap,
  //       const SimParticleContainer& particles,
  //       const FW::GeometryIdMultimap<Acts::PlanarModuleCluster>& clusters);

  /// @brief Prints out the 3 space points inide of a seed and relevant
  /// information
  ///
  /// @param seed The seed to be printed
  void printSeed(const Acts::Seed<SpacePoint>* seed) const;

  /// Construct from configuration and log level.
  TrackSeedingPerformanceWriter(Config cfg, Acts::Logging::Level lvl);
  ~TrackSeedingPerformanceWriter() override;

  /// Finalize plots.
  ProcessCode endRun() final override;

 private:
  /// @brief Calls TTrees writing function, and the writePlots function
  /// These create two files, one that has premade plots, and another that can
  /// be browsed.
  ProcessCode writeT(const AlgorithmContext& ctx,
                     const std::vector<std::vector<Acts::Seed<SpacePoint>>>&
                         seedVector) final override;
  /// @brief Used for writing the TTrees
  struct Impl;
  std::unique_ptr<Impl> m_impl;

  Config m_cfg;
  /// Mutex used to protect multi-threaded writes.
  std::mutex m_writeMutex;
  TFile* m_outputFile{nullptr};
  /// Plot tool for fakerate
  FakeRatePlotTool m_fakeRatePlotTool;
  FakeRatePlotTool::FakeRatePlotCache m_fakeRatePlotCache;
  /// Plot tool for efficiency
  EffPlotTool m_effPlotTool;
  EffPlotTool::EffPlotCache m_effPlotCache;
  /// Plot tool for track hit info
  TrackSummaryPlotTool m_trackSummaryPlotTool;
  TrackSummaryPlotTool::TrackSummaryPlotCache m_trackSummaryPlotCache;
};

}  // namespace FW
