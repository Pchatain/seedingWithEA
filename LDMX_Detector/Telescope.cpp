#include "Acts/Geometry/CuboidVolumeBuilder.hpp"
#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingGeometryBuilder.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/HomogeneousVolumeMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialProperties.hpp"
#include "Acts/Surfaces/LineSurface.hpp"
#include "Acts/Surfaces/PlanarBounds.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Units.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"

#include <limits>

// Options and printing tracking geo
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"
#include "ActsExamples/Options/CommonOptions.hpp"
#include "ActsExamples/Plugins/BField/BFieldOptions.hpp"
#include "ActsExamples/Plugins/Obj/ObjWriterOptions.hpp"
#include "ActsExamples/Propagation/PropagationOptions.hpp"

#include "SeedFinderOptions.hpp"

// Writers
#include "ActsExamples/Io/Root/RootPropagationStepsWriter.hpp"
#include "ActsExamples/Plugins/Obj/ObjPropagationStepsWriter.hpp"
#include "ActsExamples/Plugins/Obj/ObjTrackingGeometryWriter.hpp"

// Measurements
#include "Acts/EventData/Measurement.hpp"
#include "Acts/EventData/MeasurementHelpers.hpp"

// Bfield
#include "ActsExamples/Plugins/BField/ScalableBField.hpp"
#include <Acts/MagneticField/ConstantBField.hpp>
#include <Acts/MagneticField/InterpolatedBFieldMap.hpp>
#include <Acts/MagneticField/SharedBField.hpp>

// Propagator
#include "Acts/Propagator/Navigator.hpp"
#include "ActsExamples/Propagation/PropagationAlgorithm.hpp"
#include "ActsExamples/Utilities/Paths.hpp"
#include <Acts/Propagator/AtlasStepper.hpp>
#include <Acts/Propagator/EigenStepper.hpp>
#include <Acts/Propagator/Propagator.hpp>
#include <Acts/Propagator/StraightLineStepper.hpp>

// Root Reader

#include "ActsExamples/Io/Root/RootLdmxSimHitReader.hpp"

// Printers
#include "ActsExamples/Printers/PrintLdmxSpacePoints.hpp"

// Ldmx EDM
#include "Acts/LdmxEDM/LdmxSpacePoint.hpp"

// Seeder
#include "ActsExamples/TelescopeSeeding/TelescopeSeedingAlgorithm.hpp"

using namespace Acts;
using namespace UnitLiterals;

class LineSurfaceStub : public LineSurface {
 public:
  LineSurfaceStub() = delete;
  //
  LineSurfaceStub(std::shared_ptr<const Transform3D> htrans, double radius,
                  double halfz)
      : GeometryObject(), LineSurface(htrans, radius, halfz) { /* nop */
  }
  //
  LineSurfaceStub(std::shared_ptr<const Transform3D> htrans,
                  std::shared_ptr<const LineBounds> lbounds = nullptr)
      : GeometryObject(), LineSurface(htrans, lbounds) { /*nop */
  }
  //
  LineSurfaceStub(std::shared_ptr<const LineBounds> lbounds,
                  const DetectorElementBase& detelement)
      : GeometryObject(), LineSurface(lbounds, detelement) { /* nop */
  }

  //
  LineSurfaceStub(const LineSurfaceStub& ls)
      : GeometryObject(), LineSurface(ls) { /* nop */
  }

  LineSurfaceStub& operator=(const LineSurfaceStub& ls) {
    LineSurface::operator=(ls);
    return *this;
  }

  //
  LineSurfaceStub(const GeometryContext& gctx, const LineSurfaceStub& ls,
                  const Transform3D& t)
      : GeometryObject(), LineSurface(gctx, ls, t) { /* nop */
  }

  /// Return method for the Surface type to avoid dynamic casts
  SurfaceType type() const final { return Surface::Straw; }

  /// Simply return true to show object exists and is callable
  bool constructedOk() const { return true; }

  using Surface::normal;

  /// Return a Polyhedron for the surfaces
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param lseg is ignored for a perigee @note ignored
  ///
  /// @return A list of vertices and a face/facett description of it
  Polyhedron polyhedronRepresentation(const GeometryContext& /*gctx*/,
                                      size_t /*lseg*/) const final {
    return Polyhedron({}, {}, {});
  }
};

/// @class DetectorElementStub
///
/// This is a lightweight type of detector element,
/// it simply implements the base class.
class DetectorElementStub : public DetectorElementBase {
 public:
  DetectorElementStub() : DetectorElementBase() {}

  DetectorElementStub(std::shared_ptr<const Transform3D> transform)
      : DetectorElementBase(), m_elementTransform(std::move(transform)) {}

  /// Constructor for single sided detector element
  /// - bound to a Plane Surface
  ///
  /// @param transform is the transform that element the layer in 3D frame
  /// @param pBounds is the planar bounds for the planar detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  DetectorElementStub(
      std::shared_ptr<const Transform3D> transform,
      std::shared_ptr<const PlanarBounds> pBounds, double thickness,
      std::shared_ptr<const ISurfaceMaterial> material = nullptr)
      : DetectorElementBase(),
        m_elementTransform(std::move(transform)),
        m_elementThickness(thickness) {
    auto mutableSurface = Surface::makeShared<PlaneSurface>(pBounds, *this);
    mutableSurface->assignSurfaceMaterial(material);
    m_elementSurface = mutableSurface;
  }

  /// Constructor for single sided detector element
  /// - bound to a Line Surface
  ///
  /// @param transform is the transform that element the layer in 3D frame
  /// @param dBounds is the line bounds for the line like detector element
  /// @param thickness is the module thickness
  /// @param material is the (optional) Surface material associated to it
  DetectorElementStub(
      std::shared_ptr<const Transform3D> transform,
      std::shared_ptr<const LineBounds> lBounds, double thickness,
      std::shared_ptr<const ISurfaceMaterial> material = nullptr)
      : DetectorElementBase(),
        m_elementTransform(std::move(transform)),
        m_elementThickness(thickness) {
    auto mutableSurface = Surface::makeShared<LineSurfaceStub>(lBounds, *this);
    mutableSurface->assignSurfaceMaterial(material);
    m_elementSurface = mutableSurface;
  }

  ///  Destructor
  ~DetectorElementStub() override { /*nop */
  }

  /// Return local to global transform associated with this identifier
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///
  /// @note this is called from the surface().transform() in the PROXY mode
  const Transform3D& transform(const GeometryContext& gctx) const override;

  /// Return surface associated with this detector element
  const Surface& surface() const override;

  /// The maximal thickness of the detector element wrt normal axis
  double thickness() const override;

 private:
  /// the transform for positioning in 3D space
  std::shared_ptr<const Transform3D> m_elementTransform;
  /// the surface represented by it
  std::shared_ptr<const Surface> m_elementSurface{nullptr};
  /// the element thickness
  double m_elementThickness{0.};
};

inline const Transform3D& DetectorElementStub::transform(
    const GeometryContext& /*gctx*/) const {
  return *m_elementTransform;
}

inline const Surface& DetectorElementStub::surface() const {
  return *m_elementSurface;
}

inline double DetectorElementStub::thickness() const {
  return m_elementThickness;
}

int main(int argc, char** argv) {
  // std::cout << "Make default options" << std::endl;
  auto desc = ActsExamples::Options::makeDefaultOptions();
  ActsExamples::Options::addObjWriterOptions(desc);
  ActsExamples::Options::addSequencerOptions(desc);
  ActsExamples::Options::addPropagationOptions(desc);
  ActsExamples::Options::addRandomNumbersOptions(desc);
  ActsExamples::Options::addBFieldOptions(desc);
  ActsExamples::Options::addOutputOptions(desc);
  ActsExamples::Options::addSeedFinderOptions(desc);
  ActsExamples::Options::addSeedPerfOptions(desc);

  // std::cout <<s "parsing the options" << std::endl;
  auto vm = ActsExamples::Options::parse(desc, argc, argv);

  // Add the Acts::Sequencer
  ActsExamples::Sequencer sequencer(
      ActsExamples::Options::readSequencerConfig(vm));

  // Now read the standard options
  auto logLevel = ActsExamples::Options::readLogLevel(vm);

  // Setup the Event Store (WhiteBoard)

  size_t ievt = 0;
  size_t ialg = 0;

  // std::cout << "Evt store" << std::endl;
  ActsExamples::WhiteBoard eventStore(
      Acts::getDefaultLogger("EventStore#" + std::to_string(ievt), logLevel));

  // The algorithm context
  ActsExamples::AlgorithmContext context(ialg, ievt, eventStore);

  // The hit reader algorithm

  //---config
  ActsExamples::RootLdmxSimHitReader::Config readerCfg;
  // readerCfg.fileList =
  // {"/nfs/slac/g/ldmx/data/mc/v12/4gev_1e_ecal_pn_val/4gev_1e_ecal_pn_v12_ldmx-det-v12_run54838_seeds_109677_109678.root"};
  readerCfg.fileList = {
      "/nfs/slac/g/ldmx/data/mc/v12/4gev_1e_inclusive/v2.3.0-alpha/"
      "mc_v12-4GeV-1e-inclusive_run1310001_t1601628859_reco.root"};
  readerCfg.treeName = "LDMX_Events";
  readerCfg.outputCollection = "LdmxSimHitCollection";

  //---Create the reader and add it to the sequencer
  sequencer.addReader(std::make_shared<ActsExamples::RootLdmxSimHitReader>(
      readerCfg, logLevel));

  // std::cout << "Created the reader.." << std::endl;

  //---Print the hits out
  ActsExamples::PrintLdmxSpacePoints::Config printCfg;
  printCfg.inputCollection = readerCfg.outputCollection;

  bool printLdmxSPs = false;

  if (printLdmxSPs)
    sequencer.addAlgorithm(std::make_shared<ActsExamples::PrintLdmxSpacePoints>(
        printCfg, logLevel));

  //---Form seeds
  ActsExamples::TelescopeSeedingAlgorithm::Config seedAlgCfg =
      ActsExamples::Options::readSeedPerfConfig(vm);
  seedAlgCfg.seedFinderCfg = ActsExamples::Options::readSeedFinderConfig(vm);
  seedAlgCfg.inputSpacePoints = readerCfg.outputCollection;
  seedAlgCfg.outputProtoTracks = "ldmxProtoTracks";
  seedAlgCfg.outputSeeds = "ldmxSeeds";

  sequencer.addAlgorithm(
      std::make_shared<ActsExamples::TelescopeSeedingAlgorithm>(seedAlgCfg,
                                                                logLevel));

  // Create a cubic volume builder
  CuboidVolumeBuilder cvb;
  GeometryContext tgContext = GeometryContext();
  Material beryllium = {35.28_cm, 42.10_cm, 9.012, 4, 1.848_g / 1_cm3};
  Material vacuum = {0._cm, 0._cm, 0., 0., 0};
  Material silicon = {9.370_cm, 46.52_cm, 28.0855, 14, 2.329_g / 1_cm3};
  // X position
  std::vector<double> SX = {10 - 3,  10 + 3,  100 - 3, 100 + 3, 200 - 3,
                            200 + 3, 300 - 3, 300 + 3, 400 - 3, 400 + 3,
                            500 - 3, 500 + 3, 600 - 3, 600 + 3};

  // Y position
  std::vector<double> SY = {0., 0., 0., 0., 0., 0., 0.,
                            0., 0., 0., 0., 0., 0., 0.};

  // Z position
  std::vector<double> SZ = {0., 0., 0., 0., 0., 0., 0.,
                            0., 0., 0., 0., 0., 0., 0.};

  // Stereo angle
  std::vector<double> SA = {0.,   0.1, 0.,  -0.1, 0.,   0.1, 0.,
                            -0.1, 0.,  0.1, 0.,   -0.1, 0.,  0.1};

  // Surfaces configs
  std::vector<CuboidVolumeBuilder::SurfaceConfig> surfaceConfig;
  for (unsigned int i = 0; i < SX.size(); i++) {
    CuboidVolumeBuilder::SurfaceConfig cfg;
    cfg.position = {SX[i] * UnitConstants::mm, SY[i] * UnitConstants::mm,
                    SZ[i] * UnitConstants::mm};

    double rotationAngle = M_PI * 0.5;

    // 0 0 -1
    // 0 1 0
    // 1 0 0

    // This rotation is needed to have the plane orthogonal to the X direction.
    // Rotation of the surfaces
    Vector3D xPos(cos(rotationAngle), 0., sin(rotationAngle));
    Vector3D yPos(0., 1., 0.);
    Vector3D zPos(-sin(rotationAngle), 0., cos(rotationAngle));
    cfg.rotation.col(0) = xPos;
    cfg.rotation.col(1) = yPos;
    cfg.rotation.col(2) = zPos;

    double stereoAngle = SA[i];

    RotationMatrix3D stereoRotation = RotationMatrix3D::Identity();

    Vector3D xStereo(1., 0., 0.);
    Vector3D yStereo(0., cos(stereoAngle), sin(stereoAngle));
    Vector3D zStereo(0., -sin(stereoAngle), cos(stereoAngle));
    stereoRotation.row(0) = xStereo;
    stereoRotation.row(1) = yStereo;
    stereoRotation.row(2) = zStereo;

    cfg.rotation = stereoRotation * cfg.rotation;

    // std::cout<<cfg.rotation<<std::endl;

    // Boundaries of the surfaces. The bounds are in local coordinates. If the
    // cfg.rotation is identity, then they correspond to global.
    cfg.rBounds = std::make_shared<const RectangleBounds>(
        RectangleBounds(20.17_cm, 50_cm));
    // std::make_shared<const RectangleBounds>(RectangleBounds(0.1_mm, 0.2_mm));

    // Material of the surfaces - berillium
    MaterialProperties matProp(silicon, 0.320_mm);
    cfg.surMat = std::make_shared<HomogeneousSurfaceMaterial>(matProp);

    // Thickness of the detector element
    cfg.thickness = 320_um;

    cfg.detElementConstructor =
        [](std::shared_ptr<const Transform3D> trans,
           std::shared_ptr<const RectangleBounds> bounds, double thickness) {
          return new DetectorElementStub(trans, bounds, thickness);
        };
    surfaceConfig.push_back(cfg);

  }  // surfaces

  // Build the surfaces
  for (const auto& cfg : surfaceConfig) {
    std::shared_ptr<const PlaneSurface> pSur = cvb.buildSurface(tgContext, cfg);
  }

  // std::cout << "Formed " << surfaceConfig.size() << " Surface configs"
  //           << std::endl;
  // Layer Configurations
  std::vector<CuboidVolumeBuilder::LayerConfig> layerConfig;
  for (auto& sCfg : surfaceConfig) {
    CuboidVolumeBuilder::LayerConfig cfg;
    cfg.surfaceCfg = sCfg;
    layerConfig.push_back(cfg);
  }

  // std::cout << "Formed " << layerConfig.size() << " layer configs" <<
  // std::endl;

  // Test that 4 layers with surfaces can be built
  for (auto& cfg : layerConfig) {
    LayerPtr layer = cvb.buildLayer(tgContext, cfg);
  }

  for (auto& cfg : layerConfig) {
    cfg.surface = nullptr;
    cfg.active = true;
  }

  // Build volume configuration
  CuboidVolumeBuilder::VolumeConfig volumeConfig;
  volumeConfig.position = {0.5_m, 0., 0.};
  volumeConfig.length = {10_m, 10_m, 10_m};
  volumeConfig.layerCfg = layerConfig;
  volumeConfig.name = "TrackerVolume";
  volumeConfig.volumeMaterial =
      std::make_shared<HomogeneousVolumeMaterial>(vacuum);

  // Test the building
  std::shared_ptr<TrackingVolume> trVol =
      cvb.buildVolume(tgContext, volumeConfig);

  CuboidVolumeBuilder::Config config;
  config.position = {0., 0., 0.};
  config.length = {10_m, 1_m, 1_m};
  config.volumeCfg = {volumeConfig};

  cvb.setConfig(config);

  TrackingGeometryBuilder::Config tgbCfg;
  tgbCfg.trackingVolumeBuilders.push_back(
      [=](const auto& cxt, const auto& inner, const auto&) {
        return cvb.trackingVolume(cxt, inner, nullptr);
      });

  // std::cout << "Tracking Geometry Builder.." << std::endl;

  TrackingGeometryBuilder tgb(tgbCfg);

  std::unique_ptr<const TrackingGeometry> detector =
      tgb.trackingGeometry(tgContext);

  // Configure the writer

  // std::cout << "LogLevel" << std::endl;
  auto volumeLogLevel = Acts::Logging::Level(0);

  // std::cout << "Config the writer" << std::endl;
  auto tgObjWriterConfig =
      ActsExamples::Options::readObjTrackingGeometryWriterConfig(
          vm, "ObjTrackingGeometryWriter", volumeLogLevel);

  // std::cout << "Write" << std::endl;
  auto tgObjWriter = std::make_shared<ActsExamples::ObjTrackingGeometryWriter>(
      tgObjWriterConfig);

  if (vm["output-obj"].template as<bool>()) {
    // Write the tracking geometry object
    tgObjWriter->write(context, *detector);
  }

  std::shared_ptr<const TrackingGeometry> tGeometry = std::move(detector);

  // Visit the surfaces

  std::vector<const Surface*> surfaces;
  surfaces.reserve(SX.size());

  tGeometry->visitSurfaces([&](const Surface* surface) {
    if (surface and surface->associatedDetectorElement()) {
      // std::cout << "surface " << surface->geometryId() << " placed at: ("
      //           << surface->center(tgContext).transpose() << " )" <<
      //           std::endl;
      surfaces.push_back(surface);
    }
  });

  // std::cout << "There are " << surfaces.size() << " surfaces" << std::endl;

  // Create the random number engine
  auto randomNumberSvcCfg = ActsExamples::Options::readRandomNumbersConfig(vm);
  auto randomNumberSvc =
      std::make_shared<ActsExamples::RandomNumbers>(randomNumberSvcCfg);

  // Create BField service
  auto bFieldVar = ActsExamples::Options::readBField(vm);

  // Navigator
  Acts::Navigator rNavigator(tGeometry);

  rNavigator.resolvePassive = false;
  rNavigator.resolveMaterial = true;
  rNavigator.resolveSensitive = true;

  /*
  ConstantBField bField(Vector3D(0.,0.,0.));
  using RecoStepper = EigenStepper<ConstantBField>;
  RecoStepper rStepper(bField);
  using RecoPropagator = Propagator<RecoStepper, Navigator>;
  RecoPropagator rPropagator(rStepper,rNavigator);

  auto pAlgConfig =
      ActsExamples::Options::readPropagationConfig(vm, rPropagator);

  pAlgConfig.randomNumberSvc = randomNumberSvc;

  sequencer.addAlgorithm(
      std::make_shared<
      ActsExamples::PropagationAlgorithm<RecoPropagator>>(
          pAlgConfig, logLevel));
  */

  std::visit(
      [&](auto& bField) {
        // Resolve the bfield map and create the propgator
        using field_type =
            typename std::decay_t<decltype(bField)>::element_type;
        Acts::SharedBField<field_type> fieldMap(bField);

        using field_map_type = decltype(fieldMap);

        std::optional<std::variant<Acts::EigenStepper<field_map_type>,
                                   Acts::AtlasStepper<field_map_type>,
                                   Acts::StraightLineStepper>>
            var_stepper;

        // translate option to variant
        if (vm["prop-stepper"].template as<int>() == 0) {
          var_stepper = Acts::StraightLineStepper{};
        } else if (vm["prop-stepper"].template as<int>() == 1) {
          var_stepper = Acts::EigenStepper<field_map_type>{std::move(fieldMap)};
        } else if (vm["prop-stepper"].template as<int>() == 2) {
          var_stepper = Acts::AtlasStepper<field_map_type>{std::move(fieldMap)};
        }

        // resolve stepper, setup propagator
        std::visit(
            [&](auto& stepper) {
              using Stepper = std::decay_t<decltype(stepper)>;
              using Propagator = Acts::Propagator<Stepper, Acts::Navigator>;
              Propagator propagator(std::move(stepper), std::move(rNavigator));

              // Read the propagation config and create the algorithms
              auto pAlgConfig =
                  ActsExamples::Options::readPropagationConfig(vm, propagator);
              pAlgConfig.randomNumberSvc = randomNumberSvc;
              sequencer.addAlgorithm(
                  std::make_shared<
                      ActsExamples::PropagationAlgorithm<Propagator>>(
                      pAlgConfig, logLevel));
            },
            *var_stepper);
      },
      bFieldVar);

  // ---------------------------------------------------------------------------------
  // Output directory
  std::string outputDir = vm["output-dir"].template as<std::string>();
  auto psCollection = vm["prop-step-collection"].as<std::string>();

  if (vm["output-root"].template as<bool>()) {
    // Write the propagation steps as ROOT TTree
    ActsExamples::RootPropagationStepsWriter::Config pstepWriterRootConfig;
    pstepWriterRootConfig.collection = psCollection;
    pstepWriterRootConfig.filePath =
        ActsExamples::joinPaths(outputDir, psCollection + ".root");
    sequencer.addWriter(
        std::make_shared<ActsExamples::RootPropagationStepsWriter>(
            pstepWriterRootConfig));
  }

  return sequencer.run();
}
