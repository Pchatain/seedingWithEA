#include "Acts/LdmxEDM/LdmxSpacePoint.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "ActsExamples/TelescopeSeeding/TelescopeSeedingAlgorithm.hpp"
#include <Acts/Utilities/Logger.hpp>

#include <string>

#include <boost/program_options.hpp>

namespace ActsExamples {
namespace Options {

/// Add seed finder options such as min pT for seeds.
void addSeedFinderOptions(boost::program_options::options_description& opt);

/// Options to filter particles when analyzing performance of seedfinder
void addSeedPerfOptions(boost::program_options::options_description& opt);

/// Add an option to output ML friendly output
void addMLOutput(boost::program_options::options_description& opt);

/// Read the seed finder config.
Acts::SeedfinderConfig<LdmxSpacePoint> readSeedFinderConfig(
    const boost::program_options::variables_map& vm);

ActsExamples::TelescopeSeedingAlgorithm::Config readSeedPerfConfig(
    const boost::program_options::variables_map& vm);

bool readMLOutputConfig(const boost::program_options::variables_map& vm);
}  // namespace Options
}  // namespace ActsExamples
