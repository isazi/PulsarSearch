// Copyright 2014 Alessio Sclocco <a.sclocco@vu.nl>
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#include <iostream>
#include <fstream>
#include <exception>
#include <string>
#include <utility>
#include <map>
#include <iomanip>

#include <ArgumentList.hpp>
#include <utils.hpp>


int main(int argc, char * argv[]) {
  // DMs
  float firstDM = 0.0f;
  float stepDM = 0.0f;
  // Periods
  unsigned int firstPeriod = 0;
  unsigned int stepPeriod = 0;
  // Sampling
  unsigned int nrSamplesPerSecond = 0;
  // I/O
  std::string outFilename;
  std::ifstream searchFile;
  std::ofstream outFile;
  // Data
  unsigned int percentile = 0;
  std::multimap< float, std::pair< unsigned int, unsigned int > > snrList;

  isa::utils::ArgumentList args(argc, argv);
  try {
    outFilename = args.getSwitchArgument< std::string >("-output");
    firstDM = args.getSwitchArgument< float >("-dm_first");
    stepDM = args.getSwitchArgument< float >("-dm_step");
    firstPeriod = args.getSwitchArgument< unsigned int >("-period_first");
    stepPeriod = args.getSwitchArgument< unsigned int >("-period_step");
    nrSamplesPerSecond = args.getSwitchArgument< unsigned int >("-samples");
    percentile = args.getSwitchArgument< unsigned int >("-percentile");
  } catch ( isa::utils::EmptyCommandLine & err ) {
    std::cerr << args.getName() << " -output ... -dm_first ... -dm_step ... -period_first ... -period_step ... -samples ... -percentile ... input" << std::endl;
    return 1;
  } catch ( std::exception & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }

  // Read the SNR data
  try {
    while ( true ) {
      searchFile.open(args.getFirst< std::string >());

      while ( ! searchFile.eof() ) {
        std::string temp;
        unsigned int splitPoint = 0;
        unsigned int DM = 0;
        unsigned int period = 0;
        float snr = 0.0f;

        std::getline(searchFile, temp);
        if ( ! std::isdigit(temp[0]) ) {
          continue;
        }
        splitPoint = temp.find(" ");
        period = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
        temp = temp.substr(splitPoint + 1);
        splitPoint = temp.find(" ");
        DM = isa::utils::castToType< std::string, unsigned int >(temp.substr(0, splitPoint));
        temp = temp.substr(splitPoint + 1);
        splitPoint = temp.find(" ");
        snr = isa::utils::castToType< std::string, float >(temp);

        snrList.insert(std::make_pair(snr, std::make_pair(DM, period)));
      }
      searchFile.close();
    }
  } catch ( isa::utils::EmptyCommandLine & err ) {
  }

  // Print the percentile
  unsigned int lowerLimit = (percentile * snrList.size()) / 100;
  unsigned int counter = snrList.size() - 1;
  outFile.open(outFilename);
  outFile << std::fixed << std::setprecision(6);
  outFile << "# DMIndex DM periodIndex period snr" << std::endl;
  for ( std::multimap< float, std::pair< unsigned int, unsigned int> >::const_iterator item = snrList.crend(); counter >= lowerLimit; counter-- ) {
  outFile << (*item).second.first << " " << firstDM + ((*item).second.first * stepDM) << " ";
  outFile << (*item).second.second << " " << (firstPeriod + ((*item).second.second * stepPeriod)) / static_cast< float >(nrSamplesPerSecond) << " ";
  outFile << (*item).first << std::endl;
  }
  outFile.close();

  return 0;
}

