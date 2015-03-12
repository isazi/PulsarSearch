// Copyright 2015 Alessio Sclocco <a.sclocco@vu.nl>
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
#include <iomanip>
#include <vector>

#include <ArgumentList.hpp>
#include <utils.hpp>
#include <Stats.hpp>


int main(int argc, char * argv[]) {
  // DMs
  unsigned int nrDMs = 0;
  // Periods
  unsigned int nrPeriods = 0;
  // I/O
  std::string outFilename;
  std::ifstream searchFile;
  std::ofstream outFile;
  // Data
  std::vector< isa::utils::Stats< float > > snrDM;
  std::vector< isa::utils::Stats< float > > snrP;

  isa::utils::ArgumentList args(argc, argv);
  try {
    nrDMs = args.getSwitchArgument< unsigned int >("-dms");
    nrPeriods = args.getSwitchArgument< unsigned int >("-periods");
    outFilename = args.getSwitchArgument< std::string >("-output");
  } catch ( isa::utils::EmptyCommandLine & err ) {
    std::cerr << args.getName() << " -dms ... -periods ... -output ... input" << std::endl;
    return 1;
  } catch ( std::exception & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }
  snrDM.resize(nrDMs);
  for ( unsigned int dm = 0; dm < nrDMs; dm++ ) {
    snrDM[dm] = isa::utils::Stats< float >();
  }
  snrP.resize(nrPeriods);
  for ( unsigned int period = 0; period < nrPeriods; period++ ) {
    snrP[period] = isa::utils::Stats< float >();
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

        snrDM[DM].addElement(snr);
        snrP[period].addElement(snr);
      }
      searchFile.close();
    }
  } catch ( isa::utils::EmptyCommandLine & err ) {
  }

  // Print mean, minimum and maximum SNR
  outFile.open(outFilename);
  outFile << std::fixed;
  outFile << "# DM mean min max" << std::endl;
  for ( unsigned int dm = 0; dm < nrDMs; dm++ ) {
    outFile << dm << " " << snrDM[dm].getMean() << " " << snrDM[dm].getMin() << " " << snrDM[dm].getMax() << std::endl;
  }
  outFile << std::endl << std::endl;
  outFile << "# period mean min max" << std::endl;
  for ( unsigned int period = 0; period < nrPeriods; period++ ) {
    outFile << period << " " << snrP[period].getMean() << " " << snrP[period].getMin() << " " << snrP[period].getMax() << std::endl;
  }
  outFile.close();

  return 0;
}

