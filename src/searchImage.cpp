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
#include <limits>
#include <cctype>
#include <CImg.h>

#include <ArgumentList.hpp>
#include <ColorMap.hpp>
#include <utils.hpp>


int main(int argc, char * argv[]) {
  bool gray = false;
  unsigned int nrDMs = 0;
  unsigned int nrPeriods = 0;
  float minSNR = std::numeric_limits< float >::max();
  float maxSNR = std::numeric_limits< float >::min();
  float snrSpaceDim = 0.0f;
  float * snrSpace = 0;
  std::string outFilename;
  std::ifstream searchFile;

  isa::utils::ArgumentList args(argc, argv);
  try {
    gray = args.getSwitch("-gray");
    outFilename = args.getSwitchArgument< std::string >("-output");
    nrDMs = args.getSwitchArgument< unsigned int >("-dms");
    nrPeriods = args.getSwitchArgument< unsigned int >("-periods");
  } catch ( isa::utils::EmptyCommandLine & err ) {
    std::cerr << args.getName() << " [-gray] -output ... -dms ... -periods ... input" << std::endl;
    return 1;
  } catch ( std::exception & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }

  snrSpace = new float [nrDMs * nrPeriods];

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

        if ( snr > maxSNR ) {
          maxSNR = snr;
        }
        if ( snr < minSNR ) {
          minSNR = snr;
        }

        snrSpace[(period * nrDMs) + DM] = snr;
      }
      searchFile.close();
    }
  } catch ( isa::utils::EmptyCommandLine & err ) {
    snrSpaceDim = maxSNR - minSNR;
  }

  cimg_library::CImg< unsigned char > searchImage;
  AstroData::Color * colorMap;
  if ( gray ) {
    searchImage = cimg_library::CImg< unsigned char >(nrDMs, nrPeriods, 1, 1);
  } else {
    searchImage = cimg_library::CImg< unsigned char >(nrDMs, nrPeriods, 1, 3);
    colorMap = AstroData::getColorMap();
  }
  for ( unsigned int period = 0; period < nrPeriods; period++ ) {
    for ( unsigned int DM = 0; DM < nrDMs; DM++ ) {
      float snr = snrSpace[(period * nrDMs) + DM];
      if ( gray ) {
        searchImage(DM, (nrPeriods - 1) - period, 0, 0) = 255 - static_cast< unsigned char >(((snr - minSNR) * 255.0) / snrSpaceDim);
      } else {
        searchImage(DM, (nrPeriods - 1) - period, 0, 0) = (colorMap[static_cast< unsigned int >(((snr - minSNR) * 256.0) / snrSpaceDim)]).getR();
        searchImage(DM, (nrPeriods - 1) - period, 0, 1) = (colorMap[static_cast< unsigned int >(((snr - minSNR) * 256.0) / snrSpaceDim)]).getG();
        searchImage(DM, (nrPeriods - 1) - period, 0, 2) = (colorMap[static_cast< unsigned int >(((snr - minSNR) * 256.0) / snrSpaceDim)]).getB();
      }
    }
  }
  searchImage.save(outFilename.c_str());

  delete [] snrSpace;
  return 0;
}

