//
// Copyright (C) 2014
// Alessio Sclocco <a.sclocco@vu.nl>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.
//

#include <iostream>
using std::cerr;
using std::endl;
#include <fstream>
using std::ifstream;
#include <map>
using std::map;
#include <utility>
using std::make_pair;
#include <exception>
using std::exception;
#include <string>
using std::string;
using std::getline;
#include <limits>
using std::numeric_limits;
#include <cctype>
using std::isdigit;
#include <CImg.h>
using cimg_library::CImg;

#include <ArgumentList.hpp>
using isa::utils::ArgumentList;
#include <ColorMap.hpp>
using AstroData::Color;
using AstroData::getColorMap;
#include <utils.hpp>
using isa::utils::castToType;
#include <Exceptions.hpp>
using isa::Exceptions::EmptyCommandLine;


int main(int argc, char * argv[]) {
    unsigned int nrDMs = 0;
    unsigned int nrPeriods = 0;
    float minSNR = numeric_limits< float >::max();
    float maxSNR = numeric_limits< float >::min();
    float snrSpaceDim = 0.0f;
    float * snrSpace = 0;
    string outFilename;
    ifstream searchFile;

    try {
        ArgumentList args(argc, argv);

        searchFile.open(args.getSwitchArgument< string >("-input"));
        outFilename = args.getSwitchArgument< string >("-output");
        nrDMs = args.getSwitchArgument< unsigned int >("-dms");
        nrPeriods = args.getSwitchArgument< unsigned int >("-periods");
    } catch ( EmptyCommandLine err ) {
        cerr << argv[0] << " -input ... -output ..." << endl;
        return 1;
    } catch ( exception &err ) {
        cerr << err.what() << endl;
        return 1;
    }

    snrSpace = new float [nrDMs * nrPeriods];
    while ( ! searchFile.eof() ) {
        string temp;
        unsigned int splitPoint = 0;
        unsigned int DM = 0;
        unsigned int period = 0;
        float snr = 0.0f;

        getline(searchFile, temp);
        if ( ! isdigit(temp[0]) ) {
            continue;
        }
        splitPoint = temp.find(" ");
        period = castToType< string, unsigned int >(temp.substr(0, splitPoint));
        temp = temp.substr(splitPoint + 1);
        splitPoint = temp.find(" ");
        temp = temp.substr(splitPoint + 1);
        DM = castToType< string, unsigned int >(temp.substr(0, splitPoint));
        temp = temp.substr(splitPoint + 1);
        splitPoint = temp.find(" ");
        temp = temp.substr(splitPoint + 1);
        snr = castToType< string, float >(temp);

        if ( snr > maxSNR ) {
            maxSNR = snr;
        }
        if ( snr < minSNR ) {
            minSNR = snr;
        }

        snrSpace[(period * nrDMs) + DM] = snr;
    }
    searchFile.close();
    snrSpaceDim = maxSNR - minSNR;

    CImg< unsigned char > searchImage(nrDMs, nrPeriods, 1, 3);
    Color *colorMap = getColorMap();
    for ( unsigned int period = 0; period < nrPeriods; period++ ) {
        for ( unsigned int DM = 0; DM < nrDMs; DM++ ) {
            float snr = snrSpace[(period * nrDMs + DM];
            searchImage(DM, period, 0, 0) = (colorMap[static_cast< unsigned int >((snr * 257) / snrSpaceDim)]).getR();
            searchImage(DM, period, 0, 1) = (colorMap[static_cast< unsigned int >((snr * 257) / snrSpaceDim)]).getG();
            searchImage(DM, period, 0, 2) = (colorMap[static_cast< unsigned int >((snr * 257) / snrSpaceDim)]).getB();
        }
    }
    searchImage.save(outFilename.c_str());

    return 0;
}
