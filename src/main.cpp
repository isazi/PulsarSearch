/*
 * Copyright (C) 2013
 * Alessio Sclocco <a.sclocco@vu.nl>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <configuration.hpp>

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;
#include <string>
using std::string;
#include <exception>
using std::exception;
#include <fstream>
using std::ofstream;
#include <vector>
using std::vector;

#include <ArgumentList.hpp>
using isa::utils::ArgumentList;
#include <Observation.hpp>
using AstroData::Observation;
#include <ReadData.hpp>
using AstroData::readSIGPROC;
using AstroData::readLOFAR;
#include <Dedispersion.hpp>
using PulsarSearch::Dedispersion;
#include <Folding.hpp>
using PulsarSearch::Folding;
#include <SNR.hpp>
using PulsarSearch::SNR;
#include <CLData.hpp>
using isa::OpenCL::CLData;


int main(int argc, char * argv[]) {
	bool dataLOFAR = false;
	bool dataSIGPROC = false;
	unsigned int clPlatformID = 0;
	unsigned int clDeviceID = 0;
	string dataFile;
	string headerFile;
	string outputFile;
	// Observation object
	Observation< dataType > obs("PulsarSearch", dataName);

	try {
		ArgumentList args(argc, argv);

		clPlatformID = args.getSwitchArgument< unsigned int >("-opencl_platform");
		clDeviceID = args.getSwitchArgument< unsigned int >("-opencl_device");

		dataLOFAR = args.getSwitch("-lofar");
		dataSIGPROC = args.getSwitch("-sigproc");
		if ( dataLOFAR && dataSIGPROC ) {
			cerr << "-lofar and -sigproc are mutually exclusive." << endl;
			throw exception();
		} else if ( dataLOFAR ) {
			headerFile = args.getSwitchArgument< string >("-header");
			dataFile = args.getSwitchArgument< string >("-data");
		} else if ( dataSIGPROC ) {
			dataFile = args.getSwitchArgument< string >("-data");
		} else {
			cerr << "Need to specify the -header and -data arguments." << endl;
			throw exception();
		}
		outputFile = args.getSwitchArgument< string >("-output");

		obs.setFirstDM(args.getSwitchArgument< float >("-dm_first"));
		obs.setDMStep(args.getSwitchArgument< float >("-dm_step"));
		obs.setNrDMs(args.getSwitchArgument< unsigned int >("-dm_number"));
		obs.setFirstPeriod(args.getSwitchArgument< unsigned int >("-period_first"));
		obs.setNrPeriods(args.getSwitchArgument< unsigned int >("-period_number"));
		obs.setNrBins(args.getSwitchArgument< unsigned int >("-period_bins"));
		obs.setPeriodStep(obs.getNrBins());
	} catch ( exception &err ) {
		cerr << err.what() << endl;
		return 1;
	}

	// Load observation data
	vector< CLData< dataType > * > * input = new vector< CLData< dataType > * >(1);
	if ( dataLOFAR ) {
		readLOFAR(headerFile, dataFile, observation, *input);
	} else if ( SIGPROC ) {
		// TODO: implement SIGPROC pipeline
	}

	// Store output

	return 0;
}
