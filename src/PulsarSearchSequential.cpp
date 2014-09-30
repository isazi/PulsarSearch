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
#include <string>
#include <exception>
#include <fstream>
#include <vector>
#include <iomanip>
#include <algorithm>
#include <boost/mpi.hpp>

#include <configuration.hpp>

#include <ArgumentList.hpp>
#include <utils.hpp>
#include <Observation.hpp>
#include <ReadData.hpp>
#include <Timer.hpp>

#include <Shifts.hpp>
#include <Dedispersion.hpp>
#include <Folding.hpp>
#include <SNR.hpp>


int main(int argc, char * argv[]) {
	bool dataLOFAR = false;
	bool dataSIGPROC = false;
  unsigned int MPIRows = 0;
  unsigned int MPICols = 0;
	unsigned int bytesToSkip = 0;
  unsigned int secondsToBuffer = 0;
  unsigned int remainingSamples = 0;
	std::string dataFile;
	std::string headerFile;
	std::string outputFile;
  isa::utils::ArgumentList args(argc, argv);
	// Observation object
  AstroData::Observation obs;

	// Initialize MPI
	boost::mpi::environment envMPI(argc, argv);
  boost::mpi::communicator world;

	try {
		// Cols are associated with periods
		MPICols = args.getSwitchArgument< unsigned int >("-mpi_cols");
		// Rows are associated with DMs
		MPIRows = args.getSwitchArgument< unsigned int >("-mpi_rows");

		obs.setPadding(args.getSwitchArgument< unsigned int >("-padding"));

		dataLOFAR = args.getSwitch("-lofar");
		dataSIGPROC = args.getSwitch("-sigproc");
		if ( dataLOFAR && dataSIGPROC ) {
			std::cerr << "-lofar and -sigproc are mutually exclusive." << std::endl;
			throw std::exception();
		} else if ( dataLOFAR ) {
			headerFile = args.getSwitchArgument< std::string >("-header");
			dataFile = args.getSwitchArgument< std::string >("-data");
		} else if ( dataSIGPROC ) {
			bytesToSkip = args.getSwitchArgument< unsigned int >("-header");
			dataFile = args.getSwitchArgument< std::string >("-data");
			obs.setNrSeconds(args.getSwitchArgument< unsigned int >("-seconds"));
      obs.setFrequencyRange(args.getSwitchArgument< unsigned int >("-channels"), args.getSwitchArgument< float >("-min_freq"), args.getSwitchArgument< float >("-channel_bandwidth"));
			obs.setNrSamplesPerSecond(args.getSwitchArgument< unsigned int >("-samples"));
		} else {
			std::cerr << "Need to specify the -header and -data arguments." << std::endl;
			throw std::exception();
		}
		outputFile = args.getSwitchArgument< std::string >("-output");

    unsigned int tempUInts[3] = {args.getSwitchArgument< unsigned int >("-dm_node"), 0, 0};
    float tempFloats[2] = {args.getSwitchArgument< float >("-dm_first"), args.getSwitchArgument< float >("-dm_step")};
    obs.setDMRange(tempUInts[0], tempFloats[0] + ((world.rank() / MPIRows) * tempUInts[0] * tempFloats[1]), tempFloats[1]);
    tempUInts[0] = args.getSwitchArgument< unsigned int >("-period_node");
    tempUInts[1] = args.getSwitchArgument< unsigned int >("-period_first");
    tempUInts[2] = args.getSwitchArgument< unsigned int >("-period_step");
    obs.setPeriodRange(tempUInts[0], tempUInts[1] + ((world.rank() % MPICols) * tempUInts[0] * tempUInts[2]), tempUInts[2]);
		obs.setNrBins(args.getSwitchArgument< unsigned int >("-period_bins"));
	} catch ( isa::utils::EmptyCommandLine & err ) {
    std::cerr <<  args.getName() << " -mpi_cols ... -mpi_rows ... -padding ... [-lofar] [-sigproc] -output ... -dm_node ... -dm_first ... -dm_step ... -period_node ... -period_first ... -period_step ... -period_bins ..."<< std::endl;
    std::cerr << "\t -lofar -header ... -data ..." << std::endl;
    std::cerr << "\t -sigproc -header ... -data ... -seconds ... -channels ... -min_freq ... -channel_bandwidth ... -samples ..." << std::endl;
    return 1;
  } catch ( std::exception & err ) {
		std::cerr << err.what() << std::endl;
		return 1;
	}

	// Load observation data
  isa::utils::Timer loadTime;
	std::vector< std::vector< dataType > * > * input = new std::vector< std::vector< dataType > * >(1);
	loadTime.start();
	if ( dataLOFAR ) {
    AstroData::readLOFAR(headerFile, dataFile, obs, *input);
	} else if ( dataSIGPROC ) {
		input->resize(obs.getNrSeconds());
    AstroData::readSIGPROC(obs, bytesToSkip, dataFile, *input);
	}
	loadTime.stop();
	if ( DEBUG && world.rank() == 0 ) {
		std::cout << "Time to load the input: " << std::fixed << std::setprecision(6) << loadTime.getTotalTime() << " seconds." << std::endl;
	}

	// Host memory allocation
  std::vector< unsigned int > * shifts = PulsarSearch::getShifts(obs);
  obs.setNrSamplesPerDispersedChannel(obs.getNrSamplesPerSecond() + shifts->at(((obs.getNrDMs() - 1) * obs.getNrPaddedChannels())));
  secondsToBuffer = obs.getNrSamplesPerDispersedChannel() / obs.getNrSamplesPerSecond();
  remainingSamples = obs.getNrSamplesPerDispersedChannel() % obs.getNrSamplesPerSecond();
  std::vector< dataType > dispersedData(obs.getNrChannels() * obs.getNrSamplesPerDispersedChannel());
  std::vector< dataType > dedispersedData(obs.getNrDMs() * obs.getNrSamplesPerPaddedSecond());
  std::vector< dataType > foldedData(obs.getNrDMs() * obs.getNrPeriods() * obs.getNrPaddedBins());
  std::vector < unsigned int > counterData(obs.getNrDMs() * obs.getNrPeriods() * obs.getNrPaddedBins());
  std::vector< float > snrTable(obs.getNrDMs() * obs.getNrPaddedPeriods());

	// Search loop
  isa::utils::Timer searchTime;
  isa::utils::Timer inputLoadTime;
  isa::utils::Timer dedispTime;
  isa::utils::Timer foldTime;
  isa::utils::Timer snrTime;

  world.barrier();
  searchTime.start();
	if ( DEBUG && world.rank() == 0 ) {
		std::cout << "Starting the search." << std::endl;
		std::cout << "Processing seconds: ";
	}
	for ( unsigned int second = 0; second < obs.getNrSeconds() - secondsToBuffer; second++ ) {
		if ( DEBUG && world.rank() == 0 ) {
			std::cout << second << " " << std::flush;
		}
		// Load the input
		inputLoadTime.start();
		for ( unsigned int channel = 0; channel < obs.getNrChannels(); channel++ ) {
			for ( unsigned int chunk = 0; chunk < secondsToBuffer; chunk++ ) {
        memcpy(reinterpret_cast< void * >(&(dispersedData.data()[(channel * obs.getNrSamplesPerDispersedChannel()) + (chunk * obs.getNrSamplesPerSecond())])), reinterpret_cast< void * >(&((input->at(second + chunk))->at(channel * obs.getNrSamplesPerPaddedSecond()))), obs.getNrSamplesPerSecond() * sizeof(dataType));
			}
      memcpy(reinterpret_cast< void * >(&(dispersedData.data()[(channel * obs.getNrSamplesPerDispersedChannel()) + (secondsToBuffer * obs.getNrSamplesPerSecond())])), reinterpret_cast< void * >(&((input->at(second + secondsToBuffer))->at(channel * obs.getNrSamplesPerPaddedSecond()))), remainingSamples * sizeof(dataType));
		}
		inputLoadTime.stop();

		// Run the kernels
    dedispTime.start();
    PulsarSearch::dedispersion(obs, dispersedData, dedispersedData, *shifts);
    dedispTime.stop();
    foldTime.start();
    PulsarSearch::folding(second, obs, dedispersedData, foldedData, counterData);
    foldTime.stop();
	}
	if ( DEBUG && world.rank() == 0 ) {
		std::cout << "." << std::endl;
		std::cout << "Search complete." << std::endl;
	}
  searchTime.stop();

	// Store output
	if ( DEBUG && world.rank() == 0 ) {
		std::cout << "Analyzing processed data." << std::endl;
	}
  snrTime.start();
  PulsarSearch::snrFoldedTS(obs, foldedData, snrTable);
  snrTime.stop();

	if ( DEBUG && world.rank() == 0 ) {
		std::cout << "Saving output to disk." << std::endl;
	}
	std::ofstream output;
	output.open(outputFile + "_" + isa::utils::toString(world.rank()) + ".dat");
  output << "# period DM SNR" << std::endl;
  output << std::fixed << std::setprecision(6);
	for ( unsigned int period = 0; period < obs.getNrPeriods(); period++ ) {
		for ( unsigned int dm = 0; dm < obs.getNrDMs(); dm++ ) {
      output << ((world.rank() % MPICols) * obs.getNrPeriods()) + period << " ";
      output << ((world.rank() / MPIRows) * obs.getNrDMs()) + dm << " ";
      output << snrTable[(period * obs.getNrPaddedDMs()) + dm] << std::endl;
		}
	}
	output.close();
	output.open(outputFile + "_" + isa::utils::toString(world.rank()) + ".stat");
  output << "# searchTime inputLoadTotal inputLoadAvg err dedispersionTotal dedispersionvg err foldingTotal foldingAvg err snrTotal snrAvg err" << std::endl;
  output << std::fixed << std::setprecision(6);
  output << searchTime.getTotalTime() << " ";
  output << inputLoadTime.getTotalTime() << " " << inputLoadTime.getAverageTime() << " " << inputLoadTime.getStdDev() << " ";
  output << dedispTime.getTotalTime() << " " << dedispTime.getAverageTime() << " " << dedispTime.getStdDev() << " ";
  output << foldTime.getTotalTime() << " " << foldTime.getAverageTime() << " " << foldTime.getStdDev() << " ";
  output << snrTime.getTotalTime() << " " << snrTime.getAverageTime() << " " << snrTime.getStdDev() << " ";
  output << std::endl;

  if ( DEBUG && world.rank() == 0 ) {
    std::cout << "Output and statistics saved to disk." << std::endl;
  }
  world.barrier();

	return 0;
}

