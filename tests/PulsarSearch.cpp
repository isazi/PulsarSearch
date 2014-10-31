// Copyright 2013 Alessio Sclocco <a.sclocco@vu.nl>
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
#include <readConfiguration.hpp>

#include <ArgumentList.hpp>
#include <utils.hpp>
#include <Observation.hpp>
#include <Generator.hpp>
#include <InitializeOpenCL.hpp>
#include <Kernel.hpp>

#include <Shifts.hpp>
#include <Dedispersion.hpp>
#include <Transpose.hpp>
#include <Bins.hpp>
#include <Folding.hpp>
#include <SNR.hpp>


int main(int argc, char * argv[]) {
  bool random = false;
	unsigned int clPlatformID = 0;
	unsigned int clDeviceID = 0;
  unsigned int MPIRows = 0;
  unsigned int MPICols = 0;
  unsigned int secondsToBuffer = 0;
  unsigned int remainingSamples = 0;
	std::string deviceName;
	std::string outputFile;
	std::ofstream output;
  isa::utils::ArgumentList args(argc, argv);
	// Observation object
  AstroData::Observation obs;
  // Fake pulsar
  unsigned int period = 0;
  unsigned int width = 0;
  float DM = 0;

	// Initialize MPI
	boost::mpi::environment envMPI(argc, argv);
  boost::mpi::communicator world;

	try {
		// Cols are associated with periods
		MPICols = args.getSwitchArgument< unsigned int >("-mpi_cols");
		// Rows are associated with DMs
		MPIRows = args.getSwitchArgument< unsigned int >("-mpi_rows");
		clPlatformID = args.getSwitchArgument< unsigned int >("-opencl_platform");
		clDeviceID = args.getSwitchArgument< unsigned int >("-opencl_device");
		deviceName = args.getSwitchArgument< std::string >("-device_name");

		readPadding(padding, args.getSwitchArgument< std::string >("-padding_file"));
		readVectorWidth(vectorWidth, args.getSwitchArgument< std::string >("-vector_file"));
		readDedispersion(dedispersionParameters, args.getSwitchArgument< std::string >("-dedispersion_file"));
		readTranspose(transposeParameters, args.getSwitchArgument< std::string >("-transpose_file"));
		readFolding(foldingParameters, args.getSwitchArgument< std::string >("-folding_file"));
		readSNR(snrParameters, args.getSwitchArgument< std::string >("-snr_file"));

    random = args.getSwitch("-random");
		obs.setPadding(padding[deviceName]);

    period = args.getSwitchArgument< unsigned int >("-period");
    width = args.getSwitchArgument< unsigned int >("-width");
    DM = args.getSwitchArgument< float >("-dm");
    obs.setNrSeconds(args.getSwitchArgument< unsigned int >("-seconds"));
    obs.setFrequencyRange(args.getSwitchArgument< unsigned int >("-channels"), args.getSwitchArgument< float >("-min_freq"), args.getSwitchArgument< float >("-channel_bandwidth"));
    obs.setNrSamplesPerSecond(args.getSwitchArgument< unsigned int >("-samples"));
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
    std::cerr <<  args.getName() << " -mpi_cols ... -mpi_rows ... -opencl_platform ... -opencl_device ... -device_name ... -padding_file ... -vector_file ... -dedispersion_file ... -transpose_file ... -folding_file ... -snr_file ... [-random] -period ... -width ... -dm ... -seconds ... -channels ... -min_freq ... -channel_bandwidth ... -samples ... -output ... -dm_node ... -dm_first ... -dm_step ... -period_node ... -period_first ... -period_step ... -period_bins ..."<< std::endl;
    return 1;
  } catch ( std::exception & err ) {
		std::cerr << err.what() << std::endl;
		return 1;
	}

	// Load observation data
	std::vector< std::vector< dataType > * > * input = new std::vector< std::vector< dataType > * >(obs.getNrSeconds());
  AstroData::generatePulsar(period, width, DM, obs, *input, random);
	if ( world.rank() == 0 ) {
    std::cout << "Seconds: " << obs.getNrSeconds() << std::endl;
    std::cout << "Samples: " << obs.getNrSamplesPerSecond() << std::endl;
    std::cout << "Frequency range: " << obs.getMinFreq() << " " << obs.getMaxFreq() << std::endl;
    std::cout << "Channels: " << obs.getNrChannels() << " " << obs.getChannelBandwidth() << std::endl;
	}

	// Initialize OpenCL
	cl::Context * clContext = new cl::Context();
	std::vector< cl::Platform > * clPlatforms = new std::vector< cl::Platform >();
	std::vector< cl::Device > * clDevices = new std::vector< cl::Device >();
	std::vector< std::vector< cl::CommandQueue > > * clQueues = new std::vector< std::vector < cl::CommandQueue > >();

	try {
    isa::OpenCL::initializeOpenCL(clPlatformID, 1, clPlatforms, clContext, clDevices, clQueues);
	} catch ( isa::OpenCL::OpenCLError & err ) {
		std::cerr << err.what() << std::endl;
		return 1;
	}

	// Host memory allocation
  std::vector< unsigned int > * shifts = PulsarSearch::getShifts(obs);
  obs.setNrSamplesPerDispersedChannel(obs.getNrSamplesPerSecond() + (*shifts)[((obs.getNrDMs() - 1) * obs.getNrPaddedChannels())]);
  secondsToBuffer = obs.getNrSamplesPerDispersedChannel() / obs.getNrSamplesPerSecond();
  remainingSamples = obs.getNrSamplesPerDispersedChannel() % obs.getNrSamplesPerSecond();
  std::vector< unsigned int > * nrSamplesPerBin = PulsarSearch::getSamplesPerBin(obs);
  std::vector< dataType > dispersedData(obs.getNrChannels() * obs.getNrSamplesPerDispersedChannel());
  std::vector< dataType > foldedData(obs.getNrBins() * obs.getNrPeriods() * obs.getNrPaddedDMs());
  std::vector< float > snrTable(obs.getNrPeriods() * obs.getNrPaddedDMs());

  // Device memory allocation and data transfers
  cl::Buffer shifts_d, nrSamplesPerBin_d, dispersedData_d, dedispersedData_d, transposedData_d, foldedData_d, counterData0_d, counterData1_d, snrTable_d;

  try {
    shifts_d = cl::Buffer(*clContext, CL_MEM_READ_ONLY, shifts->size() * sizeof(unsigned int), 0, 0);
    nrSamplesPerBin_d = cl::Buffer(*clContext, CL_MEM_READ_ONLY, nrSamplesPerBin->size() * sizeof(unsigned int), 0, 0);
    dispersedData_d = cl::Buffer(*clContext, CL_MEM_READ_ONLY, dispersedData.size() * sizeof(dataType), 0, 0);
    dedispersedData_d = cl::Buffer(*clContext, CL_MEM_READ_WRITE, obs.getNrDMs() * obs.getNrSamplesPerPaddedSecond() * sizeof(dataType), 0, 0);
    transposedData_d = cl::Buffer(*clContext, CL_MEM_READ_WRITE, obs.getNrSamplesPerSecond() * obs.getNrPaddedDMs() * sizeof(dataType), 0, 0);
    foldedData_d = cl::Buffer(*clContext, CL_MEM_READ_WRITE, obs.getNrBins() * obs.getNrPeriods() * obs.getNrPaddedDMs() * sizeof(dataType), 0, 0);
    counterData0_d = cl::Buffer(*clContext, CL_MEM_READ_WRITE, obs.getNrBins() * obs.getNrPaddedPeriods() * sizeof(unsigned int), 0, 0);
    counterData1_d = cl::Buffer(*clContext, CL_MEM_READ_WRITE, obs.getNrBins() * obs.getNrPaddedPeriods() * sizeof(unsigned int), 0, 0);
    snrTable_d = cl::Buffer(*clContext, CL_MEM_WRITE_ONLY, obs.getNrPeriods() * obs.getNrPaddedDMs() * sizeof(float), 0, 0);

    // shifts_d
    clQueues->at(clDeviceID)[0].enqueueWriteBuffer(shifts_d, CL_TRUE, 0, shifts->size() * sizeof(unsigned int), reinterpret_cast< void * >(shifts->data()));
    // nrSamplesPerBin_d
    clQueues->at(clDeviceID)[0].enqueueWriteBuffer(nrSamplesPerBin_d, CL_TRUE, 0, nrSamplesPerBin->size() * sizeof(unsigned int), reinterpret_cast< void * >(nrSamplesPerBin->data()));
    delete nrSamplesPerBin;
    // foldedData_d
    std::vector< dataType > transferDataType(obs.getNrBins() * obs.getNrPeriods() * obs.getNrPaddedDMs());
    std::fill(transferDataType.begin(), transferDataType.end(), 0);
    clQueues->at(clDeviceID)[0].enqueueWriteBuffer(foldedData_d, CL_TRUE, 0, transferDataType.size() * sizeof(dataType), reinterpret_cast< void * >(transferDataType.data()));
    // counterData0_d and counterData1_d
    std::vector< unsigned int > transferUInt(obs.getNrBins() * obs.getNrPaddedPeriods());
    std::fill(transferUInt.begin(), transferUInt.end(), 0);
    clQueues->at(clDeviceID)[0].enqueueWriteBuffer(counterData0_d, CL_TRUE, 0, transferUInt.size() * sizeof(unsigned int), reinterpret_cast< void * >(transferUInt.data()));
    clQueues->at(clDeviceID)[0].enqueueWriteBuffer(counterData1_d, CL_TRUE, 0, transferUInt.size() * sizeof(unsigned int), reinterpret_cast< void * >(transferUInt.data()));
  } catch ( cl::Error & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }

	if ( world.rank() == 0 ) {
		double hostMemory = 0.0;
		double deviceMemory = 0.0;

    hostMemory += dispersedData.size() * sizeof(dataType);
    hostMemory += snrTable.size() * sizeof(float);
    deviceMemory += hostMemory;
    deviceMemory += shifts->size() * sizeof(unsigned int);
    deviceMemory += obs.getNrPeriods() * obs.getNrBins() * isa::utils::pad(2, obs.getPadding()) * sizeof(unsigned int);
    deviceMemory += obs.getNrDMs() * obs.getNrSamplesPerPaddedSecond() * sizeof(dataType);
    deviceMemory += obs.getNrSamplesPerSecond() * obs.getNrPaddedDMs() * sizeof(dataType);
    deviceMemory += obs.getNrBins() * obs.getNrPeriods() * obs.getNrPaddedDMs() * sizeof(dataType);
    deviceMemory += obs.getNrBins() * obs.getNrPaddedPeriods() * 2 * sizeof(unsigned int);

		std::cout << "Allocated host memory: " << std::fixed << std::setprecision(3) << isa::utils::giga(hostMemory) << " GB." << std::endl;
		std::cout << "Allocated device memory: " << std::fixed << std::setprecision(3) << isa::utils::giga(deviceMemory) << "GB." << std::endl;
	}

	// Generate OpenCL kernels
  std::string * code;
  cl::Kernel * dedispersionK, * foldingK, * transposeK, * snrK;

  code = PulsarSearch::getDedispersionOpenCL(dedispersionParameters[deviceName][obs.getNrDMs()][0], dedispersionParameters[deviceName][obs.getNrDMs()][1], dedispersionParameters[deviceName][obs.getNrDMs()][2], dedispersionParameters[deviceName][obs.getNrDMs()][3], dedispersionParameters[deviceName][obs.getNrDMs()][4], dataName, obs, *shifts);
	try {
    dedispersionK = isa::OpenCL::compile("dedispersion", *code, "-cl-mad-enable -Werror", *clContext, clDevices->at(clDeviceID));
	} catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
		return 1;
	}
  delete shifts;
  code = isa::OpenCL::getTransposeOpenCL(transposeParameters[deviceName][obs.getNrDMs()], obs.getNrSamplesPerSecond(), obs.getNrDMs(), obs.getPadding(), vectorWidth[deviceName], dataName);
  try {
    transposeK = isa::OpenCL::compile("transpose", *code, "-cl-mad-enable -Werror", *clContext, clDevices->at(clDeviceID));
  } catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }
  code = PulsarSearch::getFoldingOpenCL(foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][1], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][2], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][3], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][4], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][5], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][6], dataName, obs);
  try {
    foldingK = isa::OpenCL::compile("folding", *code, "-cl-mad-enable -Werror", *clContext, clDevices->at(clDeviceID));
  } catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }
  code = PulsarSearch::getSNROpenCL(snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0], snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][1], snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][2], snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][3], dataName, obs);
  try {
    snrK = isa::OpenCL::compile("snr", *code, "-cl-mad-enable -Werror", *clContext, clDevices->at(clDeviceID));
  } catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }

  // Set execution parameters
  cl::NDRange dedispersionGlobal(obs.getNrSamplesPerPaddedSecond() / dedispersionParameters[deviceName][obs.getNrDMs()][3], obs.getNrDMs() / dedispersionParameters[deviceName][obs.getNrDMs()][3]);
  cl::NDRange dedispersionLocal(dedispersionParameters[deviceName][obs.getNrDMs()][1], dedispersionParameters[deviceName][obs.getNrDMs()][2]);
  if ( world.rank() == 0 ) {
    std::cout << "Dedispersion global: " << obs.getNrSamplesPerPaddedSecond() / dedispersionParameters[deviceName][obs.getNrDMs()][3] << ", " << obs.getNrDMs() / dedispersionParameters[deviceName][obs.getNrDMs()][3] << std::endl;
    std::cout << "Dedispersion local: " << dedispersionParameters[deviceName][obs.getNrDMs()][1] << ", " << dedispersionParameters[deviceName][obs.getNrDMs()][2] << std::endl;
  }
  cl::NDRange transposeGlobal(obs.getNrSamplesPerPaddedSecond(), static_cast< unsigned int >(std::ceil(static_cast< double >(obs.getNrDMs()) / transposeParameters[deviceName][obs.getNrDMs()])));
  cl::NDRange transposeLocal(transposeParameters[deviceName][obs.getNrDMs()], 1);
  if ( world.rank() == 0 ) {
    std::cout << "Transpose global: " << obs.getNrSamplesPerPaddedSecond() << ", " << static_cast< unsigned int >(std::ceil(static_cast< double >(obs.getNrDMs()) / transposeParameters[deviceName][obs.getNrDMs()])) << std::endl;
    std::cout << "Transpose local: " << transposeParameters[deviceName][obs.getNrDMs()] << ", 1" << std::endl;
  }
  cl::NDRange foldingGlobal(obs.getNrPaddedDMs() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][6] / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][3], obs.getNrPeriods() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][6] / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][4], obs.getNrBins() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][6] / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][5]);
  cl::NDRange foldingLocal(foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][1], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][2]);
  if ( world.rank() == 0 ) {
    std::cout << "Folding global: " << obs.getNrPaddedDMs() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][6] / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][3] << ", " << obs.getNrPeriods() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][6] / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][4] << ", " << obs.getNrBins() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][6] / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][5] << std::endl;
    std::cout << "Folding local: " << foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0] << ", " << foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][1] << ", " << foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][2] << std::endl;
  }
  cl::NDRange snrGlobal(obs.getNrPaddedDMs() / snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][2], obs.getNrPeriods() / snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][3]);
  cl::NDRange snrLocal(snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0], snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][1]);
  if ( world.rank() == 0 ) {
    std::cout << "SNR global: " << obs.getNrPaddedDMs() / snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][2] << ", " << obs.getNrPeriods() / snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][3] << std::endl;
    std::cout << "SNR local: " << snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0] << ", " << snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][1] << std::endl;
  }

  dedispersionK->setArg(0, dispersedData_d);
  dedispersionK->setArg(1, dedispersedData_d);
  dedispersionK->setArg(2, shifts_d);
  transposeK->setArg(0, dedispersedData_d);
  transposeK->setArg(1, transposedData_d);
  foldingK->setArg(1, transposedData_d);
  foldingK->setArg(2, foldedData_d);
  foldingK->setArg(5, nrSamplesPerBin_d);
  snrK->setArg(0, foldedData_d);
  snrK->setArg(1, snrTable_d);

	// Search loop
  cl::Event syncPoint;

  world.barrier();
	if ( world.rank() == 0 ) {
		std::cout << "Processing seconds: ";
	}
	for ( unsigned int second = 0; second < obs.getNrSeconds() - secondsToBuffer; second++ ) {
		if ( world.rank() == 0 ) {
			std::cout << second << " " << std::flush;
		}
		// Load the input
		for ( unsigned int channel = 0; channel < obs.getNrChannels(); channel++ ) {
			for ( unsigned int chunk = 0; chunk < secondsToBuffer; chunk++ ) {
        memcpy(reinterpret_cast< void * >(&(dispersedData.data()[(channel * obs.getNrSamplesPerDispersedChannel()) + (chunk * obs.getNrSamplesPerSecond())])), reinterpret_cast< void * >(&((input->at(second + chunk))->at(channel * obs.getNrSamplesPerPaddedSecond()))), obs.getNrSamplesPerSecond() * sizeof(dataType));
			}
      memcpy(reinterpret_cast< void * >(&(dispersedData.data()[(channel * obs.getNrSamplesPerDispersedChannel()) + (secondsToBuffer * obs.getNrSamplesPerSecond())])), reinterpret_cast< void * >(&((input->at(second + secondsToBuffer))->at(channel * obs.getNrSamplesPerPaddedSecond()))), remainingSamples * sizeof(dataType));
		}
    try {
      clQueues->at(clDeviceID)[0].enqueueWriteBuffer(dispersedData_d, CL_TRUE, 0, dispersedData.size() * sizeof(dataType), reinterpret_cast< void * >(dispersedData.data()));
    } catch ( cl::Error & err ) {
      std::cerr << err.what() << std::endl;
      return 1;
    }

		// Run the kernels
		try {
      clQueues->at(clDeviceID)[0].enqueueNDRangeKernel(*dedispersionK, cl::NullRange, dedispersionGlobal, dedispersionLocal, 0, &syncPoint);
      syncPoint.wait();
      clQueues->at(clDeviceID)[0].enqueueNDRangeKernel(*transposeK, cl::NullRange, transposeGlobal, transposeLocal, 0, &syncPoint);
      syncPoint.wait();
      foldingK->setArg(0, second);
			if ( second % 2 == 0 ) {
        foldingK->setArg(3, counterData0_d);
        foldingK->setArg(4, counterData1_d);
			} else {
        foldingK->setArg(3, counterData1_d);
        foldingK->setArg(4, counterData0_d);
			}
      clQueues->at(clDeviceID)[0].enqueueNDRangeKernel(*foldingK, cl::NullRange, foldingGlobal, foldingLocal, 0, &syncPoint);
      syncPoint.wait();
      // Save the folded output every second
      output.open(outputFile + "_" + isa::utils::toString(world.rank()) + ".fold" + isa::utils::toString(second));
      output << "# bin SNR" << std::endl;
      output << std::fixed << std::setprecision(6);
      for ( unsigned int dm = 0; dm < obs.getNrDMs(); dm++ ) {
        for ( unsigned int period = 0; period < obs.getNrPeriods(); period++ ) {
          output << "# DM: " << dm << " period: " << period << std::endl;
          for ( unsigned int bin = 0; bin < obs.getNrBins(); bin++ ) {
            output << bin << " " << foldedData[(bin * obs.getNrPeriods() * obs.getNrPaddedDMs()) + (period * obs.getNrPaddedDMs()) + dm] << std::endl;
          }
          output << std::endl << std:: endl;
        }
      }
      output.close();
		} catch ( cl::Error & err ) {
			std::cerr << err.what() << std::endl;
			return 1;
		}
	}
	if ( world.rank() == 0 ) {
		std::cout << "." << std::endl;
		std::cout << "Search complete." << std::endl;
	}

	// Store output
	if ( world.rank() == 0 ) {
		std::cout << "Analyzing processed data." << std::endl;
	}
	try {
    clQueues->at(clDeviceID)[0].enqueueNDRangeKernel(*snrK, cl::NullRange, snrGlobal, snrLocal, 0, &syncPoint);
    syncPoint.wait();
    clQueues->at(clDeviceID)[0].enqueueReadBuffer(snrTable_d, CL_TRUE, 0, snrTable.size() * sizeof(float), reinterpret_cast< void * >(snrTable.data()));
    clQueues->at(clDeviceID)[0].enqueueReadBuffer(foldedData_d, CL_TRUE, 0, foldedData.size() * sizeof(dataType), reinterpret_cast< void * >(foldedData.data()));
  } catch ( cl::Error & err ) {
		std::cerr << err.what() << std::endl;
		return 1;
	}

	if ( world.rank() == 0 ) {
		std::cout << "Saving output to disk." << std::endl;
	}
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

  if ( world.rank() == 0 ) {
    std::cout << "Output saved to disk." << std::endl;
  }
  world.barrier();

	return 0;
}
