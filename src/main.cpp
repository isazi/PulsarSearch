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
#include <ReadData.hpp>
#include <Timer.hpp>
#include <InitializeOpenCL.hpp>
#include <Kernel.hpp>

#include <Shifts.hpp>
#include <Dedispersion.hpp>
#include <Transpose.hpp>
#include <Bins.hpp>
#include <Folding.hpp>
#include <SNR.hpp>


int main(int argc, char * argv[]) {
	bool dataLOFAR = false;
	bool dataSIGPROC = false;
	unsigned int clPlatformID = 0;
	unsigned int clDeviceID = 0;
	unsigned int bytesToSkip = 0;
  unsigned int secondsToBuffer = 0;
  unsigned int remainingSamples = 0;
	std::string deviceName;
	std::string dataFile;
	std::string headerFile;
	std::string outputFile;
	// Observation object
  AstroData::Observation obs;

	// Initialize MPI
	boost::mpi::environment envMPI(argc, argv);
  boost::mpi::communicator world;

	try {
    isa::utils::ArgumentList args(argc, argv);

		// Cols are associated with periods
		unsigned int MPICols = args.getSwitchArgument< unsigned int >("-mpi_cols");
		// Rows are associated with DMs
		unsigned int MPIRows = args.getSwitchArgument< unsigned int >("-mpi_rows");
		clPlatformID = args.getSwitchArgument< unsigned int >("-opencl_platform");
		clDeviceID = args.getSwitchArgument< unsigned int >("-opencl_device");
		deviceName = args.getSwitchArgument< std::string >("-device_name");

		readPadding(padding, args.getSwitchArgument< std::string >("-padding_file"));
		readVectorWidth(vectorWidth, args.getSwitchArgument< std::string >("-vector_file"));
		readDedispersion(dedispersionParameters, args.getSwitchArgument< std::string >("-dedispersion_file"));
		readTranspose(transposeParameters, args.getSwitchArgument< std::string >("-transpose_file"));
		readFolding(foldingParameters, args.getSwitchArgument< std::string >("-folding_file"));
		readSNR(snrParameters, args.getSwitchArgument< std::string >("-snr_file"));

		obs.setPadding(padding[deviceName]);

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
			obs.setNrChannels(args.getSwitchArgument< unsigned int >("-channels"));
			obs.setNrSamplesPerSecond(args.getSwitchArgument< unsigned int >("-samples"));
			obs.setMinFreq(args.getSwitchArgument< float >("-low_freq"));
			obs.setChannelBandwidth(args.getSwitchArgument< float >("-channel_band"));
			obs.setMaxFreq(obs.getMinFreq() + ((obs.getNrChannels() - 1) * obs.getChannelBandwidth()));
		} else {
			std::cerr << "Need to specify the -header and -data arguments." << std::endl;
			throw std::exception();
		}
		outputFile = args.getSwitchArgument< std::string >("-output");

    obs.setDMRange(args.getSwitchArgument< unsigned int >("-dm_node"), args.getSwitchArgument< float >("-dm_first") + ((world.rank() / MPIRows) * obs.getNrDMs() * obs.getDMStep()), args.getSwitchArgument< float >("-dm_step"));
    obs.setPeriodRange(args.getSwitchArgument< unsigned int >("-period_node"), args.getSwitchArgument< unsigned int >("-period_first") + ((world.rank() % MPICols) * obs.getPeriodStep() * obs.getNrPeriods()), args.getSwitchArgument< unsigned int >("-period_step"));
		obs.setNrBins(args.getSwitchArgument< unsigned int >("-period_bins"));
	} catch ( isa::utils::EmptyCommandLine & err ) {
    // TODO: usage string
    std::cerr << std::endl;
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

	// Initialize OpenCL
	cl::Context * clContext = new cl::Context();
	std::vector< cl::Platform > * clPlatforms = new std::vector< cl::Platform >();
	std::vector< cl::Device > * clDevices = new std::vector< cl::Device >();
	std::vector< std::vector< cl::CommandQueue > > * clQueues = new std::vector< std::vector < cl::CommandQueue > >();

	try {
		initializeOpenCL(clPlatformID, 1, clPlatforms, clContext, clDevices, clQueues);
	} catch ( isa::OpenCL::OpenCLError & err ) {
		std::cerr << err.what() << std::endl;
		return 1;
	}

	// Host memory allocation
  std::vector< unsigned int > * shifts = PulsarSearch::getShifts(obs);
  obs.setNrSamplesPerDispersedChannel(obs.getNrSamplesPerSecond() + (*shifts)[((obs.getNrDMs() - 1) * obs.getNrPaddedChannels())]);
  secondsToBuffer = std::ceil(obs.getNrSamplesPerDispersedChannel() / static_cast< float >(obs.getNrSamplesPerSecond()));
  remainingSamples = obs.getNrSamplesPerDispersedChannel() % obs.getNrSamplesPerSecond();
  std::vector< unsigned int > * nrSamplesPerBin = PulsarSearch::getSamplesPerBin(obs);
  std::vector< dataType > dispersedData(obs.getNrChannels() * obs.getNrSamplesPerDispersedChannel());
  std::vector< float > snrTable(obs.getNrPeriods() * obs.getNrPaddedDMs());

  // Device memory allocation and data transfers
  cl::Buffer shifts_d, nrSamplesPerBin_d, dispersedData_d, dedispersedData_d, transposedData_d, foldedData_d, counterData0_d, counterData1_d, snrTable_d;

  try {
    shifts_d = cl::Buffer(*clContext, CL_MEM_READ_ONLY, shifts->size() * sizeof(unsigned int), 0, 0);
    samplesPerBin_d = cl::Buffer(*clContext, CL_MEM_READ_ONLY, samplesPerBin->size() * sizeof(unsigned int), 0, 0);
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
    delete samplesPerBin;
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

	if ( DEBUG && world.rank() == 0 ) {
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
  cl::Kernel * dedispersionK, * foldingK, transposeK, snrK;

  code = PulsarSearch::getDedispersionOpenCL< dataType >(dedispersionParameters[deviceName][obs.getNrDMs()][0], dedispersionParameters[deviceName][obs.getNrDMs()][1], dedispersionParameters[deviceName][obs.getNrDMs()][2], dedispersionParameters[deviceName][obs.getNrDMs()][3], dedispersionParameters[deviceName][obs.getNrDMs()][4], typeName, obs, *shifts);
	try {
    dedispersionK = isa::OpenCL::compile("dedispersion", *code, "-cl-mad-enable -Werror", *clContext, clDevices->at(clDeviceID));
	} catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
		return 1;
	}
  delete shifts;
  code = isa::OpenCL::getTransposeOpenCL(transposeParameters[deviceName][obs.getNrDMs()], obs.getNrSamplesPerSecond(), obs.getNrDMs(), obs.getPadding(), vectorWidth[deviceName], typeName);
  try {
    transposeK = isa::OpenCL::compile("transpose", *code, "-cl-mad-enable -Werror", *clContext, clDevices->at(clDeviceID));
  } catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }
  code = PulsarSearch::getFoldingOpenCL(foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][1], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][2], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][3], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][4], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][5], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][6], typeName, obs);
  try {
    foldingK = isa::OpenCL::compile("folding", *code, "-cl-mad-enable -Werror", *clContext, clDevices->at(clDeviceID));
  } catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }
  code = PulsarSearch::getSNROpenCL(snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0], snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][1], snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][2], snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][3], typeName, obs);
  try {
    snrK = isa::OpenCL::compile("snr", *code, "-cl-mad-enable -Werror", *clContext, clDevices->at(clDeviceID));
  } catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }

  // Set execution parameters
  cl::NDRange dedispersionGlobal(obs.getNrSamplesPerPaddedSecond() / dedispersionParameters[deviceName][obs.getNrDMs()][3], obs.getNrDMs() / dedispersionParameters[deviceName][obs.getNrDMs()][3]);
  cl::NDRange dedispersionLocal(dedispersionParameters[deviceName][obs.getNrDMs()][1], dedispersionParameters[deviceName][obs.getNrDMs()][2]);
  cl::NDRange transposeGlobal(obs.getNrSamplesPerSecond(), std::ceil(static_cast< double >(obs.getNrDMs()) / transposeParameters[deviceName][obs.getNrDMs()]));
  cl::NDRange transposeLocal(transposeParameters[deviceName][obs.getNrDMs()], 1);
  cl::NDRange foldingGlobal(obs.getNrPaddedDMs() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][6] / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][3], obs.getNrPeriods() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][6] / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][4], obs.getNrBins() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][6] / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][5]);
  cl::NDRange foldingLocal(foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][1], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][2]);
  cl::NDRange snrGlobal(obs.getNrPaddedDMs() / snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][2], obs.getNrPeriods() / snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][3]);
  cl::NDRange snrLocal(snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0], snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][1]);

  dedispersionK->setArg(0, dispersedData_d);
  dedispersionK->setArg(1, dedispersedData_d);
  dedispersionK->setArg(2, shifts_d);
  transposeK->setArg(0, dedispersedData_d);
  transposeK->setArg(1, transposedData_d);
  foldingK->setArg(1, transposedData_d);
  foldingK->setArg(2, foldedData_d);
  foldingK->setArg(5, samplesPerBin_d);
  snrK->setArg(0, foldedData_d);
  snrK->setArg(1, snrTable_d);

	// Search loop
  cl::Event syncPoint;
  isa::utils::Timer searchTime;
  isa::utils::Timer inputLoadTime;
  isa::utils::Timer dedispTime;
  isa::utils::Timer transTime;
  isa::utils::Timer foldTime;
  isa::utils::Timer snrTime;
  isa::utils::Timer outputStoreTime;

  searchTime.start();
	if ( DEBUG && world.rank() == 0 ) {
		std::cout << "Starting the search." << std::endl;
		std::cout << "Processing seconds: ";
	}
	for ( unsigned int second = 0; second <= obs.getNrSeconds() - secondsToBuffer; second++ ) {
		if ( DEBUG && world.rank() == 0 ) {
			std::cout << second << " " << std::flush;
		}
		// Load the input
		inputLoadTime.start();
		for ( unsigned int channel = 0; channel < obs.getNrChannels(); channel++ ) {
			for ( unsigned int chunk = 0; chunk < secondsToBuffer - 1; chunk++ ) {
        memcpy(reinterpret_cast< void * >(&(dispersedData.data()[(channel * obs.getNrSamplesPerDispersedChannel()) + (chunk * obs.getNrSamplesPerSecond())])), (input->at(second + chunk))[channel * obs.getNrSamplesPerPaddedSecond()], obs.getNrSamplesPerSecond() * sizeof(dataType));
			}
      memcpy(reinterpret_cast< void * >(&(dispersedData.data()[(channel * obs.getNrSamplesPerDispersedChannel()) + ((secondsToBuffer - 1) * obs.getNrSamplesPerSecond())])), (input->at(second + (secondsToBuffer - 1)))[channel * obs.getNrSamplesPerPaddedSecond()], remainingSamples * sizeof(dataType));
		}
    try {
      clQueues->at(clDeviceID)[0].enqueueWriteBuffer(dispersedData_d, CL_TRUE, 0, dispersedData.size() * sizeof(dataType), reinterpret_cast< void * >(dispersedData.data()));
    } catch ( cl::Error & err ) {
      std::cerr << err.what() << std::endl;
      return 1;
    }
		inputLoadTime.stop();

		// Run the kernels
		try {
      dedispTime.start();
      clQueues->at(clDeviceID)[0].enqueueNDRangeKernel(*dedispersionK, cl::NullRange, dedispersionGlobal, dedispersionLocal, 0, &syncPoint);
      syncPoint.wait();
      dedispTime.stop();
      transTime.start();
      clQueues->at(clDeviceID)[0].enqueueNDRangeKernel(*transposeK, cl::NullRange, transposeGlobal, transposeLocal, 0, &syncPoint);
      syncPoint.wait();
      transTime.stop();
      foldTime.start();
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
      foldTime.stop();
		} catch ( cl::Error & err ) {
			std::cerr << err.what() << std::endl;
			return 1;
		}
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
	try {
    snrTime.start();
    clQueues->at(clDeviceID)[0].enqueueNDRangeKernel(*snrK, cl::NullRange, snrGlobal, snrLocal, 0, &syncPoint);
    syncPoint.wait();
    snrTime.stop();
    outputStoreTime.start();
    clQueues->at(clDeviceID)[0].enqueueReadBuffer(snrTable_d, CL_TRUE, 0, snrTable.size() * sizeof(float), reinterpret_cast< void * >(snrTable.data()));
    outputStoreTime.stop();
  } catch ( cl::Error & err ) {
		std::cerr << err.what() << std::endl;
		return 1;
	}
	if ( DEBUG && world.rank() == 0 ) {
		std::cout << "Saving output to disk." << std::endl;
	}
	std::ofstream output;
  output << std::fixed << std::setprecision(3);
	output.open(outputFile + "_" + isa::utils::toString(world.rank()));
	for ( unsigned int period = 0; period < obs.getNrPeriods(); period++ ) {
		for ( unsigned int dm = 0; dm < obs.getNrDMs(); dm++ ) {
			output << (world.rank() * obs.getNrPeriods()) + period << " ";
      output << (obs.getFirstPeriod() + (period * obs.getPeriodStep())) / static_cast< float > (obs.getNrSamplesPerSecond()) << " ";
      output << dm << " ";
      output << obs.getFirstDM() + (dm * obs.getDMStep()) << " ";
      output << snrTable[(period * obs.getNrPaddedDMs()) + dm] << std::endl;
		}
	}
	output.close();

	// Wait for all MPI processes
	double maxTime = 0.0;
	double maxKernel = 0.0;
	std::vector< double > nodeSearchTimes(world.size());
	std::vector< double > nodeMainLoopTimes(world.size());
	gather(world, searchTime.getTotalTime() + outputTime.getTotalTime(), nodeSearchTimes, 0);
	gather(world, (clDedisperse.getTimer()).getTotalTime() + (clFold.getTimer()).getTotalTime(), nodeMainLoopTimes, 0);
	for ( int node = 0; node < world.size(); node++ ) {
		if ( nodeSearchTimes[node] > maxTime ) {
			maxTime = nodeSearchTimes[0];
		}
		if ( nodeMainLoopTimes[node] > maxKernel ) {
			maxKernel = nodeMainLoopTimes[node];
		}
	}

	if ( world.rank() == 0 ) {
		std::cout << " " << std::endl;
	}

	return 0;
}

