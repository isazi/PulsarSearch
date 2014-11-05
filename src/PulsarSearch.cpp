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
  unsigned int MPIRows = 0;
  unsigned int MPICols = 0;
	unsigned int bytesToSkip = 0;
  unsigned int secondsToBuffer = 0;
  unsigned int remainingSamples = 0;
	std::string deviceName;
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
    std::cerr <<  args.getName() << " -mpi_cols ... -mpi_rows ... -opencl_platform ... -opencl_device ... -device_name ... -padding_file ... -vector_file ... -dedispersion_file ... -transpose_file ... -folding_file ... -snr_file ... [-lofar] [-sigproc] -output ... -dm_node ... -dm_first ... -dm_step ... -period_node ... -period_first ... -period_step ... -period_bins ..."<< std::endl;
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
    std::cout << "Device: " << deviceName << std::endl;
    std::cout << "Padding: " << padding[deviceName] << std::endl;
    std::cout << "Vector: " << vectorWidth[deviceName] << std::endl;
    std::cout << std::endl;
    std::cout << "Seconds: " << obs.getNrSeconds() << std::endl;
    std::cout << "Samples: " << obs.getNrSamplesPerSecond() << std::endl;
    std::cout << "Frequency range: " << obs.getMinFreq() << " MHz, " << obs.getMaxFreq() << " MHz" << std::endl;
    std::cout << "Channels: " << obs.getNrChannels() << " (" << obs.getChannelBandwidth() << " MHz)" << std::endl;
    std::cout << std::endl;
		std::cout << "Time to load the input: " << std::fixed << std::setprecision(6) << loadTime.getTotalTime() << " seconds." << std::endl;
    std::cout << std::endl;
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
  std::vector< dataType > maxDedispersedTable(obs.getNrPaddedDMs());
  std::vector< float > meanDedispersedTable(obs.getNrPaddedDMs());
  std::vector< float > rmsDedispersedTable(obs.getNrPaddedDMs());
  std::vector< float > snrFoldedTable(obs.getNrPeriods() * obs.getNrPaddedDMs());

  // Device memory allocation and data transfers
  cl::Buffer shifts_d, nrSamplesPerBin_d, dispersedData_d, dedispersedData_d, transposedData_d, foldedData_d, counterData0_d, counterData1_d, maxDedispersedTable_d, meanDedispersedTable_d, rmsDedispersedTable_d, snrFoldedTable_d;

  try {
    shifts_d = cl::Buffer(*clContext, CL_MEM_READ_ONLY, shifts->size() * sizeof(unsigned int), 0, 0);
    nrSamplesPerBin_d = cl::Buffer(*clContext, CL_MEM_READ_ONLY, nrSamplesPerBin->size() * sizeof(unsigned int), 0, 0);
    dispersedData_d = cl::Buffer(*clContext, CL_MEM_READ_ONLY, dispersedData.size() * sizeof(dataType), 0, 0);
    dedispersedData_d = cl::Buffer(*clContext, CL_MEM_READ_WRITE, obs.getNrDMs() * obs.getNrSamplesPerPaddedSecond() * sizeof(dataType), 0, 0);
    transposedData_d = cl::Buffer(*clContext, CL_MEM_READ_WRITE, obs.getNrSamplesPerSecond() * obs.getNrPaddedDMs() * sizeof(dataType), 0, 0);
    foldedData_d = cl::Buffer(*clContext, CL_MEM_READ_WRITE, obs.getNrBins() * obs.getNrPeriods() * obs.getNrPaddedDMs() * sizeof(dataType), 0, 0);
    counterData0_d = cl::Buffer(*clContext, CL_MEM_READ_WRITE, obs.getNrBins() * obs.getNrPaddedPeriods() * sizeof(unsigned int), 0, 0);
    counterData1_d = cl::Buffer(*clContext, CL_MEM_READ_WRITE, obs.getNrBins() * obs.getNrPaddedPeriods() * sizeof(unsigned int), 0, 0);
    maxDedispersedTable_d = cl::Buffer(*clContext, CL_MEM_READ_WRITE, obs.getNrPaddedDMs() * sizeof(dataType), 0, 0);
    meanDedispersedTable_d = cl::Buffer(*clContext, CL_MEM_READ_WRITE, obs.getNrPaddedDMs() * sizeof(float), 0, 0);
    rmsDedispersedTable_d = cl::Buffer(*clContext, CL_MEM_READ_WRITE, obs.getNrPaddedDMs() * sizeof(float), 0, 0);
    snrFoldedTable_d = cl::Buffer(*clContext, CL_MEM_WRITE_ONLY, obs.getNrPeriods() * obs.getNrPaddedDMs() * sizeof(float), 0, 0);

    // shifts_d
    clQueues->at(clDeviceID)[0].enqueueWriteBuffer(shifts_d, CL_TRUE, 0, shifts->size() * sizeof(unsigned int), reinterpret_cast< void * >(shifts->data()));
    // nrSamplesPerBin_d
    clQueues->at(clDeviceID)[0].enqueueWriteBuffer(nrSamplesPerBin_d, CL_TRUE, 0, nrSamplesPerBin->size() * sizeof(unsigned int), reinterpret_cast< void * >(nrSamplesPerBin->data()));
    delete nrSamplesPerBin;
    // foldedData_d and maxDedispersedTable_d
    std::vector< dataType > transferDataType(obs.getNrBins() * obs.getNrPeriods() * obs.getNrPaddedDMs());
    std::fill(transferDataType.begin(), transferDataType.end(), 0);
    clQueues->at(clDeviceID)[0].enqueueWriteBuffer(foldedData_d, CL_TRUE, 0, transferDataType.size() * sizeof(dataType), reinterpret_cast< void * >(transferDataType.data()));
    // maxDedispersedTable_d
    std::fill(maxDedispersedTable.begin(), maxDedispersedTable.end(), 0);
    clQueues->at(clDeviceID)[0].enqueueWriteBuffer(maxDedispersedTable_d, CL_TRUE, 0, maxDedispersedTable.size() * sizeof(dataType), reinterpret_cast< void * >(maxDedispersedTable.data()));
    // meanDedispersedTable_d and rmsDedispersedTable_d
    std::fill(meanDedispersedTable.begin(), meanDedispersedTable.end(), 0);
    clQueues->at(clDeviceID)[0].enqueueWriteBuffer(meanDedispersedTable_d, CL_TRUE, 0, meanDedispersedTable.size() * sizeof(float), reinterpret_cast< void * >(meanDedispersedTable.data()));
    std::fill(rmsDedispersedTable.begin(), rmsDedispersedTable.end(), 0);
    clQueues->at(clDeviceID)[0].enqueueWriteBuffer(rmsDedispersedTable_d, CL_TRUE, 0, rmsDedispersedTable.size() * sizeof(float), reinterpret_cast< void * >(rmsDedispersedTable.data()));
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
    hostMemory += maxDedispersedTable.size() * sizeof(dataType);
    hostMemory += meanDedispersedTable.size() * sizeof(float);
    hostMemory += rmsDedispersedTable.size() * sizeof(float);
    hostMemory += snrFoldedTable.size() * sizeof(float);
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
  cl::Kernel * dedispersionK, * foldingK, * transposeK, * snrDedispersedK, * snrFoldedK;

  code = PulsarSearch::getDedispersionOpenCL(dedispersionParameters[deviceName][obs.getNrDMs()][0], dedispersionParameters[deviceName][obs.getNrDMs()][1], dedispersionParameters[deviceName][obs.getNrDMs()][2], dedispersionParameters[deviceName][obs.getNrDMs()][3], dedispersionParameters[deviceName][obs.getNrDMs()][4], dataName, obs, *shifts);
	try {
    dedispersionK = isa::OpenCL::compile("dedispersion", *code, "-cl-mad-enable -Werror", *clContext, clDevices->at(clDeviceID));
	} catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
		return 1;
	}
  delete shifts;
  code = isa::OpenCL::getTransposeOpenCL(transposeParameters[deviceName][obs.getNrDMs()], obs.getNrDMs(), obs.getNrSamplesPerSecond(), obs.getPadding(), vectorWidth[deviceName], dataName);
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
  code = PulsarSearch::getSNRDedispersedOpenCL(snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0], snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][2], dataName, obs);
  try {
    snrDedispersedK = isa::OpenCL::compile("snrDedispersed", *code, "-cl-mad-enable -Werror", *clContext, clDevices->at(clDeviceID));
  } catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }
  code = PulsarSearch::getSNRFoldedOpenCL(snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0], snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][1], snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][2], snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][3], dataName, obs);
  try {
    snrFoldedK = isa::OpenCL::compile("snrFolded", *code, "-cl-mad-enable -Werror", *clContext, clDevices->at(clDeviceID));
  } catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }

  // Set execution parameters
  cl::NDRange dedispersionGlobal(obs.getNrSamplesPerPaddedSecond() / dedispersionParameters[deviceName][obs.getNrDMs()][3], obs.getNrDMs() / dedispersionParameters[deviceName][obs.getNrDMs()][3]);
  cl::NDRange dedispersionLocal(dedispersionParameters[deviceName][obs.getNrDMs()][1], dedispersionParameters[deviceName][obs.getNrDMs()][2]);
  cl::NDRange transposeGlobal(obs.getNrSamplesPerPaddedSecond(), static_cast< unsigned int >(std::ceil(static_cast< double >(obs.getNrDMs()) / transposeParameters[deviceName][obs.getNrDMs()])));
  cl::NDRange transposeLocal(transposeParameters[deviceName][obs.getNrDMs()], 1);
  cl::NDRange foldingGlobal(obs.getNrPaddedDMs() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][6] / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][3], obs.getNrPeriods() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][6] / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][4], obs.getNrBins() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][6] / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][5]);
  cl::NDRange foldingLocal(foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][1], foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][2]);
  cl::NDRange snrDedispersedGlobal(obs.getNrPaddedDMs() / snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][2]);
  cl::NDRange snrDedispersedLocal(snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0]);
  cl::NDRange snrFoldedGlobal(obs.getNrPaddedDMs() / snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][2], obs.getNrPeriods() / snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][3]);
  cl::NDRange snrFoldedLocal(snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0], snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][1]);

  dedispersionK->setArg(0, dispersedData_d);
  dedispersionK->setArg(1, dedispersedData_d);
  dedispersionK->setArg(2, shifts_d);
  transposeK->setArg(0, dedispersedData_d);
  transposeK->setArg(1, transposedData_d);
  foldingK->setArg(1, transposedData_d);
  foldingK->setArg(2, foldedData_d);
  foldingK->setArg(5, nrSamplesPerBin_d);
  snrDedispersedK->setArg(1, transposedData_d);
  snrDedispersedK->setArg(2, maxDedispersedTable_d);
  snrDedispersedK->setArg(3, meanDedispersedTable_d);
  snrDedispersedK->setArg(4, rmsDedispersedTable_d);
  snrFoldedK->setArg(0, foldedData_d);
  snrFoldedK->setArg(1, snrFoldedTable_d);

	// Search loop
  cl::Event syncPoint;
  isa::utils::Timer searchTime;
  isa::utils::Timer inputLoadTime;
  isa::utils::Timer dedispTime;
  isa::utils::Timer transTime;
  isa::utils::Timer foldTime;
  isa::utils::Timer snrDedispersedTime;
  isa::utils::Timer snrFoldedTime;
  isa::utils::Timer outputStoreTime;

  world.barrier();
  searchTime.start();
	if ( DEBUG && world.rank() == 0 ) {
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
      snrDedispersedTime.start();
      snrDedispersedK->setArg(0, static_cast< float >(second));
      clQueues->at(clDeviceID)[0].enqueueNDRangeKernel(*snrDedispersedK, cl::NullRange, snrDedispersedGlobal, snrDedispersedLocal, 0, &syncPoint);
      syncPoint.wait();
      snrDedispersedTime.stop();
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
    snrFoldedTime.start();
    clQueues->at(clDeviceID)[0].enqueueNDRangeKernel(*snrFoldedK, cl::NullRange, snrFoldedGlobal, snrFoldedLocal, 0, &syncPoint);
    syncPoint.wait();
    snrFoldedTime.stop();
    outputStoreTime.start();
    clQueues->at(clDeviceID)[0].enqueueReadBuffer(maxDedispersedTable_d, CL_TRUE, 0, maxDedispersedTable.size() * sizeof(dataType), reinterpret_cast< void * >(maxDedispersedTable.data()));
    clQueues->at(clDeviceID)[0].enqueueReadBuffer(meanDedispersedTable_d, CL_TRUE, 0, meanDedispersedTable.size() * sizeof(float), reinterpret_cast< void * >(meanDedispersedTable.data()));
    clQueues->at(clDeviceID)[0].enqueueReadBuffer(rmsDedispersedTable_d, CL_TRUE, 0, rmsDedispersedTable.size() * sizeof(float), reinterpret_cast< void * >(rmsDedispersedTable.data()));
    clQueues->at(clDeviceID)[0].enqueueReadBuffer(snrFoldedTable_d, CL_TRUE, 0, snrFoldedTable.size() * sizeof(float), reinterpret_cast< void * >(snrFoldedTable.data()));
    clQueues->at(clDeviceID)[0].enqueueReadBuffer(foldedData_d, CL_TRUE, 0, foldedData.size() * sizeof(dataType), reinterpret_cast< void * >(foldedData.data()));
    outputStoreTime.stop();
  } catch ( cl::Error & err ) {
		std::cerr << err.what() << std::endl;
		return 1;
	}

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
      output << snrFoldedTable[(period * obs.getNrPaddedDMs()) + dm] << std::endl;
		}
	}
	output.close();
  output.open(outputFile + "_" + isa::utils::toString(world.rank()) + ".fold");
  output << std::fixed << std::setprecision(6);
  for ( unsigned int dm = 0; dm < obs.getNrDMs(); dm++ ) {
    for ( unsigned int period = 0; period < obs.getNrPeriods(); period++ ) {
      output << "# " << dm << " " << period << std::endl;
      for ( unsigned int bin = 0; bin < obs.getNrBins(); bin++ ) {
        output << bin << " " << foldedData[(bin * obs.getNrPeriods() * obs.getNrPaddedDMs()) + (period * obs.getNrPaddedDMs()) + dm] << std::endl;
      }
      output << std::endl << std:: endl;
    }
  }
  output.close();
	output.open(outputFile + "_" + isa::utils::toString(world.rank()) + ".dedi");
  output << "# DM SNR" << std::endl;
  for ( unsigned int dm = 0; dm < obs.getNrDMs(); dm++ ) {
    output << dm << " " << (maxDedispersedTable[dm] - meanDedispersedTable[dm]) / std::sqrt(rmsDedispersedTable[dm]) << std::endl;
  }
  output.close();
	output.open(outputFile + "_" + isa::utils::toString(world.rank()) + ".stat");
  output << "# searchTime inputLoadTotal inputLoadAvg err dedispersionTotal dedispersionvg err transposeTotal transposeAvg err snrDedispersedTotal snrDedispersedAvg err foldingTotal foldingAvg err snrFoldedTotal snrFoldedAvg err outputStoreTotal outputStoreAvg err" << std::endl;
  output << std::fixed << std::setprecision(6);
  output << searchTime.getTotalTime() << " ";
  output << inputLoadTime.getTotalTime() << " " << inputLoadTime.getAverageTime() << " " << inputLoadTime.getStandardDeviation() << " ";
  output << dedispTime.getTotalTime() << " " << dedispTime.getAverageTime() << " " << dedispTime.getStandardDeviation() << " ";
  output << transTime.getTotalTime() << " " << transTime.getAverageTime() << " " << transTime.getStandardDeviation() << " ";
  output << snrDedispersedTime.getTotalTime() << " " << snrDedispersedTime.getAverageTime() << " " << snrDedispersedTime.getStandardDeviation() << " ";
  output << foldTime.getTotalTime() << " " << foldTime.getAverageTime() << " " << foldTime.getStandardDeviation() << " ";
  output << snrFoldedTime.getTotalTime() << " " << snrFoldedTime.getAverageTime() << " " << snrFoldedTime.getStandardDeviation() << " ";
  output << outputStoreTime.getTotalTime() << " " << outputStoreTime.getAverageTime() << " " << outputStoreTime.getStandardDeviation() << " ";
  output << std::endl;
  output.close();

  if ( DEBUG && world.rank() == 0 ) {
    std::cout << "Output and statistics saved to disk." << std::endl;
  }
  world.barrier();

	return 0;
}

