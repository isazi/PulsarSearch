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
#include <cmath>

#include <configuration.hpp>

#include <ArgumentList.hpp>
#include <utils.hpp>
#include <Timer.hpp>
#include <Observation.hpp>
#include <Platform.hpp>
#include <ReadData.hpp>
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
  bool print = false;
  bool saveOutput = false;
  bool noData = false;
  bool random = false;
	bool dataLOFAR = false;
	bool dataSIGPROC = false;
  bool limit = false;
	unsigned int clPlatformID = 0;
	unsigned int clDeviceID = 0;
	unsigned int bytesToSkip = 0;
  unsigned int secondsToBuffer = 0;
  unsigned int nrThreads = 0;
  unsigned int remainingSamples = 0;
	std::string deviceName;
	std::string dataFile;
	std::string headerFile;
	std::string outputFile;
	std::ofstream output;
  isa::utils::ArgumentList args(argc, argv);
	// Observation object
  AstroData::Observation obs;
  // Fake pulsar
  unsigned int period = 0;
  unsigned int width = 0;
  float DM = 0;
  // Configurations
  AstroData::paddingConf padding;
  AstroData::vectorWidthConf vectorWidth;
  PulsarSearch::tunedDedispersionConf dedispersionParameters;
  isa::OpenCL::tunedTransposeConf transposeParameters;
  PulsarSearch::tunedFoldingConf foldingParameters;
  PulsarSearch::tunedSNRDedispersedConf snrDParameters;
  PulsarSearch::tunedSNRFoldedConf snrFParameters;

	try {
		clPlatformID = args.getSwitchArgument< unsigned int >("-opencl_platform");
		clDeviceID = args.getSwitchArgument< unsigned int >("-opencl_device");
		deviceName = args.getSwitchArgument< std::string >("-device_name");

    AstroData::readPaddingConf(padding, args.getSwitchArgument< std::string >("-padding_file"));
    AstroData::readVectorWidthConf(vectorWidth, args.getSwitchArgument< std::string >("-vector_file"));
    PulsarSearch::readTunedDedispersionConf(dedispersionParameters, args.getSwitchArgument< std::string >("-dedispersion_file"));
    isa::OpenCL::readTunedTransposeConf(transposeParameters, args.getSwitchArgument< std::string >("-transpose_file"));
    PulsarSearch::readTunedFoldingConf(foldingParameters, args.getSwitchArgument< std::string >("-folding_file"));
    PulsarSearch::readTunedSNRDedispersedConf(snrDParameters, args.getSwitchArgument< std::string >("-snrd_file"));
		PulsarSearch::readTunedSNRFoldedConf(snrFParameters, args.getSwitchArgument< std::string >("-snrf_file"));

    print = args.getSwitch("-print");
    saveOutput = args.getSwitch("-save_output");
		obs.setPadding(padding[deviceName]);

		dataLOFAR = args.getSwitch("-lofar");
		dataSIGPROC = args.getSwitch("-sigproc");
		if ( dataLOFAR && dataSIGPROC ) {
			std::cerr << "-lofar and -sigproc are mutually exclusive." << std::endl;
			throw std::exception();
		} else if ( dataLOFAR ) {
			headerFile = args.getSwitchArgument< std::string >("-header");
			dataFile = args.getSwitchArgument< std::string >("-data");
      limit = args.getSwitch("-limit");
      if ( limit ) {
        obs.setNrSeconds(args.getSwitchArgument< unsigned int >("-seconds"));
      }
		} else if ( dataSIGPROC ) {
			bytesToSkip = args.getSwitchArgument< unsigned int >("-header");
			dataFile = args.getSwitchArgument< std::string >("-data");
			obs.setNrSeconds(args.getSwitchArgument< unsigned int >("-seconds"));
      obs.setFrequencyRange(args.getSwitchArgument< unsigned int >("-channels"), args.getSwitchArgument< float >("-min_freq"), args.getSwitchArgument< float >("-channel_bandwidth"));
			obs.setNrSamplesPerSecond(args.getSwitchArgument< unsigned int >("-samples"));
		} else {
      noData = args.getSwitch("-no_data");
      if ( !noData ) {
        random = args.getSwitch("-random");
        period = args.getSwitchArgument< unsigned int >("-period");
        width = args.getSwitchArgument< unsigned int >("-width");
        DM = args.getSwitchArgument< float >("-dm");
      }
      obs.setNrSeconds(args.getSwitchArgument< unsigned int >("-seconds"));
      obs.setFrequencyRange(args.getSwitchArgument< unsigned int >("-channels"), args.getSwitchArgument< float >("-min_freq"), args.getSwitchArgument< float >("-channel_bandwidth"));
      obs.setNrSamplesPerSecond(args.getSwitchArgument< unsigned int >("-samples"));
		}
		outputFile = args.getSwitchArgument< std::string >("-output");

    unsigned int tempUInts[3] = {args.getSwitchArgument< unsigned int >("-dm_node"), 0, 0};
    float tempFloats[2] = {args.getSwitchArgument< float >("-dm_first"), args.getSwitchArgument< float >("-dm_step")};
    obs.setDMRange(tempUInts[0], tempFloats[0], tempFloats[1]);
    tempUInts[0] = args.getSwitchArgument< unsigned int >("-period_node");
    tempUInts[1] = args.getSwitchArgument< unsigned int >("-period_first");
    tempUInts[2] = args.getSwitchArgument< unsigned int >("-period_step");
    obs.setPeriodRange(tempUInts[0], tempUInts[1], tempUInts[2]);
		obs.setNrBins(args.getSwitchArgument< unsigned int >("-period_bins"));
	} catch ( isa::utils::EmptyCommandLine & err ) {
    std::cerr <<  args.getName() << " -opencl_platform ... -opencl_device ... -device_name ... -padding_file ... -vector_file ... -dedispersion_file ... -transpose_file ... -folding_file ... -snrd_file ... -snrf_file [-print] [-save_output] [-lofar] [-sigproc] -output ... -dm_node ... -dm_first ... -dm_step ... -period_node ... -period_first ... -period_step ... -period_bins ..."<< std::endl;
    std::cerr << "\t -lofar -header ... -data ... [-limit]" << std::endl;
    std::cerr << "\t\t -limit -seconds ..." << std::endl;
    std::cerr << "\t -sigproc -header ... -data ... -seconds ... -channels ... -min_freq ... -channel_bandwidth ... -samples ..." << std::endl;
    std::cerr << "\t [-random] -period ... -width ... -dm ... -seconds ... -channels ... -min_freq ... -channel_bandwidth ... -samples ..." << std::endl;
    std::cerr << "\t -no_data -seconds ... -channels ... -min_freq ... -channel_bandwidth ... -samples ..." << std::endl;
    return 1;
  } catch ( std::exception & err ) {
		std::cerr << err.what() << std::endl;
		return 1;
	}

	// Load observation data
  isa::utils::Timer loadTime;
	std::vector< std::vector< dataType > * > * input = new std::vector< std::vector< dataType > * >(obs.getNrSeconds());
	if ( dataLOFAR ) {
    loadTime.start();
    if ( limit ) {
      AstroData::readLOFAR(headerFile, dataFile, obs, *input, obs.getNrSeconds());
    } else {
      AstroData::readLOFAR(headerFile, dataFile, obs, *input);
    }
    loadTime.stop();
	} else if ( dataSIGPROC ) {
    loadTime.start();
		input->resize(obs.getNrSeconds());
    AstroData::readSIGPROC(obs, bytesToSkip, dataFile, *input);
    loadTime.stop();
	} else {
    if ( noData ) {
      input->at(0) = new std::vector< dataType >(obs.getNrChannels() * obs.getNrSamplesPerPaddedSecond());
      std::fill(input->at(0)->begin(), input->at(0)->end(), 42);
      for ( unsigned int second = 1; second < obs.getNrSeconds(); second++ ) {
        input->at(second) = input->at(0);
      }
    } else {
      AstroData::generatePulsar(period, width, DM, obs, *input, random);
    }
  }
	if ( DEBUG == 0 ) {
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
  std::vector< float > * shifts = PulsarSearch::getShifts(obs);
  obs.setNrSamplesPerDispersedChannel(obs.getNrSamplesPerSecond() + static_cast< unsigned int >(shifts->at(0) * (obs.getFirstDM() + ((obs.getNrDMs() - 1) * obs.getDMStep()))));
  secondsToBuffer = obs.getNrSamplesPerDispersedChannel() / obs.getNrSamplesPerSecond();
  remainingSamples = obs.getNrSamplesPerDispersedChannel() % obs.getNrSamplesPerSecond();
  std::vector< unsigned int > * nrSamplesPerBin = PulsarSearch::getSamplesPerBin(obs);
  std::vector< dataType > dispersedData(obs.getNrChannels() * obs.getNrSamplesPerDispersedChannel());
  std::vector< dataType > dedispersedData(obs.getNrDMs() * obs.getNrSamplesPerPaddedSecond());
  std::vector< dataType > transposedData(obs.getNrSamplesPerSecond() * obs.getNrPaddedDMs());
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
    // maxDedispersedTable_d, meanDedispersedTable_d, rmsDedispersedTable_d
    std::fill(maxDedispersedTable.begin(), maxDedispersedTable.end(), 0);
    clQueues->at(clDeviceID)[0].enqueueWriteBuffer(maxDedispersedTable_d, CL_TRUE, 0, maxDedispersedTable.size() * sizeof(dataType), reinterpret_cast< void * >(maxDedispersedTable.data()));
    std::fill(meanDedispersedTable.begin(), meanDedispersedTable.end(), 0);
    clQueues->at(clDeviceID)[0].enqueueWriteBuffer(meanDedispersedTable_d, CL_TRUE, 0, meanDedispersedTable.size() * sizeof(float), reinterpret_cast< void * >(meanDedispersedTable.data()));
    std::fill(rmsDedispersedTable.begin(), rmsDedispersedTable.end(), 0);
    clQueues->at(clDeviceID)[0].enqueueWriteBuffer(rmsDedispersedTable_d, CL_TRUE, 0, rmsDedispersedTable.size() * sizeof(float), reinterpret_cast< void * >(rmsDedispersedTable.data()));
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

	if ( DEBUG == 0 ) {
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
    std::cout << std::endl;
	}

	// Generate OpenCL kernels
  std::string * code;
  cl::Kernel * dedispersionK, * foldingK, * transposeK, * snrDedispersedK, * snrFoldedK;

  code = PulsarSearch::getDedispersionOpenCL(dedispersionParameters[deviceName][obs.getNrDMs()], dataName, obs, *shifts);
	try {
    dedispersionK = isa::OpenCL::compile("dedispersion", *code, "-cl-mad-enable -Werror", *clContext, clDevices->at(clDeviceID));
	} catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
		return 1;
	}
  delete shifts;
  delete code;
  code = isa::OpenCL::getTransposeOpenCL(transposeParameters[deviceName][obs.getNrDMs()], obs.getNrDMs(), obs.getNrSamplesPerSecond(), obs.getPadding(), vectorWidth[deviceName], dataName);
  try {
    transposeK = isa::OpenCL::compile("transpose", *code, "-cl-mad-enable -Werror", *clContext, clDevices->at(clDeviceID));
  } catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }
  delete code;
  code = PulsarSearch::getFoldingOpenCL(foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()], dataName, obs);
  try {
    foldingK = isa::OpenCL::compile("folding", *code, "-cl-mad-enable -Werror", *clContext, clDevices->at(clDeviceID));
  } catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }
  delete code;
  code = PulsarSearch::getSNRDedispersedOpenCL(snrDParameters[deviceName][obs.getNrDMs()], dataName, obs);
  try {
    snrDedispersedK = isa::OpenCL::compile("snrDedispersed", *code, "-cl-mad-enable -Werror", *clContext, clDevices->at(clDeviceID));
  } catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }
  delete code;
  code = PulsarSearch::getSNRFoldedOpenCL(snrFParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()], dataName, obs);
  try {
    snrFoldedK = isa::OpenCL::compile("snrFolded", *code, "-cl-mad-enable -Werror", *clContext, clDevices->at(clDeviceID));
  } catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }
  delete code;

  // Set execution parameters
  if ( obs.getNrSamplesPerSecond() % (dedispersionParameters[deviceName][obs.getNrDMs()].getNrSamplesPerBlock() * dedispersionParameters[deviceName][obs.getNrDMs()].getNrSamplesPerThread()) == 0 ) {
    nrThreads = obs.getNrSamplesPerSecond() / dedispersionParameters[deviceName][obs.getNrDMs()].getNrSamplesPerThread();
  } else {
    nrThreads = obs.getNrSamplesPerPaddedSecond() / dedispersionParameters[deviceName][obs.getNrDMs()].getNrSamplesPerThread();
  }
  cl::NDRange dedispersionGlobal(nrThreads, obs.getNrDMs() / dedispersionParameters[deviceName][obs.getNrDMs()].getNrDMsPerThread());
  cl::NDRange dedispersionLocal(dedispersionParameters[deviceName][obs.getNrDMs()].getNrSamplesPerBlock(), dedispersionParameters[deviceName][obs.getNrDMs()].getNrDMsPerBlock());
  if ( DEBUG == 0 ) {
    std::cout << "Dedispersion" << std::endl;
    std::cout << "Global: " << nrThreads << ", " << obs.getNrDMs() / dedispersionParameters[deviceName][obs.getNrDMs()].getNrDMsPerThread() << std::endl;
    std::cout << "Local: " << dedispersionParameters[deviceName][obs.getNrDMs()].getNrSamplesPerBlock() << ", " << dedispersionParameters[deviceName][obs.getNrDMs()].getNrDMsPerBlock() << std::endl;
    std::cout << "Parameters: ";
    std::cout << dedispersionParameters[deviceName][obs.getNrDMs()].print() << std::endl;
    std::cout << std::endl;
  }
  cl::NDRange transposeGlobal(obs.getNrPaddedDMs(), static_cast< unsigned int >(std::ceil(static_cast< double >(obs.getNrSamplesPerSecond()) / transposeParameters[deviceName][obs.getNrDMs()].getNrItemsPerBlock())));
  cl::NDRange transposeLocal(transposeParameters[deviceName][obs.getNrDMs()].getNrItemsPerBlock(), 1);
  if ( DEBUG == 0 ) {
    std::cout << "Transpose" << std::endl;
    std::cout << "Global: " << obs.getNrPaddedDMs() << ", " << static_cast< unsigned int >(std::ceil(static_cast< double >(obs.getNrSamplesPerSecond()) / transposeParameters[deviceName][obs.getNrDMs()].getNrItemsPerBlock())) << std::endl;
    std::cout << "Local: " << transposeParameters[deviceName][obs.getNrDMs()].getNrItemsPerBlock() << ", 1" << std::endl;
    std::cout << "Parameters: " << transposeParameters[deviceName][obs.getNrDMs()].print() << std::endl;
    std::cout << std::endl;
  }
  if ( obs.getNrDMs() % (foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrDMsPerBlock() * foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrDMsPerThread() * foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getVector()) == 0 ) {
    nrThreads = obs.getNrDMs() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrDMsPerThread() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getVector();
  } else {
    nrThreads = obs.getNrPaddedDMs() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrDMsPerThread() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getVector();
  }
  cl::NDRange foldingGlobal(nrThreads, obs.getNrPeriods() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrPeriodsPerThread(), obs.getNrBins() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrBinsPerThread());
  cl::NDRange foldingLocal(foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrDMsPerBlock(), foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrPeriodsPerBlock(), foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrBinsPerBlock());
  if ( DEBUG == 0 ) {
    std::cout << "Folding" << std::endl;
    std::cout << "Global: " << nrThreads << ", " << obs.getNrPeriods() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrPeriodsPerThread() << ", " << obs.getNrBins() / foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrBinsPerThread() << std::endl;
    std::cout << "Local: " << foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrDMsPerBlock() << ", " << foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrPeriodsPerBlock() << ", " << foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrBinsPerBlock() << std::endl;
    std::cout << "Parameters: ";
    std::cout << foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].print() << std::endl;
    std::cout << std::endl;
  }
  if ( obs.getNrDMs() % (snrDParameters[deviceName][obs.getNrDMs()].getNrDMsPerBlock() * snrDParameters[deviceName][obs.getNrDMs()].getNrDMsPerThread()) == 0 ) {
    nrThreads = obs.getNrDMs() / snrDParameters[deviceName][obs.getNrDMs()].getNrDMsPerThread();
  } else {
    nrThreads = obs.getNrPaddedDMs() / snrDParameters[deviceName][obs.getNrDMs()].getNrDMsPerThread();
  }
  cl::NDRange snrDedispersedGlobal(nrThreads);
  cl::NDRange snrDedispersedLocal(snrDParameters[deviceName][obs.getNrDMs()].getNrDMsPerBlock());
  if ( DEBUG == 0 ) {
    std::cout << "SNRDedispersed" << std::endl;
    std::cout << "Global: " << nrThreads << std::endl;
    std::cout << "Local: " << snrDParameters[deviceName][obs.getNrDMs()].getNrDMsPerBlock() << std::endl;
    std::cout << "Parameters: ";
    std::cout << snrDParameters[deviceName][obs.getNrDMs()].print() << std::endl;
    std::cout << std::endl;
  }
  if ( obs.getNrDMs() % (snrFParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrDMsPerBlock() * snrFParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrDMsPerThread()) == 0 ) {
    nrThreads = obs.getNrDMs() / snrFParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrDMsPerThread();
  } else {
    nrThreads = obs.getNrPaddedDMs() / snrFParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrDMsPerThread();
  }
  cl::NDRange snrFoldedGlobal(nrThreads, obs.getNrPeriods() / snrFParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrPeriodsPerThread());
  cl::NDRange snrFoldedLocal(snrFParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrDMsPerBlock(), snrFParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrPeriodsPerBlock());
  if ( DEBUG == 0 ) {
    std::cout << "SNRFolded" << std::endl;
    std::cout << "Global: " << nrThreads << ", " << obs.getNrPeriods() / snrFParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrPeriodsPerThread() << std::endl;
    std::cout << "Local: " << snrFParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrDMsPerBlock() << ", " << snrFParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].getNrPeriodsPerBlock() << std::endl;
    std::cout << "Parameters: ";
    std::cout << snrFParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()].print() << std::endl;
    std::cout << std::endl;
  }

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
  isa::utils::Timer inputHandlingTime;
  isa::utils::Timer inputCopyTime;
  isa::utils::Timer dedispTime;
  isa::utils::Timer transTime;
  isa::utils::Timer foldTime;
  isa::utils::Timer snrDedispersedTime;
  isa::utils::Timer snrFoldedTime;
  isa::utils::Timer outputCopyTime;
  isa::utils::Timer foldedTSCopyTime;

  searchTime.start();
	for ( unsigned int second = 0; second < obs.getNrSeconds() - secondsToBuffer; second++ ) {
		// Load the input
    inputHandlingTime.start();
		for ( unsigned int channel = 0; channel < obs.getNrChannels(); channel++ ) {
			for ( unsigned int chunk = 0; chunk < secondsToBuffer; chunk++ ) {
        memcpy(reinterpret_cast< void * >(&(dispersedData.data()[(channel * obs.getNrSamplesPerDispersedChannel()) + (chunk * obs.getNrSamplesPerSecond())])), reinterpret_cast< void * >(&((input->at(second + chunk))->at(channel * obs.getNrSamplesPerPaddedSecond()))), obs.getNrSamplesPerSecond() * sizeof(dataType));
			}
      memcpy(reinterpret_cast< void * >(&(dispersedData.data()[(channel * obs.getNrSamplesPerDispersedChannel()) + (secondsToBuffer * obs.getNrSamplesPerSecond())])), reinterpret_cast< void * >(&((input->at(second + secondsToBuffer))->at(channel * obs.getNrSamplesPerPaddedSecond()))), remainingSamples * sizeof(dataType));
		}
    try {
      inputCopyTime.start();
      clQueues->at(clDeviceID)[0].enqueueWriteBuffer(dispersedData_d, CL_TRUE, 0, dispersedData.size() * sizeof(dataType), reinterpret_cast< void * >(dispersedData.data()), 0, &syncPoint);
      syncPoint.wait();
      inputCopyTime.stop();
      if ( DEBUG ) {
        if ( print == 0 ) {
          std::cout << std::fixed << std::setprecision(3);
          for ( unsigned int channel = 0; channel < obs.getNrChannels(); channel++ ) {
            std::cout << channel << " : ";
            for ( unsigned int sample = 0; sample < obs.getNrSamplesPerDispersedChannel(); sample++ ) {
              std::cout << dispersedData[(channel * obs.getNrSamplesPerDispersedChannel()) + sample] << " ";
            }
            std::cout << std::endl;
          }
          std::cout << std::endl;
        }
      }
    } catch ( cl::Error & err ) {
      std::cerr << err.what() << std::endl;
      return 1;
    }
    inputHandlingTime.stop();

		// Run the kernels
		try {
      dedispTime.start();
      clQueues->at(clDeviceID)[0].enqueueNDRangeKernel(*dedispersionK, cl::NullRange, dedispersionGlobal, dedispersionLocal, 0, &syncPoint);
      syncPoint.wait();
      dedispTime.stop();
      if ( DEBUG ) {
        if ( print == 0 ) {
          clQueues->at(clDeviceID)[0].enqueueReadBuffer(dedispersedData_d, CL_TRUE, 0, dedispersedData.size() * sizeof(dataType), reinterpret_cast< void * >(dedispersedData.data()));
          std::cout << std::fixed << std::setprecision(3);
          for ( unsigned int dm = 0; dm < obs.getNrDMs(); dm++ ) {
            std::cout << dm << " : ";
            for ( unsigned int sample = 0; sample < obs.getNrSamplesPerSecond(); sample++ ) {
              std::cout << dedispersedData[(dm * obs.getNrSamplesPerPaddedSecond()) + sample] << " ";
            }
            std::cout << std::endl;
          }
          std::cout << std::endl;
        }
      }
      transTime.start();
      clQueues->at(clDeviceID)[0].enqueueNDRangeKernel(*transposeK, cl::NullRange, transposeGlobal, transposeLocal, 0, &syncPoint);
      syncPoint.wait();
      transTime.stop();
      if ( DEBUG ) {
        if ( print == 0 ) {
          clQueues->at(clDeviceID)[0].enqueueReadBuffer(transposedData_d, CL_TRUE, 0, transposedData.size() * sizeof(dataType), reinterpret_cast< void * >(transposedData.data()));
          std::cout << std::fixed << std::setprecision(3);
          for ( unsigned int sample = 0; sample < obs.getNrSamplesPerSecond(); sample++ ) {
            std::cout << sample << " : ";
            for ( unsigned int dm = 0; dm < obs.getNrDMs(); dm++ ) {
              std::cout << transposedData[(sample * obs.getNrPaddedDMs()) + dm] << " ";
            }
            std::cout << std::endl;
          }
          std::cout << std::endl;
        }
      }
      snrDedispersedK->setArg(0, static_cast< float >(second));
      snrDedispersedTime.start();
      clQueues->at(clDeviceID)[0].enqueueNDRangeKernel(*snrDedispersedK, cl::NullRange, snrDedispersedGlobal, snrDedispersedLocal, 0, &syncPoint);
      syncPoint.wait();
      snrDedispersedTime.stop();
      if ( DEBUG ) {
        if ( print == 0 ) {
          clQueues->at(clDeviceID)[0].enqueueReadBuffer(maxDedispersedTable_d, CL_TRUE, 0, maxDedispersedTable.size() * sizeof(dataType), reinterpret_cast< void * >(maxDedispersedTable.data()));
          clQueues->at(clDeviceID)[0].enqueueReadBuffer(meanDedispersedTable_d, CL_TRUE, 0, meanDedispersedTable.size() * sizeof(float), reinterpret_cast< void * >(meanDedispersedTable.data()));
          clQueues->at(clDeviceID)[0].enqueueReadBuffer(rmsDedispersedTable_d, CL_TRUE, 0, rmsDedispersedTable.size() * sizeof(float), reinterpret_cast< void * >(rmsDedispersedTable.data()));
          std::cout << std::fixed << std::setprecision(6);
          for ( unsigned int dm = 0; dm < obs.getNrDMs(); dm++ ) {
            std::cout << dm << ": " << maxDedispersedTable[dm] << " " << meanDedispersedTable[dm] << " " << std::sqrt(rmsDedispersedTable[dm]) << " " << (maxDedispersedTable[dm] - meanDedispersedTable[dm]) / std::sqrt(rmsDedispersedTable[dm]) << std::endl;
          }
          std::cout << std::endl;
        }
      }
      foldingK->setArg(0, second);
			if ( second % 2 == 0 ) {
        foldingK->setArg(3, counterData0_d);
        foldingK->setArg(4, counterData1_d);
			} else {
        foldingK->setArg(3, counterData1_d);
        foldingK->setArg(4, counterData0_d);
			}
      foldTime.start();
      clQueues->at(clDeviceID)[0].enqueueNDRangeKernel(*foldingK, cl::NullRange, foldingGlobal, foldingLocal, 0, &syncPoint);
      syncPoint.wait();
      foldTime.stop();
      if ( DEBUG ) {
        if ( print == 0 ) {
          clQueues->at(clDeviceID)[0].enqueueReadBuffer(foldedData_d, CL_TRUE, 0, foldedData.size() * sizeof(dataType), reinterpret_cast< void * >(foldedData.data()));
          std::cout << std::fixed << std::setprecision(3);
          for ( unsigned int dm = 0; dm < obs.getNrDMs(); dm++ ) {
            for ( unsigned int period = 0; period < obs.getNrPeriods(); period++ ) {
              std::cout << dm << " " << period << " : ";
              for ( unsigned int bin = 0; bin < obs.getNrBins(); bin++ ) {
                std::cout << foldedData[(bin * obs.getNrPeriods() * obs.getNrPaddedDMs()) + (period * obs.getNrPaddedDMs()) + dm] << " ";
              }
              std::cout << std::endl;
            }
            std::cout << std::endl;
          }
          std::cout << std::endl;
        }
      }
      if ( second == (obs.getNrSeconds() - secondsToBuffer) - 1 ) {
        snrFoldedTime.start();
        clQueues->at(clDeviceID)[0].enqueueNDRangeKernel(*snrFoldedK, cl::NullRange, snrFoldedGlobal, snrFoldedLocal, 0, &syncPoint);
        syncPoint.wait();
        snrFoldedTime.stop();
      }
		} catch ( cl::Error & err ) {
			std::cerr << err.what() <<" "  << err.err() << std::endl;
			return 1;
		}
	}
  searchTime.stop();

  // Copy output from device to host
	try {
    outputCopyTime.start();
    clQueues->at(clDeviceID)[0].enqueueReadBuffer(snrFoldedTable_d, CL_TRUE, 0, snrFoldedTable.size() * sizeof(float), reinterpret_cast< void * >(snrFoldedTable.data()), 0, &syncPoint);
    syncPoint.wait();
    clQueues->at(clDeviceID)[0].enqueueReadBuffer(maxDedispersedTable_d, CL_TRUE, 0, maxDedispersedTable.size() * sizeof(dataType), reinterpret_cast< void * >(maxDedispersedTable.data()), 0, &syncPoint);
    syncPoint.wait();
    clQueues->at(clDeviceID)[0].enqueueReadBuffer(meanDedispersedTable_d, CL_TRUE, 0, meanDedispersedTable.size() * sizeof(float), reinterpret_cast< void * >(meanDedispersedTable.data()), 0, &syncPoint);
    syncPoint.wait();
    clQueues->at(clDeviceID)[0].enqueueReadBuffer(rmsDedispersedTable_d, CL_TRUE, 0, rmsDedispersedTable.size() * sizeof(float), reinterpret_cast< void * >(rmsDedispersedTable.data()), 0, &syncPoint);
    syncPoint.wait();
    outputCopyTime.stop();
    foldedTSCopyTime.start();
    clQueues->at(clDeviceID)[0].enqueueReadBuffer(foldedData_d, CL_TRUE, 0, foldedData.size() * sizeof(dataType), reinterpret_cast< void * >(foldedData.data()), 0, &syncPoint);
    syncPoint.wait();
    foldedTSCopyTime.stop();
  } catch ( cl::Error & err ) {
		std::cerr << err.what() << std::endl;
		return 1;
	}

  output.sync_with_stdio(false);
	// Store output
  if ( saveOutput ) {
    output.open(outputFile + ".fold");
    output << std::fixed << std::setprecision(6);
    output << "# bin SNR" << std::endl;
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
    output.open(outputFile + ".foldSNR");
    output << "# period DM SNR" << std::endl;
    output << std::fixed << std::setprecision(6);
    for ( unsigned int period = 0; period < obs.getNrPeriods(); period++ ) {
      for ( unsigned int dm = 0; dm < obs.getNrDMs(); dm++ ) {
        output << period << " ";
        output << dm << " ";
        output << snrFoldedTable[(period * obs.getNrPaddedDMs()) + dm] << std::endl;
      }
    }
    output.close();
    output.open(outputFile + ".dediSNR");
    output << "# DM SNR" << std::endl;
    for ( unsigned int dm = 0; dm < obs.getNrDMs(); dm++ ) {
      output << dm << " " << (maxDedispersedTable[dm] - meanDedispersedTable[dm]) / std::sqrt(rmsDedispersedTable[dm]) << std::endl;
    }
    output.close();
  }
  // Store statistics
	output.open(outputFile + ".stats");
  output << "# nrDMs nrPeriods searchTime inputHandlingTotal inputHandlingAvg err inputCopyTotal inputCopyAvg err dedispersionTotal dedispersionvg err transposeTotal transposeAvg err snrDedispersedTotal snrDedispersedAvg err foldingTotal foldingAvg err snrFoldedTotal snrFoldedAvg err outputCopyTotal outputCopyAvg err foldedTSCopyTotal foldedTSCopyAvg err" << std::endl;
  output << std::fixed << std::setprecision(6);
  output << obs.getNrDMs() << " " << obs.getNrPeriods() << " ";
  output << searchTime.getTotalTime() << " ";
  output << inputHandlingTime.getTotalTime() << " " << inputHandlingTime.getAverageTime() << " " << inputHandlingTime.getStandardDeviation() << " ";
  output << inputCopyTime.getTotalTime() << " " << inputCopyTime.getAverageTime() << " " << inputCopyTime.getStandardDeviation() << " ";
  output << dedispTime.getTotalTime() << " " << dedispTime.getAverageTime() << " " << dedispTime.getStandardDeviation() << " ";
  output << transTime.getTotalTime() << " " << transTime.getAverageTime() << " " << transTime.getStandardDeviation() << " ";
  output << snrDedispersedTime.getTotalTime() << " " << snrDedispersedTime.getAverageTime() << " " << snrDedispersedTime.getStandardDeviation() << " ";
  output << foldTime.getTotalTime() << " " << foldTime.getAverageTime() << " " << foldTime.getStandardDeviation() << " ";
  output << snrFoldedTime.getTotalTime() << " " << snrFoldedTime.getAverageTime() << " " << snrFoldedTime.getStandardDeviation() << " ";
  output << outputCopyTime.getTotalTime() << " " << outputCopyTime.getAverageTime() << " " << outputCopyTime.getStandardDeviation() << " ";
  output << foldedTSCopyTime.getTotalTime() << " " << foldedTSCopyTime.getAverageTime() << " " << foldedTSCopyTime.getStandardDeviation() << " ";
  output << std::endl;
  output.close();

	return 0;
}

