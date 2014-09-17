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
	std::string deviceName;
	std::string dataFile;
	std::string headerFile;
	std::string outputFile;
	// Observation object
  AstroData::Observation< dataType > obs("PulsarSearch", dataName);

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

		obs.setNrDMs(args.getSwitchArgument< unsigned int >("-dm_node"));
		obs.setDMStep(args.getSwitchArgument< float >("-dm_step"));
		obs.setFirstDM(args.getSwitchArgument< float >("-dm_first") + ((world.rank() / MPIRows) * obs.getNrDMs() * obs.getDMStep()));
		obs.setNrPeriods(args.getSwitchArgument< unsigned int >("-period_node"));
		obs.setPeriodStep(args.getSwitchArgument< unsigned int >("-period_step"));
		obs.setFirstPeriod(args.getSwitchArgument< unsigned int >("-period_first") + ((world.rank() % MPICols) * obs.getPeriodStep() * obs.getNrPeriods()));
		obs.setNrBins(args.getSwitchArgument< unsigned int >("-period_bins"));
	} catch ( std::exception & err ) {
		std::cerr << err.what() << std::endl;
		return 1;
	}

	// Load observation data
  isa::utils::Timer loadTime("LoadInputTimer");
	std::vector< CLData< dataType > * > * input = new std::vector< CLData< dataType > * >(1);
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
	} catch ( OpenCLError & err ) {
		std::cerr << err.what() << std::endl;
		return 1;
	}

	// Host memory allocation
  std::vector< unsigned int > * shifts = PulsarSearch::getShifts(obs);
  obs.setNrSamplesPerDispersedChannel(obs.getNrSamplesPerSecond() + (*shifts)[((obs.getNrDMs() - 1) * obs.getNrPaddedChannels())]);
  std::vector< unsigned int > * nrSamplesPerBin = PulsarSearch::getSamplesPerBin(obs);
  std::vector< dataType > dispersedData(obs.getNrChannels() * obs.getNrSamplesPerDispersedChannel());
  std::vector< dataType > snrTable(obs.getNrPeriods() * obs.getNrPaddedDMs());

  // Device memory allocation and data transfers
  cl::Buffer shifts_d, nrSamplesPerBin_d, dispersedData_d, dedispersedData_d, transposedData_d, foldedData_d, counterData0_d, counterData1_d, snrTable_d;

  try {
    shifts_d = cl::Buffer(*clContext, CL_MEM_READ_ONLY, shifts->size() * sizeof(unsigned int), 0, 0);
    delete shifts;
    samplesPerBin_d = cl::Buffer(*clContext, CL_MEM_READ_ONLY, samplesPerBin->size() * sizeof(unsigned int), 0, 0);
    delete samplesPerBin;
    dispersedData_d = cl::Buffer(*clContext, CL_MEM_READ_ONLY, dispersedData.size() * sizeof(dataType), 0, 0);
    dedispersedData_d = cl::Buffer(*clContext, CL_MEM_READ_WRITE, obs.getNrDMs() * obs.getNrSamplesPerPaddedSecond() * sizeof(dataType), 0, 0);
    transposedData_d = cl::Buffer(*clContext, CL_MEM_READ_WRITE, obs.getNrSamplesPerSecond() * obs.getNrPaddedDMs() * sizeof(dataType), 0, 0);
  } catch ( isa::OpenCL::OpenCLError & err ) {
    std::cerr << err.what() << std::endl;
    return 1;
  }

	try {
		// Shifts
		shifts->setCLContext(clContext);
		shifts->setCLQueue(&((clQueues->at(clDeviceID)).at(0)));
		shifts->allocateDeviceData();
		shifts->copyHostToDevice();
		// nrSamplesPerBin
		nrSamplesPerBin.allocateHostData(*(getNrSamplesPerBin(obs)));
		nrSamplesPerBin.setCLContext(clContext);
		nrSamplesPerBin.setCLQueue(&((clQueues->at(clDeviceID)).at(0)));
		nrSamplesPerBin.setDeviceReadOnly();
		nrSamplesPerBin.allocateDeviceData();
		nrSamplesPerBin.copyHostToDevice();
		nrSamplesPerBin.deleteHostData();
		// DispersedData
		dispersedData.allocateHostData(secondsToBuffer * obs.getNrChannels() * obs.getNrSamplesPerPaddedSecond());
		dispersedData.setCLContext(clContext);
		dispersedData.setCLQueue(&((clQueues->at(clDeviceID)).at(0)));
		dispersedData.setDeviceReadOnly();
		dispersedData.allocateDeviceData();
		// DedispersedData
		dedispersedData.setCLContext(clContext);
		dedispersedData.setCLQueue(&((clQueues->at(clDeviceID)).at(0)));
		dedispersedData.allocateDeviceData(obs.getNrDMs() * obs.getNrSamplesPerPaddedSecond());
		// TransposedData
		transposedData.setCLContext(clContext);
		transposedData.setCLQueue(&((clQueues->at(clDeviceID)).at(0)));
		transposedData.allocateDeviceData(obs.getNrSamplesPerSecond() * obs.getNrPaddedDMs());
		// FoldedData
		foldedData.setCLContext(clContext);
		foldedData.setCLQueue(&((clQueues->at(clDeviceID)).at(0)));
		foldedData.allocateHostData(obs.getNrBins() * obs.getNrPeriods() * obs.getNrPaddedDMs());
		foldedData.allocateDeviceData();
		foldedData.blankHostData();
		foldedData.copyHostToDevice();
		foldedData.deleteHostData();
		// CounterData0
		counterData0.setCLContext(clContext);
		counterData0.setCLQueue(&((clQueues->at(clDeviceID)).at(0)));
		counterData0.allocateHostData(obs.getNrPeriods() * obs.getNrPaddedBins());
		counterData0.blankHostData();
		counterData0.allocateDeviceData();
		counterData0.copyHostToDevice();
		counterData0.deleteHostData();
		// CounterData1
		counterData1.setCLContext(clContext);
		counterData1.setCLQueue(&((clQueues->at(clDeviceID)).at(0)));
		counterData1.allocateHostData(obs.getNrPeriods() * obs.getNrPaddedBins());
		counterData1.blankHostData();
		counterData1.allocateDeviceData();
		counterData1.copyHostToDevice();
		counterData1.deleteHostData();
		// SNRData
		snrTable.allocateHostData(obs.getNrPeriods() * obs.getNrPaddedDMs());
		snrTable.setCLContext(clContext);
		snrTable.setCLQueue(&((clQueues->at(clDeviceID)).at(0)));
		snrTable.setDeviceWriteOnly();
		snrTable.allocateDeviceData();
	} catch ( OpenCLError &err ) {
		std::cerr << err.what() << std::endl;
		return 1;
	}

	if ( DEBUG && world.rank() == 0 ) {
		double hostMemory = 0.0;
		double deviceMemory = 0.0;
		hostMemory += dispersedData.getHostDataSize() + snrTable.getHostDataSize();
		deviceMemory += shifts->getDeviceDataSize() + nrSamplesPerBin.getDeviceDataSize() + dispersedData.getDeviceDataSize() + dedispersedData.getDeviceDataSize() + transposedData.getDeviceDataSize() + foldedData.getDeviceDataSize() + (2 * counterData0.getDeviceDataSize()) + snrTable.getDeviceDataSize();

		std::cout << "Allocated host memory: " << std::fixed << std::setprecision(3) << isa::utils::giga(hostMemory) << " GB." << std::endl;
		std::cout << "Allocated device memory: " << std::fixed << std::setprecision(3) << isa::utils::giga(deviceMemory) << "GB." << std::endl;
	}

	// Generate OpenCL kernels
	Dedispersion< dataType > clDedisperse("clDedisperse", dataName);
	Transpose< dataType > clTranspose("clTranspose", dataName);
	Folding< dataType > clFold("clFold", dataName);
	SNR< dataType > clSNR("clSNR", dataName);
	try {
		// Dedispersion
		clDedisperse.bindOpenCL(clContext, &(clDevices->at(clDeviceID)), &((clQueues->at(clDeviceID)).at(0)));
		clDedisperse.setNrSamplesPerBlock(dedispersionParameters[deviceName][obs.getNrDMs()][0]);
		clDedisperse.setNrDMsPerBlock(dedispersionParameters[deviceName][obs.getNrDMs()][1]);
		clDedisperse.setNrSamplesPerThread(dedispersionParameters[deviceName][obs.getNrDMs()][2]);
		clDedisperse.setNrDMsPerThread(dedispersionParameters[deviceName][obs.getNrDMs()][3]);
		clDedisperse.setNrSamplesPerDispersedChannel(secondsToBuffer * obs.getNrSamplesPerPaddedSecond());
		clDedisperse.setObservation(&obs);
		clDedisperse.setShifts(shifts);
		clDedisperse.generateCode();
		shifts->deleteHostData();
		// Transposition
		clTranspose.bindOpenCL(clContext, &(clDevices->at(clDeviceID)), &((clQueues->at(clDeviceID)).at(0)));
		clTranspose.setNrThreadsPerBlock(transposeParameters[deviceName][obs.getNrDMs()]);
		clTranspose.setDimensions(obs.getNrDMs(), obs.getNrSamplesPerSecond());
		clTranspose.setPaddingFactor(padding[deviceName]);
		clTranspose.setVectorWidth(std::vectorWidth[deviceName]);
		clTranspose.generateCode();
		// Folding
		clFold.bindOpenCL(clContext, &(clDevices->at(clDeviceID)), &((clQueues->at(clDeviceID)).at(0)));
		clFold.setObservation(&obs);
		clFold.setNrSamplesPerBin(&nrSamplesPerBin);
		clFold.setNrDMsPerBlock(foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0]);
		clFold.setNrPeriodsPerBlock(foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][1]);
		clFold.setNrBinsPerBlock(foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][2]);
		clFold.setNrDMsPerThread(foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][3]);
		clFold.setNrPeriodsPerThread(foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][4]);
		clFold.setNrBinsPerThread(foldingParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][5]);
		clFold.generateCode();
		// SNR
		clSNR.bindOpenCL(clContext, &(clDevices->at(clDeviceID)), &((clQueues->at(clDeviceID)).at(0)));
		clSNR.setObservation(&obs);
		clSNR.setNrDMsPerBlock(snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][0]);
		clSNR.setNrPeriodsPerBlock(snrParameters[deviceName][obs.getNrDMs()][obs.getNrPeriods()][1]);
		clSNR.setPulsarPipeline();
		clSNR.generateCode();
	} catch ( OpenCLError &err ) {
		std::cerr << err.what() << std::endl;
		return 1;
	}

	// Search loop
  isa::utils::Timer searchTime("SearchTimer");
  isa::utils::Timer inputPreTime("InputPreProcessingTimer");
	if ( DEBUG && world.rank() == 0 ) {
		std::cout << "Starting the search." << std::endl;
		std::cout << "Processing seconds: ";
	}
	for ( unsigned int second = 0; second <= obs.getNrSeconds() - secondsToBuffer; second++ ) {
		if ( DEBUG && world.rank() == 0 ) {
			std::cout << second << " " << std::flush;
		}
		searchTime.start();
		// Prepare the input
		inputPreTime.start();
		for ( unsigned int channel = 0; channel < obs.getNrChannels(); channel++ ) {
			for ( unsigned int chunk = 0; chunk < secondsToBuffer; chunk++ ) {
				memcpy(dispersedData.getRawHostDataAt((channel * secondsToBuffer * obs.getNrSamplesPerPaddedSecond()) + (chunk * obs.getNrSamplesPerSecond())), (input->at(second + chunk))->getRawHostDataAt(channel * obs.getNrSamplesPerPaddedSecond()), obs.getNrSamplesPerSecond() * sizeof(dataType));
			}
		}
		inputPreTime.stop();

		// Run the kernels
		try {
			dispersedData.copyHostToDevice();
			clDedisperse(&dispersedData, &dedispersedData);
			clTranspose(&dedispersedData, &transposedData);
			if ( second % 2 == 0 ) {
				clFold(second, &transposedData, &foldedData, &counterData0, &counterData1);
			} else {
				clFold(second, &transposedData, &foldedData, &counterData1, &counterData0);
			}
		} catch ( OpenCLError &err ) {
			std::cerr << err.what() << std::endl;
			return 1;
		}
		searchTime.stop();
	}
	if ( DEBUG && world.rank() == 0 ) {
		std::cout << "." << std::endl;
		std::cout << "Search complete." << std::endl;
	}

	// Release unnecessary memory
	delete input;
	try {
		dispersedData.deleteDeviceData();
		dispersedData.deleteHostData();
		dedispersedData.deleteDeviceData();
		transposedData.deleteDeviceData();
		counterData0.deleteDeviceData();
		counterData1.deleteDeviceData();
	} catch ( OpenCLError &err ) {
		std::cerr << err.what() << std::endl;
	}

	// Store output
  isa::utils::Timer outputTime("OutputTimer");
	if ( DEBUG && world.rank() == 0 ) {
		std::cout << "Analyzing processed data." << std::endl;
	}
	try {
		outputTime.start();
		clSNR(&foldedData, &snrTable);
		snrTable.copyDeviceToHost();
		outputTime.stop();
		foldedData.deleteDeviceData();
		foldedData.deleteHostData();
	} catch ( OpenCLError &err ) {
		std::cerr << err.what() << std::endl;
		return 1;
	}
	if ( DEBUG && world.rank() == 0 ) {
		std::cout << "Saving output to disk." << std::endl;
	}
	std::ofstream output;
	output.open(outputFile + "_" + isa::utils::toString< unsigned int >(world.rank()));
	for ( unsigned int period = 0; period < obs.getNrPeriods(); period++ ) {
		for ( unsigned int dm = 0; dm < obs.getNrDMs(); dm++ ) {
			output << std::fixed << std::setprecision(6) << (world.rank() * obs.getNrPeriods()) + period << " " << (obs.getFirstPeriod() + (period * obs.getPeriodStep())) / static_cast< float > (obs.getNrSamplesPerSecond()) << " " << dm << " " << obs.getFirstDM() + (dm * obs.getDMStep()) << " " << snrTable[(period * obs.getNrPaddedDMs()) + dm] << std::endl;
		}
	}
	output.close();

	try {
		snrTable.deleteDeviceData();
		snrTable.deleteHostData();
	} catch ( OpenCLError &err ) {
		std::cerr << err.what() << std::endl;
	}

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
		std::cout << "# nodes accelerator processedSeconds nrDMs nrPeriods nrBins nrSamplesPerSecond searchTime mainLoopTime mainLoopAverageTime searchGFLOPs mainLoopGFLOPs" << std::endl;
		std::cout << std::fixed << world.size() << " " << deviceName << " " << 1 + obs.getNrSeconds() - secondsToBuffer << " " << obs.getNrDMs() << " " << obs.getNrPeriods() << " " << obs.getNrBins() << " " << obs.getNrSamplesPerSecond() << " " << std::setprecision(6) << maxTime << " " << maxKernel << " " << maxKernel / (1 + obs.getNrSeconds() - secondsToBuffer) << " " << std::setprecision(3) << (world.size() * (((1 + obs.getNrSeconds() - secondsToBuffer) * (clDedisperse.getGFLOP() + clFold.getGFLOP())) + clSNR.getGFLOP())) / maxTime << " " << (world.size() * (((1 + obs.getNrSeconds() - secondsToBuffer) * (clDedisperse.getGFLOP() + clFold.getGFLOP())) + clSNR.getGFLOP())) / maxKernel  << std::endl;
	}

	return 0;
}

