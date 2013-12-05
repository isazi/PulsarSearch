//
// Copyright (C) 2013
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
using std::cout;
using std::cerr;
using std::endl;
using std::fixed;
#include <string>
using std::string;
#include <exception>
using std::exception;
#include <fstream>
using std::ofstream;
#include <vector>
using std::vector;
#include <iomanip>
using std::setprecision;
#include <boost/mpi.hpp>
using boost::mpi::environment;
using boost::mpi::communicator;
#include <map>
using std::map;

#include <configuration.hpp>
#include <readConfiguration.hpp>

#include <ArgumentList.hpp>
using isa::utils::ArgumentList;
#include <utils.hpp>
using isa::utils::giga;
using isa::utils::toStringValue;
#include <Observation.hpp>
using AstroData::Observation;
#include <ReadData.hpp>
using AstroData::readSIGPROC;
using AstroData::readLOFAR;
#include <Dedispersion.hpp>
using PulsarSearch::Dedispersion;
#include <Shifts.hpp>
using PulsarSearch::getShifts;
#include <Bins.hpp>
using PulsarSearch::getNrSamplesPerBin;
#include <Folding.hpp>
using PulsarSearch::Folding;
#include <Transpose.hpp>
using isa::OpenCL::Transpose;
#include <SNR.hpp>
using PulsarSearch::SNR;
#include <CLData.hpp>
using isa::OpenCL::CLData;
#include <Exceptions.hpp>
using isa::Exceptions::OpenCLError;
#include <InitializeOpenCL.hpp>
using isa::OpenCL::initializeOpenCL;
#include <Timer.hpp>
using isa::utils::Timer;


int main(int argc, char * argv[]) {
	bool dataLOFAR = false;
	bool dataSIGPROC = false;
	unsigned int clPlatformID = 0;
	unsigned int clDeviceID = 0;
	string deviceName;
	string dataFile;
	string headerFile;
	string outputFile;
	// Observation object
	Observation< dataType > obs("PulsarSearch", dataName);

	// Initialize MPI
	environment envMPI(argc, argv);
	communicator world;

	try {
		ArgumentList args(argc, argv);

		clPlatformID = args.getSwitchArgument< unsigned int >("-opencl_platform");
		clDeviceID = args.getSwitchArgument< unsigned int >("-opencl_device");
		deviceName = args.getSwitchArgument< string >("-device_name");

		readPadding(padding, args.getSwitchArgument< string >("-padding_file"));
		readVectorWidth(vectorWidth, args.getSwitchArgument< string >("-vector_file"));
		readDedispersion(dedispersionParameters, args.getSwitchArgument< string >("-dedispersion_file"));
		readTranspose(transposeParameters, args.getSwitchArgument< string >("-transpose_file"));
		readFolding(foldingParameters, args.getSwitchArgument< string >("-folding_file"));
		readSNR(snrParameters, args.getSwitchArgument< string >("-snr_file"));

		obs.setPadding(padding[deviceName]);

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

		obs.setNrDMs(args.getSwitchArgument< unsigned int >("-dm_number"));
		obs.setFirstDM(args.getSwitchArgument< float >("-dm_first"));
		obs.setDMStep(args.getSwitchArgument< float >("-dm_step"));
		obs.setNrPeriods(args.getSwitchArgument< unsigned int >("-period_number"));
		obs.setFirstPeriod(args.getSwitchArgument< unsigned int >("-period_first") + (world.rank() * obs.getNrPeriods()));
		obs.setPeriodStep(args.getSwitchArgument< unsigned int >("-period_step"));
		obs.setNrBins(args.getSwitchArgument< unsigned int >("-period_bins"));
	} catch ( exception &err ) {
		cerr << err.what() << endl;
		return 1;
	}

	// Load observation data
	Timer loadTime("LoadInputTimer");
	vector< CLData< dataType > * > * input = new vector< CLData< dataType > * >(1);
	loadTime.start();
	if ( dataLOFAR ) {
		readLOFAR(headerFile, dataFile, obs, *input);
	} else if ( dataSIGPROC ) {
		// TODO: implement SIGPROC pipeline
	}
	loadTime.stop();
	if ( DEBUG && world.rank() == 0 ) {
		cout << "Time to load the input: " << fixed << setprecision(6) << loadTime.getTotalTime() << " seconds." << endl;
	}

	// Initialize OpenCL
	cl::Context *clContext = new cl::Context();
	vector< cl::Platform > *clPlatforms = new vector< cl::Platform >();
	vector< cl::Device > *clDevices = new vector< cl::Device >();
	vector< vector< cl::CommandQueue > > *clQueues = new vector< vector < cl::CommandQueue > >();
	
	try {
		initializeOpenCL(clPlatformID, 1, clPlatforms, clContext, clDevices, clQueues);
	} catch ( OpenCLError err ) {
		cerr << err.what() << endl;
		return 1;
	}

	// Memory allocation
	unsigned int nrSamplesPerChannel = 0;
	unsigned int secondsToBuffer = 0;
	CLData< unsigned int > * shifts = getShifts(obs);
	CLData< unsigned int > nrSamplesPerBin("nrSamplesPerBin", true);
	CLData< dataType > dispersedData("DispersedData", true);
	CLData< dataType > dedispersedData("DedispersedData", true);
	CLData< dataType > transposedData("TransposedData", true);
	CLData< dataType > foldedData("FoldedData", true);
	CLData< unsigned int > counterData0("CounterData0", true);
	CLData< unsigned int > counterData1("CounterData1", true);
	CLData< dataType > snrTable("SNRData", true);

	if ( ((obs.getNrSamplesPerSecond() + (*shifts)[((obs.getNrDMs() - 1) * obs.getNrPaddedChannels())]) % obs.getPadding()) != 0 ) {
		nrSamplesPerChannel = (obs.getNrSamplesPerSecond() + (*shifts)[((obs.getNrDMs() - 1) * obs.getNrPaddedChannels())]) + (obs.getPadding() - ((obs.getNrSamplesPerSecond() + (*shifts)[((obs.getNrDMs() - 1) * obs.getNrPaddedChannels())]) % obs.getPadding()));
	} else {
		nrSamplesPerChannel = (obs.getNrSamplesPerSecond() + (*shifts)[((obs.getNrDMs() - 1) * obs.getNrPaddedChannels())]);
	}
	secondsToBuffer = static_cast< unsigned int >(ceil(static_cast< float >(nrSamplesPerChannel) / obs.getNrSamplesPerPaddedSecond()));

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
		nrSamplesPerBin.allocateDeviceData();
		nrSamplesPerBin.copyHostToDevice();
		nrSamplesPerBin.deleteHostData();
		// DispersedData
		dispersedData.allocateHostData(secondsToBuffer * obs.getNrChannels() * obs.getNrSamplesPerPaddedSecond());
		dispersedData.setCLContext(clContext);
		dispersedData.setCLQueue(&((clQueues->at(clDeviceID)).at(0)));
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
		snrTable.allocateDeviceData();
	} catch ( OpenCLError err ) {
		cerr << err.what() << endl;
		return 1;
	}

	if ( DEBUG && world.rank() == 0 ) {
		double hostMemory = 0.0;
		double deviceMemory = 0.0;
		hostMemory += dispersedData.getHostDataSize() + snrTable.getHostDataSize();
		deviceMemory += shifts->getDeviceDataSize() + nrSamplesPerBin.getDeviceDataSize() + dispersedData.getDeviceDataSize() + dedispersedData.getDeviceDataSize() + transposedData.getDeviceDataSize() + foldedData.getDeviceDataSize() + (2 * counterData0.getDeviceDataSize()) + snrTable.getDeviceDataSize();

		cout << "Allocated host memory: " << fixed << setprecision(3) << giga(hostMemory) << " GB." << endl;
		cout << "Allocated device memory: " << fixed << setprecision(3) << giga(deviceMemory) << "GB." << endl;
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
		clTranspose.setVectorWidth(vectorWidth[deviceName]);
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
	} catch ( OpenCLError err ) {
		cerr << err.what() << endl;
		return 1;
	}

	// Search loop
	Timer searchTime("SearchTimer");
	Timer inputPreTime("InputPreProcessingTimer");
	if ( DEBUG && world.rank() == 0 ) {
		cout << "Starting the search." << endl;
		cout << "Processing seconds: ";
	}
	for ( unsigned int second = 0; second <= obs.getNrSeconds() - secondsToBuffer; second++ ) {
		if ( DEBUG && world.rank() == 0 ) {
			cout << second << " ";
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
		} catch ( OpenCLError err ) {
			cerr << err.what() << endl;
			return 1;
		}
		searchTime.stop();
	}
	if ( DEBUG && world.rank() == 0 ) {
		cout << "." << endl;
		cout << "Search complete." << endl;
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
	} catch ( OpenCLError err ) {
		cerr << err.what() << endl;
	}

	// Store output
	Timer outputTime("OutputTimer");
	if ( DEBUG && world.rank() == 0 ) {
		cout << "Analyzing processed data." << endl;
	}
	try {
		outputTime.start();
		clSNR(&foldedData, &snrTable);
		snrTable.copyDeviceToHost();
		outputTime.stop();
		foldedData.deleteDeviceData();
		foldedData.deleteHostData();
	} catch ( OpenCLError err ) {
		cerr << err.what() << endl;
		return 1;
	}
	if ( DEBUG && world.rank() == 0 ) {
		cout << "Saving output to disk." << endl;
	}
	ofstream output;
	output.open(outputFile + "_" + toStringValue< unsigned int >(world.rank()));
	for ( unsigned int period = 0; period < obs.getNrPeriods(); period++ ) {
		for ( unsigned int dm = 0; dm < obs.getNrDMs(); dm++ ) {
			output << obs.getFirstPeriod() + (period * obs.getPeriodStep()) << " " << (obs.getFirstPeriod() + (period * obs.getPeriodStep())) / static_cast< float > (obs.getNrSamplesPerSecond()) << " " << dm << " " << obs.getFirstDM() + (dm * obs.getDMStep()) << " " << fixed << setprecision(3) << snrTable[(period * obs.getNrPaddedDMs()) + dm] << endl;
		}
	}
	output.close();

	try {
		snrTable.deleteDeviceData();
		snrTable.deleteHostData();
	} catch ( OpenCLError err ) {
		cerr << err.what() << endl;
	}

	// Wait for all MPI processes
	world.barrier();

	if ( world.rank() == 0 ) {
		cout << "# processedSeconds nrDMs firstDM DMStep nrPeriods firstPeriod periodStep nrBins nrSamplesPerSecond searchTime searchAverageTime err inputPreProcessingTime inputPreProcessingAverageTime err outputTime H2DTime H2DAverageTime err D2HTime dedispersionTime dedispersionAverageTime err transposeTime transposeAverageTime err foldingTime foldingAverageTime err snrTime" << endl;
		cout << fixed << 1 + obs.getNrSeconds() - secondsToBuffer << " " << obs.getNrDMs() << " " << obs.getFirstDM() << " " << obs.getDMStep() << " " << obs.getNrPeriods() << " " << obs.getFirstPeriod() << " " << obs.getPeriodStep() << " " << obs.getNrBins() << " " << obs.getNrSamplesPerSecond() << " " << setprecision(6) << searchTime.getTotalTime() << " " << searchTime.getAverageTime() << " " << searchTime.getStdDev() << " " << inputPreTime.getTotalTime() << " " << inputPreTime.getAverageTime() << " " << inputPreTime.getStdDev() << " " << outputTime.getTotalTime() << " " << dispersedData.getH2DTimer().getAverageTime() << " " << dispersedData.getH2DTimer().getStdDev() << " " << snrTable.getD2HTimer().getAverageTime() << " " << snrTable.getD2HTimer().getStdDev() << " " << clDedisperse.getTimer().getTotalTime() << " " << clDedisperse.getTimer().getAverageTime() << " " << clDedisperse.getTimer().getStdDev() << " " << clTranspose.getTimer().getTotalTime() << " " << clTranspose.getTimer().getAverageTime() << " " << clTranspose.getTimer().getStdDev() << " " << clFold.getTimer().getTotalTime() << " " << clFold.getTimer().getAverageTime() << " " << clFold.getTimer().getStdDev() << " " << clSNR.getTimer().getTotalTime() << endl;
	}

	return 0;
}
