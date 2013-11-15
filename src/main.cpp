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

#include <configuration.hpp>

#include <ArgumentList.hpp>
using isa::utils::ArgumentList;
#include <utils.hpp>
using isa::utils::giga;
#include <Observation.hpp>
using AstroData::Observation;
#include <ReadData.hpp>
using AstroData::readSIGPROC;
using AstroData::readLOFAR;
#include <Dedispersion.hpp>
using PulsarSearch::Dedispersion;
#include <Shifts.hpp>
using PulsarSearch::getShifts;
#include <Folding.hpp>
using PulsarSearch::Folding;
#include <Transpose.hpp>
using PulsarSearch::Transpose;
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
		readLOFAR(headerFile, dataFile, obs, *input);
	} else if ( dataSIGPROC ) {
		// TODO: implement SIGPROC pipeline
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
	CLData< dataType > dispersedData("DispersedData", true);
	CLData< dataType > dedispersedData("DedispersedData", true);
	CLData< dataType > transposedData("TransposedData", true);
	CLData< dataType > foldedData("FoldedData", true);
	CLData< unsigned int > counterData("CounterData", true);
	CLData< dataType > snrTable("SNRTable", true);

	if ( ((obs.getNrSamplesPerSecond() + (*shifts)[((obs.getNrDMs() - 1) * obs.getNrPaddedChannels())]) % obs.getPadding()) != 0 ) {
		nrSamplesPerChannel = (obs.getNrSamplesPerSecond() + (*shifts)[((obs.getNrDMs() - 1) * obs.getNrPaddedChannels())]) + (obs.getPadding() - ((obs.getNrSamplesPerSecond() + (*shifts)[((obs.getNrDMs() - 1) * obs.getNrPaddedChannels())]) % obs.getPadding()));
	} else {
		nrSamplesPerChannel = (obs.getNrSamplesPerSecond() + (*shifts)[((obs.getNrDMs() - 1) * obs.getNrPaddedChannels())]);
	}
	secondsToBuffer = static_cast< unsigned int >(ceil(static_cast< float >(nrSamplesPerChannel) / obs.getNrSamplesPerPaddedSecond()));

	dispersedData.allocateHostData(secondsToBuffer * obs.getNrChannels() * obs.getNrSamplesPerPaddedSecond());
	snrTable.allocateHostData(obs.getNrPeriods() * obs.getNrPaddedDMs());

	try {
		shifts->setCLContext(clContext);
		shifts->setCLQueue(&((clQueues->at(clDeviceID)).at(0)));
		shifts->allocateDeviceData();
		shifts->copyHostToDevice(true);
		dispersedData.setCLContext(clContext);
		dispersedData.setCLQueue(&((clQueues->at(clDeviceID)).at(0)));
		dispersedData.allocateDeviceData();
		dedispersedData.setCLContext(clContext);
		dedispersedData.setCLQueue(&((clQueues->at(clDeviceID)).at(0)));
		dedispersedData.allocateDeviceData(obs.getNrDMs() * obs.getNrSamplesPerPaddedSecond());
		transposedData.setCLContext(clContext);
		transposedData.setCLQueue(&((clQueues->at(clDeviceID)).at(0)));
		transposedData.allocateDeviceData(obs.getNrSamplesPerSecond() * obs.getNrPaddedDMs());
		foldedData.setCLContext(clContext);
		foldedData.setCLQueue(&((clQueues->at(clDeviceID)).at(0)));
		foldedData.allocateDeviceData(obs.getNrBins() * obs.getNrPeriods() * obs.getNrPaddedDMs());
		foldedData.blankHostData();
		foldedData.allocateDeviceData();
		foldedData.copyHostToDevice();
		foldedData.deleteHostData();
		counterData.setCLContext(clContext);
		counterData.setCLQueue(&((clQueues->at(clDeviceID)).at(0)));
		counterData.allocateDeviceData(obs.getNrBins() * obs.getNrPeriods() * obs.getNrPaddedDMs());
		counterData.blankHostData();
		counterData.allocateDeviceData();
		counterData.copyHostToDevice();
		counterData.deleteHostData();
	} catch ( OpenCLError err ) {
		cerr << err.what() << endl;
		return 1;
	}

	if ( DEBUG ) {
		double hostMemory = 0.0;
		double deviceMemory = 0.0;
		hostMemory += shifts->getHostDataSize() + dispersedData.getHostDataSize() + snrTable.getHostDataSize();
		deviceMemory += shifts->getDeviceDataSize() + dispersedData.getDeviceDataSize() + dedispersedData.getDeviceDataSize() + transposedData.getDeviceDataSize() + foldedData.getDeviceDataSize() + counterData.getDeviceDataSize() + snrTable.getDeviceDataSize();

		cout << "Allocated host memory: " << fixed << setprecision(3) << giga(hostMemory) << endl;
		cout << "Allocated device memory: " << fixed << setprecision(3) << giga(deviceMemory) << endl;
	}

	// Generate OpenCL kernels
	Dedispersion< dataType > clDedisperse("clDedisperse", dataName);
	Transpose< dataType > clTranspose("clTranspose", dataName);
	Folding< dataType > clFold("clFold", dataName);
	SNR< dataType > clSNR("clSNR", dataName);
	try {
		clDedisperse.bindOpenCL(clContext, &(clDevices->at(clDeviceID)), &((clQueues->at(clDeviceID)).at(0)));
		// TODO: tuned parameters
		clDedisperse.setNrSamplesPerBlock();
		clDedisperse.setNrDMsPerBlock();
		clDedisperse.setNrSamplesPerThread();
		clDedisperse.setNrDMsPerThread();
		clDedisperse.setNrSamplesPerDispersedChannel(secondsToBuffer * observation.getNrSamplesPerPaddedSecond());
		clDedisperse.setObservation(&obs);
		clDedisperse.setShifts(shifts);
		clDedisperse.generateCode();
		clTranspose.bindOpenCL(clContext, &(clDevices->at(clDeviceID)), &((clQueues->at(clDeviceID)).at(0)));
		clTranspose.setObservation(&observation);
		// TODO: tuned parameters
		clTranspose.setNrThreadsPerBlock();
		clTranspose.generateCode();
		clFold.bindOpenCL(clContext, &(clDevices->at(clDeviceID)), &((clQueues->at(clDeviceID)).at(0)));
		clFold.setObservation(&observation);
		// TODO: tuned parameters
		clFold.setNrDMsPerBlock();
		clFold.setNrPeriodsPerBlock();
		clFold.setNrBinsPerBlock();
		clFold.setNrDMsPerThread();
		clFold.setNrPeriodsPerThread();
		clFold.setNrBinsPerThread();
		clFold.generateCode();
		clSNR.bindOpenCL(clContext, &(clDevices->at(clDeviceID)), &((clQueues->at(clDeviceID)).at(0)));
		clSNR.setObservation(&observation);
		// TODO: tuned parameters
		clSNR.setNrDMsPerBlock();
		clSNR.setNrPeriodsPerBlock();
		clSNR.setPulsarPipeline();
		clSNR.generateCode();
	} catch ( OpenCLError err ) {
		cerr << err.what() << endl;
		return 1;
	}

	// Search loop
	Timer searchTime("SearchTimer");
	for ( unsigned int second = 0; second <= obs.getNrSeconds() - secondsToBuffer; second++ ) {
		if ( DEBUG ) {
			cout << "Processing second " << second << endl;
		}
		searchTime.start();
		// Prepare the input
		for ( unsigned int channel = 0; channel < obs.getNrChannels(); channel++ ) {
			for ( unsigned int chunk = 0; chunk < secondsToBuffer; chunk++ ) {
				memcpy(dispersedData.getRawHostDataAt((channel * secondsToBuffer * obs.getNrSamplesPerPaddedSecond()) + (chunk * obs.getNrSamplesPerSecond())), (input.at(second + chunk))->getRawHostDataAt(channel * obs.getNrSamplesPerPaddedSecond()), obs.getNrSamplesPerSecond() * sizeof(T));
			}
		}

		// Run the kernels
		try {
			dispersedData.copyHostToDevice();
			clDedisperse(&dispersedData, &dedispersedData);
			clTranspose(&dedispersedData, &transposedData);
			clFold(&transposedData, &foldedData, &counterData);
		} catch ( OpenCLError err ) {
			cerr << err.what() << endl;
			return 1;
		}
		searchTime.stop();
	}

	// Release unnecessary memory
	delete [] input;
	try {
		dispersedData.deleteDeviceData();
		dispersedData.deleteHostData();
		dedispersedData.deleteDeviceData();
		dedispersedData.deleteHostData();
		transposedData.deleteDeviceData();
		transposedData.deleteHostData();
		counterData.deleteDeviceData();
		counterData.deleteHostData();
	} catch ( OpenCLError err ) {
		cerr << err.what() << endl;
	}

	// Store output
	try {
		clSNR(&foldedData, &snrTable);
		snrTable.copyDeviceToHost();
		foldedData.deleteHostData();
		foldedData.deleteDeviceData();
	} catch ( OpenCLError err ) {
		cerr << err.what() << endl;
		return 1;
	}
	ofstream output;
	output.open(outputFile);
	for ( unsigned int period = 0; period < obs.getNrPeriods(); period++ ) {
		for ( unsigned int dm = 0; dm < obs.getNrDMs(); dm++ ) {
			output << period << " " << dm << " " << fixed << setprecision(3) << snrTable[(period * obs.getNrPaddedDMs()) + dm] << endl;
		}
	}
	output.close();

	if ( DEBUG ) {
		cout << "# processedSeconds nrDMs nrPeriods nrBins nrSamplesPerSecond totalTime averageTime err inputAverageTime err outputAverageTime err" << endl;
		cout << 1 + obs.getNrSeconds() - secondsToBuffer << " " << obs.getNrDMs() << " " << obs.getNrPeriods() << " " << obs.getNrBins() << " " << obs.getNrSamplesPerSecond() << " " << searchTime.getTotalTime() << " " << searchTime.getAverageTime() << " " searchTime.getStdDev() << " " << dispersedData.getTimer().getAverageTime() << " " << dispersedData.getTimer().getStdDev() << " " << snrTable.getTimer().getAverageTime() << " " << snrTable.getTimer().getStdDev() << endl;
	}

	return 0;
}
